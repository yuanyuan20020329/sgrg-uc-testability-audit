# -*- coding: utf-8 -*-
"""Part D（单细胞方案）公共模块：h5ad 读取、归一化、基因集打分、pseudobulk"""
import os
import numpy as np
import scipy.sparse as sp
import pandas as pd

RAW_DIR = "data/raw/sc"
GS_DIR = "data/raw/genesets"


# ---------------------------------------------------------------- 读取
def load_h5ad(path, obs_cols=None):
    """返回 (X: csr float32, obs, var, state)
    state='counts'  : X 为原始整数 counts，需自行归一化
    state='lognorm' : X 已是 log 标准化矩阵，直接使用
    CELLxGENE 的 X 有时是 log 标准化、raw/X 才是 counts，这里自动判定。
    """
    import anndata as ad
    A = ad.read_h5ad(path)
    X = A.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.tocsr().astype(np.float32)
    var_src = A.var
    state = "counts" if _is_int(X) else "lognorm"

    if state == "lognorm" and A.raw is not None and A.raw.X is not None:
        RX = A.raw.X
        if not sp.issparse(RX):
            RX = sp.csr_matrix(RX)
        RX = RX.tocsr().astype(np.float32)
        if _is_int(RX):
            X = RX
            var_src = A.raw.var
            state = "counts"

    obs = A.obs.copy()
    var = _mk_var(var_src)
    if obs_cols:
        obs = obs[[c for c in obs_cols if c in obs.columns]]
    return X, obs, var, state


def _mk_var(vdf):
    """统一 var：symbol 优先 gene_symbols -> feature_name -> 索引（Ensembl ID）"""
    for col in ("gene_symbols", "feature_name"):
        if col in vdf.columns:
            sym = vdf[col].astype(str).values
            if not sym[0].startswith("ENSG"):
                return pd.DataFrame({"symbol": sym, "gene_id": vdf.index.astype(str)})
    return pd.DataFrame({"symbol": vdf.index.astype(str), "gene_id": vdf.index.astype(str)})


def _is_int(X, n=2000000):
    d = X.data[:n] if X.data.size > n else X.data
    return bool(np.all(d == np.floor(d)))


def normalize_log1p(X, target=1e4):
    """CP10K + log1p；输入应为原始 counts（CSR）"""
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    D = sp.diags(target / lib)
    X = (D @ X).tocsr()
    X.data = np.log1p(X.data)
    return X


def is_integer_counts(X):
    d = X.data
    return float(np.mean(d == np.floor(d)))


# ---------------------------------------------------------------- 打分
def _mean_expr(Xn):
    return np.asarray(Xn.mean(axis=0)).ravel()


def control_matrix(mean_expr, gene_mask_exclude, n_bins=25, ctrl_size=50, n_perm=200, seed=0):
    """按表达量分箱，为每个箱抽 ctrl_size 个对照基因；返回 V (n_genes, n_perm) 的稀疏指示矩阵（每列归一化）"""
    rng = np.random.default_rng(seed)
    ng = len(mean_expr)
    ok = ~gene_mask_exclude
    idx_ok = np.where(ok)[0]
    v = mean_expr[idx_ok]
    # 分箱（按表达量分位）
    qs = np.quantile(v, np.linspace(0, 1, n_bins + 1))
    qs[0] -= 1e-9
    bins = np.digitize(v, qs) - 1
    bins = np.clip(bins, 0, n_bins - 1)
    bin_idx = {b: idx_ok[bins == b] for b in range(n_bins)}
    bin_idx = {b: a for b, a in bin_idx.items() if len(a) > 0}
    cols, rows = [], []
    for p in range(n_perm):
        chosen = []
        for b, arr in bin_idx.items():
            k = min(ctrl_size, len(arr))
            chosen.append(rng.choice(arr, size=k, replace=False))
        chosen = np.concatenate(chosen)
        cols.append(np.full(len(chosen), p, dtype=np.int32))
        rows.append(chosen.astype(np.int32))
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    V = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=(ng, n_perm))
    # 每列归一化为均值（CSR：按 indices 逐元素除）
    colsum = np.asarray(V.sum(axis=0)).ravel().astype(np.float32)
    V.data = V.data / colsum[V.indices]
    return V.tocsr()


def score_gene_set(Xn, mean_expr, gene_idx, n_bins=25, ctrl_size=50, n_perm=200, seed=0):
    """scanpy 风格 score_genes，返回 (score, zscore, null_mean, null_sd)
    gene_idx 中表达全为 0 的基因会被剔除。"""
    nz = np.asarray(Xn[:, gene_idx].sum(axis=1)).ravel()   # 每细胞在该集上的总表达（判断可用性）
    tot = np.asarray(Xn[:, gene_idx].sum(axis=0)).ravel()
    gene_idx = np.asarray(gene_idx)[tot > 0]
    if len(gene_idx) == 0:
        z = np.zeros(Xn.shape[0], dtype=np.float32)
        return z, z, z, z
    mask = np.zeros(Xn.shape[1], dtype=bool)
    mask[gene_idx] = True
    V = control_matrix(mean_expr, mask, n_bins, ctrl_size, n_perm, seed)
    # 观测集均值
    obs = np.asarray(Xn[:, gene_idx].mean(axis=1)).ravel()
    # 零分布：每个随机集的均值
    null = np.asarray((Xn @ V).todense(), dtype=np.float32)   # (n_cells, n_perm)
    nulmean = null.mean(axis=1)
    nulsd = null.std(axis=1, ddof=1)
    score = obs - nulmean
    z = score / np.where(nulsd > 0, nulsd, np.nan)
    return score.astype(np.float32), z.astype(np.float32), nulmean.astype(np.float32), nulsd.astype(np.float32)


# ---------------------------------------------------------------- 基因集
def load_gmt(path):
    sets = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            sets[p[0]] = [g for g in p[2:] if g]
    return sets


def resolve(genes, sym2idx):
    return np.array(sorted({sym2idx[g] for g in genes if g in sym2idx}), dtype=np.int64)


# ---------------------------------------------------------------- 统计
def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(o[::-1]):
        r = n - rank
        prev = min(prev, p[i] * n / r)
        q[i] = prev
    return q


def perm_test_paired_diff(a, b, n_perm=20000, seed=0):
    """两组独立样本均值差置换检验（双侧）"""
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan, np.nan
    obs = a.mean() - b.mean()
    pool = np.concatenate([a, b])
    n1 = len(a)
    rng = np.random.default_rng(seed)
    cnt = 0
    for _ in range(n_perm):
        rng.shuffle(pool)
        if abs(pool[:n1].mean() - pool[n1:].mean()) >= abs(obs) - 1e-12:
            cnt += 1
    return obs, (cnt + 1) / (n_perm + 1), 1.0 * cnt / n_perm


def cohen_d(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan
    va, vb = a.var(ddof=1), b.var(ddof=1)
    s = np.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    return (a.mean() - b.mean()) / s if s > 0 else np.nan


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)
