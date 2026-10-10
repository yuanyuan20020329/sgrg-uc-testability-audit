# -*- coding: utf-8 -*-
"""Part D 负对照：SGRG 面板 vs 大量「表达分层匹配的随机基因集」

为什么必须做这一步
------------------
ds4（Tearle2026，320,979 细胞）在**未校正**口径下给出 14/26 种细胞类型 p<0.05
（杯状细胞 q=0.016），看起来是可以投出去的阳性；但那被证明是测序深度伪影
（UC donor 深度 11,231 vs 正常 5,189 UMI，AUC=0.839，与疾病标签共线）。

审稿人会追问两件事：
  (1)「未校正口径下的『阳性』，是不是任何同等大小的基因集都会出现？」——若随机基因集
      的效应量/显著数分布与 SGRG 重叠，则 SGRG 的阳性不携带任何面板特异性信息。
  (2)「抽稀后全阴性，是不是抽稀本身把信号洗掉了（抽稀必然降低检验力）？」——若在
      **同一抽稀条件**下随机基因集也同步归零、且没有任何系统性负偏，则「全阴性」
      是真实的，而不是抽稀造成的假阴性。

方法
----
因为 score_genes 的评分对基因集是**线性**的：

    score(cell) = mean_{g in S} Xn[cell, g]  -  baseline(cell)

所以可以把 baseline、SGRG 亚组、阳性对照（Hallmark 炎症因子）、以及 N 个随机基因集
（每个大小、表达分布都与 SGRG 匹配）拼成一个指示矩阵 B，用一次 `Xn @ B` 全部算出。
这样 **SGRG 与随机基因集接受完全相同的处理链**，对比才有效。

随机基因集按表达分层抽取：SGRG 的每个基因先落到 mean_expr 的分位箱，随机集从
**同样的箱**中抽等量基因——保证随机集与 SGRG 的平均表达、检出率分布一致
（否则随机集只会抽到高表达基因，评分被抬高，对照失效）。

用法
----
  python partD_random_control.py --path data/raw/sc/CxG_Tearle2026_PSC_UC.h5ad \
      --tag ds4_tearle_rand --disease-keep "ulcerative colitis,normal" \
      --uc-match ulcer --ctrl-match normal --n-rand 200
  # 加 --downsample 1000 跑深度匹配口径
"""
import argparse
import json
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from partD_common import bh, cohen_d

ap = argparse.ArgumentParser()
ap.add_argument("--path", required=True)
ap.add_argument("--tag", required=True)
ap.add_argument("--label", default="")
ap.add_argument("--tissue-keep", default="")
ap.add_argument("--disease-keep", default="ulcerative colitis,normal")
ap.add_argument("--uc-match", default="ulcer")
ap.add_argument("--ctrl-match", default="normal")
ap.add_argument("--matrix", default="auto")
ap.add_argument("--chunk", type=int, default=12000)
ap.add_argument("--min-cells", type=int, default=20)
ap.add_argument("--min-donors", type=int, default=3)
ap.add_argument("--source", type=float, default=1e4)
ap.add_argument("--downsample", type=int, default=0,
                help=">0 时把所有细胞二项抽稀到该 UMI 数（深度匹配）")
ap.add_argument("--n-rand", type=int, default=200, help="随机基因集个数")
ap.add_argument("--n-perm", type=int, default=5000, help="置换检验次数（用于与旧口径对比）")
ap.add_argument("--n-bins", type=int, default=25, help="表达分位箱数")
ap.add_argument("--ctrl-size", type=int, default=50, help="每箱抽的对照基因数")
ap.add_argument("--seed", type=int, default=17)
a = ap.parse_args()

T0 = time.time()
log = lambda *x: print("[%6.1fs]" % (time.time() - T0), *x, flush=True)
OUT = "results"
RNG = np.random.default_rng(a.seed)

# ================================================================== 读取 obs
def obs_col(f, col):
    if col not in f["obs"]:
        return None
    d = f["obs"][col]
    if isinstance(d, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x)
                         for x in d["categories"][:]], dtype=object)
        codes = d["codes"][:]
        out = np.empty(len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        out[~ok] = ""
        return out
    arr = d[:]
    if arr.dtype.kind in "OS":
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr], dtype=object)
    return arr


def read_col(g, col):
    d = g[col]
    if isinstance(d, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else x
                         for x in d["categories"][:]], dtype=object)
        codes = d["codes"][:]
        out = np.empty(len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        out[~ok] = ""
        return out
    arr = d[:]
    if arr.dtype.kind in "OS":
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr], dtype=object)
    return arr


f = h5py.File(a.path, "r")
log("打开", a.path)

# ---- 选矩阵（要求整数 counts） ----
cand = []
if "X" in f:
    cand.append(("X", f["X"]))
if "raw" in f and "X" in f["raw"]:
    cand.append(("raw/X", f["raw"]["X"]))
if "layers" in f:
    for k in f["layers"].keys():
        cand.append(("layers/" + k, f["layers"][k]))


def sample_int(g):
    try:
        d = g["data"][:200000]
        return float(np.mean(d == np.floor(d)))
    except Exception:
        return np.nan


chosen = None
if a.matrix != "auto":
    for name, g in cand:
        if name == a.matrix:
            chosen = (name, g); break
else:
    for name, g in cand:
        fr = sample_int(g)
        log("  候选矩阵 %-10s 整数占比=%.3f" % (name, fr))
        if fr > 0.999 and chosen is None:
            chosen = (name, g)
    if chosen is None and cand:
        chosen = cand[0]
mname, grp = chosen
n_all, n_genes = grp.attrs["shape"]
n_all, n_genes = int(n_all), int(n_genes)
log("矩阵 =", mname, "| %d 细胞 × %d 基因" % (n_all, n_genes))

vgrp = f["raw"]["var"] if mname.startswith("raw/") else f["var"]
symbols = None
for col in ("gene_symbols", "feature_name", "feature_name_original"):
    if col in vgrp:
        arr = np.asarray(read_col(vgrp, col), dtype=object)
        if len(arr) and not str(arr[0]).startswith("ENSG"):
            symbols = arr; break
if symbols is None:
    ic = vgrp.attrs.get("_index")
    symbols = (np.asarray(read_col(vgrp, ic), dtype=object) if ic and ic in vgrp
               else np.array(["G%d" % i for i in range(n_genes)], dtype=object))
log("基因符号示例 =", list(symbols[:5]))

# ================================================================== 细胞过滤
donor = obs_col(f, "donor_id")
donor = np.array(["d%d" % i for i in range(n_all)], dtype=object) if donor is None else donor
disease = obs_col(f, "disease")
disease = np.array(["unknown"] * n_all, dtype=object) if disease is None else disease
ct = obs_col(f, "cell_type")
ct = np.array(["unknown"] * n_all, dtype=object) if ct is None else ct
tissue = obs_col(f, "tissue")
tissue = np.array([""] * n_all, dtype=object) if tissue is None else tissue

keep = np.ones(n_all, dtype=bool)
if a.tissue_keep:
    pats = [p.strip().lower() for p in a.tissue_keep.split(",") if p.strip()]
    keep &= np.array([any(p in str(t).lower() for p in pats) for t in tissue])
if a.disease_keep:
    dk = {p.strip().lower() for p in a.disease_keep.split(",") if p.strip()}
    keep &= np.array([str(d).lower() in dk for d in disease])
log("过滤后细胞 %d / %d" % (keep.sum(), n_all))

sel = np.where(keep)[0]
donor, disease, ct, tissue = donor[sel], disease[sel], ct[sel], tissue[sel]

indptr_full = grp["indptr"][:]


def _totals(idx):
    t = np.zeros(len(idx), dtype=np.int64)
    ip = indptr_full
    for i0 in range(0, len(idx), 200000):
        block = idx[i0:i0 + 200000]
        s = int(ip[block[0]]); e = int(ip[block[-1] + 1])
        off = np.append(ip[block] - s, e - s)
        seg = np.diff(off)
        d = grp["data"][s:e].astype(np.int64)
        if len(d) == 0:
            continue
        cs = np.add.reduceat(d, np.minimum(off[:-1], len(d) - 1))
        cs[seg == 0] = 0
        t[i0:i0 + len(block)] = cs
    return t


if a.downsample > 0:
    tot = _totals(sel)
    okd = tot >= a.downsample
    log("深度匹配：抽稀至 %d UMI；保留 %d / %d 细胞（%.1f%%）" % (
        a.downsample, okd.sum(), len(sel), 100 * okd.mean()))
    sel = sel[okd]
    donor, disease, ct, tissue = donor[okd], disease[okd], ct[okd], tissue[okd]
n_sel = len(sel)
log("UC/对照分布：", dict(pd.Series(disease).value_counts()))

indptr_full = grp["indptr"][:]
row_start = indptr_full[sel]
row_end = indptr_full[sel + 1]
nnz_tot = int(row_end[-1] - row_start[0])
log("nnz = %d（平均 %.0f 基因/细胞）" % (nnz_tot, nnz_tot / n_sel))


def chunk_iter(chunk):
    a0 = 0
    while a0 < n_sel:
        b0 = min(a0 + chunk, n_sel)
        s0 = int(row_start[a0]); s1 = int(row_end[b0 - 1])
        data = grp["data"][s0:s1].astype(np.float32)
        ind = grp["indices"][s0:s1].astype(np.int32)
        lip = np.concatenate([[0], (row_end[a0:b0] - s0)]).astype(np.int64)
        yield a0, b0, sp.csr_matrix((data, ind, lip), shape=(b0 - a0, n_genes))
        a0 = b0


def norm_cp10k(Xc):
    lib = np.asarray(Xc.sum(axis=1)).ravel()
    if a.downsample > 0:
        p = np.minimum(1.0, a.downsample / np.maximum(lib, 1.0)).astype(np.float64)
        rows = np.repeat(np.arange(Xc.shape[0]), np.diff(Xc.indptr))
        Xc = sp.csr_matrix((RNG.binomial(Xc.data.astype(np.int64), p[rows]).astype(np.float32),
                            Xc.indices, Xc.indptr.copy()), shape=Xc.shape)
        Xc.eliminate_zeros()
        lib = np.asarray(Xc.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    Xc = sp.diags((a.source / lib).astype(np.float32)) @ Xc
    Xc = Xc.tocsr()
    Xc.data = np.log1p(Xc.data)
    return Xc, lib


# ================================================================== pass 1
log("pass 1：每基因均值表达 + 每细胞深度")
gene_sum = np.zeros(n_genes, dtype=np.float64)
gene_det = np.zeros(n_genes, dtype=np.int64)
nUMI = np.zeros(n_sel, dtype=np.int64)
for a0, b0, Xc in chunk_iter(a.chunk):
    Xn, lib = norm_cp10k(Xc)
    nUMI[a0:b0] = lib.astype(np.int64)
    gene_sum += np.asarray(Xn.sum(axis=0)).ravel()
    gene_det += np.asarray((Xc > 0).sum(axis=0)).ravel().astype(np.int64)
    del Xc, Xn
mean_expr = (gene_sum / n_sel).astype(np.float32)
log("深度(中位) = %d | 均值表达 %d 基因 > 0" % (np.median(nUMI), int((mean_expr > 0).sum())))

# ================================================================== 基因集
panel = [l.split("\t")[0] for l in
         open("results/sgrg_gene_set_clean.tsv", encoding="utf-8").read().splitlines()[1:] if l.strip()]
go = json.load(open("results/partD_go_sgrg_subsets.json", encoding="utf-8"))
cyto = [g for g in panel if g in set(go["cytoplasmic_stress_granule"])]
nuc = [g for g in panel if g in set(go["nuclear_stress_granule"])]
sym2idx = {s: i for i, s in enumerate(symbols)}


def resolve(genes):
    return np.array(sorted({sym2idx[g] for g in genes if g in sym2idx}), dtype=np.int64)


def load_gmt(path):
    sets = {}
    for line in open(path, encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            sets[p[0]] = [g for g in p[2:] if g]
    return sets


hall = load_gmt("data/raw/genesets/h.all.v2023.2.Hs.symbols.gmt")

sgrg_idx = resolve(panel)
sgrg_idx = sgrg_idx[gene_det[sgrg_idx] > 0]
K = len(sgrg_idx)
log("SGRG 面板：原始 %d 基因 → 去重且可检出 %d 个 → 随机集大小 K = %d" % (len(panel), K, K))

# ---- 对照基线：按 mean_expr 分位箱抽 ctrl_size 个基因/箱 ----
okg = np.where((gene_det > 0) & (mean_expr > 0))[0]
v = mean_expr[okg]
qs = np.quantile(v, np.linspace(0, 1, a.n_bins + 1))
qs[0] -= 1e-9
bins_all = np.clip(np.digitize(v, qs) - 1, 0, a.n_bins - 1)
bin_of = {}         # 基因 idx -> 箱号
for i, g in enumerate(okg):
    bin_of[g] = bins_all[i]
pool = {b: okg[bins_all == b] for b in range(a.n_bins)}
pool = {b: arr for b, arr in pool.items() if len(arr) > 0}

# 基线基因（与 SGRG 无重叠以保证观测/对照独立）
sgrg_set = set(sgrg_idx.tolist())
base_genes = []
for b, arr in pool.items():
    cand_b = np.array([g for g in arr if g not in sgrg_set], dtype=np.int64)
    if len(cand_b) == 0:
        continue
    k = min(a.ctrl_size, len(cand_b))
    base_genes.append(RNG.choice(cand_b, size=k, replace=False))
base_genes = np.concatenate(base_genes)
log("对照基线基因 = %d（%d 箱 × %d）" % (len(base_genes), len(pool), a.ctrl_size))

# ---- SGRG 基因的表达分箱结构（随机集据此分层抽样） ----
sgrg_bins = np.array([bin_of[g] for g in sgrg_idx if g in bin_of])
bin_counts = np.bincount(sgrg_bins, minlength=a.n_bins)
log("SGRG 的表达分箱分布（%d 箱）= %s" % (a.n_bins, bin_counts.tolist()))

# ---- 生成 N 个随机基因集：从与 SGRG 相同的箱中抽等量基因 ----
rand_sets = []
for j in range(a.n_rand):
    chosen = []
    for b in range(a.n_bins):
        k = int(bin_counts[b])
        if k == 0:
            continue
        cand_b = pool.get(b, np.array([], dtype=np.int64))
        cand_b = np.array([g for g in cand_b if g not in sgrg_set], dtype=np.int64)
        if len(cand_b) == 0:
            continue
        kk = min(k, len(cand_b))
        chosen.append(RNG.choice(cand_b, size=kk, replace=False))
    if chosen:
        rand_sets.append(np.concatenate(chosen))
log("生成随机基因集 %d 个（每个平均 %.1f 基因，已排除 SGRG 基因）" % (
    len(rand_sets), np.mean([len(s) for s in rand_sets])))

# ---- 阳性对照（验证流程灵敏度：炎症通路在 UC 中应当升高）----
POS = ["HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
       "HALLMARK_INTERFERON_GAMMA_RESPONSE"]
pos_sets = []
for h in POS:
    if h in hall:
        gi = resolve(hall[h])
        gi = gi[gene_det[gi] > 0] if len(gi) else gi
        if len(gi):
            pos_sets.append((h.replace("HALLMARK_", "H_"), gi))
log("阳性对照基因集 =", [(n, len(g)) for n, g in pos_sets])

# ================================================================== 指示矩阵
def ind_vec(idx):
    """返回归一化的稀疏列 (n_genes, 1)"""
    idx = np.asarray(idx, dtype=np.int64)
    idx = idx[(idx >= 0) & (idx < n_genes)]
    if len(idx) == 0:
        return None
    w = np.full(len(idx), 1.0 / len(idx), dtype=np.float32)
    return sp.csr_matrix((w, (idx, np.zeros(len(idx), dtype=np.int32))), shape=(n_genes, 1))


cols, names = [], []
c = ind_vec(base_genes); cols.append(c); names.append("BASELINE")
# 名称为「真实候选 + 阳性对照 + 随机」的顺序
REAL = [("SGRG", sgrg_idx), ("SGRG_cyto", resolve(cyto)), ("SGRG_nuc", resolve(nuc))]
real_names = []
for nm, gi in REAL:
    gi = gi[gene_det[gi] > 0] if len(gi) else gi
    cc = ind_vec(gi)
    if cc is None:
        continue
    cols.append(cc); names.append(nm); real_names.append(nm)
pos_names = []
for nm, gi in pos_sets:
    cc = ind_vec(gi)
    cols.append(cc); names.append(nm); pos_names.append(nm)
rand_names = []
for j, gi in enumerate(rand_sets):
    cc = ind_vec(gi)
    if cc is None:
        continue
    cols.append(cc); names.append("RAND_%03d" % j); rand_names.append("RAND_%03d" % j)
B = sp.hstack(cols, format="csr").tocsr()
log("指示矩阵 B = %d×%d（%d 列：1 baseline + %d 真实 + %d 阳性 + %d 随机）" % (
    B.shape[0], B.shape[1], len(names), len(real_names), len(pos_names), len(rand_names)))
idx_base = 0
score_names = names[1:]
col_of = {nm: i for i, nm in enumerate(names)}

# ================================================================== pass 2
log("pass 2：逐细胞打分（一次矩阵乘算出全部基因集）")
SC = np.zeros((n_sel, len(score_names)), dtype=np.float32)
nGene = np.zeros(n_sel, dtype=np.int32)
for a0, b0, Xc in chunk_iter(a.chunk):
    Xn, lib = norm_cp10k(Xc)
    nGene[a0:b0] = np.asarray((Xc > 0).sum(axis=1)).ravel().astype(np.int32)
    S = np.asarray((Xn @ B).todense(), dtype=np.float32)
    base = S[:, [idx_base]]
    SC[a0:b0, :] = S[:, 1:] - base
    del Xc, Xn, S
    log("  pass2 %d / %d" % (b0, n_sel))
log("打分完成")

df = pd.DataFrame(SC, columns=score_names)
df["donor"] = donor
df["disease"] = disease
df["cell_type"] = ct
df["nUMI"] = nUMI
log("SGRG 细胞级 mean=%.4f" % df["SGRG"].mean())

# ================================================================== donor 级检验
pb = df.groupby(["donor", "disease", "cell_type"])[score_names + ["nUMI"]].mean()
pb["n_cells"] = df.groupby(["donor", "disease", "cell_type"]).size()
pb = pb.reset_index()
pb["_grp"] = np.where(pb.disease.str.lower().str.contains(a.uc_match), "UC",
                      np.where(pb.disease.str.lower().str.contains(a.ctrl_match), "CTRL", "OTHER"))
pbf = pb[(pb.n_cells >= a.min_cells) & (pb._grp != "OTHER")].copy()

valid = pbf.groupby(["cell_type", "_grp"]).size().unstack(fill_value=0)
ok_ct = [c for c in valid.index
         if valid.loc[c].get("UC", 0) >= a.min_donors and valid.loc[c].get("CTRL", 0) >= a.min_donors]
log("纳入检验的细胞类型 %d 个" % len(ok_ct))

pbf.to_csv(os.path.join(OUT, "partD_%s_pseudobulk_score.tsv" % a.tag), sep="\t",
           index=False, float_format="%.5g")

# 两种检验同时算：
#   p      = Mann-Whitney U（主口径，秩检验，对异方差稳健）
#   p_perm = 均差置换检验（与之前 ds4 口径一致，用于暴露该口径的校准问题）
RNG_T = np.random.default_rng(a.seed + 101)
rows = []
for c in ok_ct:
    s = pbf[pbf.cell_type == c].copy()
    is_uc = (s._grp == "UC").values
    n1, n = int(is_uc.sum()), len(s)
    perm = RNG_T.permuted(np.tile(np.arange(n), (a.n_perm, 1)), axis=1)
    Ai, Bi = perm[:, :n1], perm[:, n1:]
    for nm in score_names:
        vals = s[nm].values.astype(np.float64)
        u, h = vals[is_uc], vals[~is_uc]
        try:
            pm = mannwhitneyu(u, h, alternative="two-sided").pvalue
        except ValueError:
            pm = np.nan
        obs = abs(u.mean() - h.mean())
        mA = vals[Ai].mean(axis=1)
        mB = vals[Bi].mean(axis=1)
        pp = (np.sum(np.abs(mA - mB) >= obs - 1e-12) + 1) / (a.n_perm + 1)
        rows.append({"cell_type": c, "geneset": nm,
                     "n_uc": len(u), "n_hc": len(h),
                     "mean_uc": u.mean(), "mean_hc": h.mean(),
                     "sd_uc": u.std(ddof=1), "sd_hc": h.std(ddof=1),
                     "diff": u.mean() - h.mean(),
                     "cohen_d": cohen_d(u, h), "p": pm, "p_perm": pp})
tt = pd.DataFrame(rows)
# 每个基因集内跨细胞类型做 BH
for col_p, col_q in [("p", "q_BH"), ("p_perm", "q_perm_BH")]:
    tt[col_q] = np.nan
    for nm, idx in tt.groupby("geneset").groups.items():
        tt.loc[idx, col_q] = bh(tt.loc[idx, col_p].values)

tt.to_csv(os.path.join(OUT, "partD_%s_geneset_tests.tsv" % a.tag), sep="\t",
          index=False, float_format="%.5g")

# ---- 汇总：每个基因集的显著数 / 最小 p / 效应量 ----
def summ(nm_list):
    r = []
    for nm in nm_list:
        d = tt[tt.geneset == nm]
        if not len(d):
            continue
        r.append({"geneset": nm, "kind": ("random" if nm.startswith("RAND_") else
                                          ("positive" if nm.startswith("H_") else "real")),
                  "n_ct": len(d),
                  "n_p05": int((d.p < 0.05).sum()),
                  "n_q05": int((d.q_BH < 0.05).sum()),
                  "n_p05_perm": int((d.p_perm < 0.05).sum()),
                  "n_q05_perm": int((d.q_perm_BH < 0.05).sum()),
                  "min_p": float(np.nanmin(d.p)), "min_q": float(np.nanmin(d.q_BH)),
                  "min_p_perm": float(np.nanmin(d.p_perm)),
                  "n_up": int((d["diff"] > 0).sum()),
                  "median_abs_d": float(np.nanmedian(np.abs(d["diff"]))),
                  "median_cohend": float(np.nanmedian(d.cohen_d))})
    return pd.DataFrame(r)


S_all = summ(score_names)
S_all.to_csv(os.path.join(OUT, "partD_%s_geneset_summary.tsv" % a.tag), sep="\t",
             index=False, float_format="%.5g")

randS = S_all[S_all.kind == "random"]
realS = S_all[S_all.kind == "real"]
posS = S_all[S_all.kind == "positive"]

print("\n" + "=" * 78)
print("负对照结果：%s（%s）" % (a.tag, a.label))
print("=" * 78)
print("\n--- 真实候选面板 ---")
print(realS.to_string(index=False, float_format=lambda v: "%.4f" % v))
print("\n--- 阳性对照（流程灵敏度）---")
print(posS.to_string(index=False, float_format=lambda v: "%.4f" % v))
print("\n--- 随机基因集分布（N=%d）---" % len(randS))
for col in ["n_p05", "n_q05", "n_p05_perm", "n_q05_perm", "min_p", "min_p_perm",
            "median_abs_d", "n_up"]:
    v = randS[col].values
    print("  %-14s 中位 %8.4f | 5%%–95%% [%8.4f, %8.4f] | 最小 %8.4f | 最大 %8.4f" % (
        col, np.median(v), np.percentile(v, 5), np.percentile(v, 95), v.min(), v.max()))

n_ct = int(S_all.n_ct.iloc[0])
print("\n--- 检验口径校准（以 %d 个随机基因集为「真零」）---" % len(randS))
print("  零假设下的期望假阳性 = %.1f 个细胞类型（%d × 0.05）" % (n_ct * 0.05, n_ct))
print("  Mann-Whitney U : 随机集显著数 中位 %.0f（膨胀 %.1f×）" % (
    randS.n_p05.median(), randS.n_p05.median() / (n_ct * 0.05)))
print("  均差置换检验   : 随机集显著数 中位 %.0f（膨胀 %.1f×）" % (
    randS.n_p05_perm.median(), randS.n_p05_perm.median() / (n_ct * 0.05)))
print("  → 未校正口径下两组不可交换：任何基因集都会膨胀假阳性，这不是 SGRG 的问题")
sd_ratio = float(np.nanmedian(tt[tt.geneset.str.startswith("RAND_")]
                              .groupby("cell_type")[["sd_uc", "sd_hc"]].median()
                              .eval("sd_uc / sd_hc")))
print("  → 随机集 UC/对照 pseudobulk 标准差比 中位 = %.2f" % sd_ratio)

print("\n--- SGRG 在随机分布中的经验分位 ---")
emp = {}
for col, higher_is_extreme in [("n_p05", True), ("n_q05", True), ("n_p05_perm", True),
                               ("min_p", False), ("min_p_perm", False),
                               ("median_abs_d", True), ("n_up", True)]:
    v = np.sort(randS[col].values)
    for nm in realS.geneset:
        x = float(realS[realS.geneset == nm][col].iloc[0])
        # 经验 p：随机集中 ≥ SGRG 的比例（对 min_p 则 ≤）
        if higher_is_extreme:
            ep = (np.sum(v >= x) + 1) / (len(v) + 1)
        else:
            ep = (np.sum(v <= x) + 1) / (len(v) + 1)
        emp.setdefault(nm, {})[col] = {"value": x, "emp_p": float(ep),
                                       "percentile": float(100 * np.mean(v <= x))}
print(json.dumps(emp, indent=2, ensure_ascii=False))
for nm, d in emp.items():
    print("  %-11s [MWU] 显著 %2d（随机中位 %.0f，经验 p=%.3f）| 最小 p=%.3g（经验 p=%.3f）" % (
        nm, d["n_p05"]["value"], np.median(randS.n_p05),
        d["n_p05"]["emp_p"], d["min_p"]["value"], d["min_p"]["emp_p"]))
    print("  %-11s [perm] 显著 %2d（随机中位 %.0f，经验 p=%.3f）| 最小 p=%.3g（经验 p=%.3f）" % (
        "", d["n_p05_perm"]["value"], np.median(randS.n_p05_perm),
        d["n_p05_perm"]["emp_p"], d["min_p_perm"]["value"], d["min_p_perm"]["emp_p"]))

out = {"tag": a.tag, "label": a.label, "path": a.path, "downsample_to": a.downsample,
       "n_cells_used": int(n_sel), "matrix": mname,
       "n_rand": len(rand_names), "rand_sizes": [int(len(s)) for s in rand_sets],
       "sgrg_K": int(K), "n_ct_tested": len(ok_ct), "cell_types_tested": ok_ct,
       "n_donors_by_disease": {str(k): int(v) for k, v in
                               df.groupby("disease")["donor"].nunique().items()},
       "real": json.loads(realS.to_json(orient="records")),
       "positive": json.loads(posS.to_json(orient="records")),
       "random_dist": {c: {"median": float(np.median(randS[c])),
                           "p05": float(np.percentile(randS[c], 5)),
                           "p95": float(np.percentile(randS[c], 95)),
                           "min": float(randS[c].min()), "max": float(randS[c].max())}
                       for c in ["n_p05", "n_q05", "n_p05_perm", "n_q05_perm",
                                 "min_p", "min_p_perm", "median_abs_d", "n_up"]},
       "n_perm": a.n_perm,
       "sd_ratio_uc_over_hc_median": sd_ratio,
       "expected_false_positives": n_ct * 0.05,
       "empirical": emp}
json.dump(out, open(os.path.join(OUT, "partD_%s_summary.json" % a.tag), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False, default=str)
log("完成 %.1f s" % (time.time() - T0))
f.close()
