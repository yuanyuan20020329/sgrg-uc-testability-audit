# -*- coding: utf-8 -*-
"""Part D / 按「数据来源」重做 ds4 的疾病对比（含同 donor 炎症配对）

背景（本轮审计结论）
--------------------
ds4（Tearle2026, CELLxGENE 5bd4d642）里 320,979 个细胞中：
  * 只有 10,519 个 normal 细胞是**本研究自己**的对照（HC-1 / HC-2，is_primary_data=True）；
  * 另外 177,982 个 normal 细胞（94.4%）是 **is_primary_data=False 的复用数据**
    （411c/414c/417c/432c/440c，批次 A3x-*/CRC-*，切除标本，10x 5' v2，中位深度低到 289–9769）；
  * 病例侧（UC 62,968 / PSC 69,510）全部 is_primary_data=True。

后果：上一轮把「正常」当成一个组来比，得到的
  UC 中位深度 5674 vs 正常 1305（4.3×）
其实是**对照组成**造成的，而不是疾病造成的。只取本研究自己的对照：
  UC 5674 vs HC 5440（1.04×）——深度差消失。

所以本脚本在**同一遍扫描**里同时算两套口径，直接对比：
  A. naive  : UC vs 全部 normal（= 上一轮口径，含 94.4% 外部复用数据）
  B. primary: UC vs 本研究自己的 normal（is_primary_data=True，n=2）
并对 B 之外再补一个**设计上真正受控**的对比：
  C. within-donor: 同一 donor 内「炎症黏膜 vs 非炎症黏膜」（配对，去除个体与批次差异）

三者都用同一套「表达分层匹配的随机基因集」做零分布校准。

用法
----
  python partD_provenance_analysis.py --path data/raw/sc/CxG_Tearle2026_PSC_UC.h5ad \
      --tag ds4_prov --n-rand 200 [--downsample 1000]
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
from scipy.stats import mannwhitneyu, wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from partD_common import bh, cohen_d

ap = argparse.ArgumentParser()
ap.add_argument("--path", required=True)
ap.add_argument("--tag", required=True)
ap.add_argument("--label", default="")
ap.add_argument("--disease-keep", default="ulcerative colitis,normal,primary sclerosing cholangitis")
ap.add_argument("--uc-match", default="ulcer")
ap.add_argument("--ctrl-match", default="normal")
ap.add_argument("--chunk", type=int, default=12000)
ap.add_argument("--min-cells", type=int, default=20)
ap.add_argument("--min-donors", type=int, default=2)
ap.add_argument("--source", type=float, default=1e4)
ap.add_argument("--downsample", type=int, default=0)
ap.add_argument("--n-rand", type=int, default=200)
ap.add_argument("--n-perm", type=int, default=2000)
ap.add_argument("--n-bins", type=int, default=25)
ap.add_argument("--ctrl-size", type=int, default=50)
ap.add_argument("--seed", type=int, default=17)
a = ap.parse_args()

T0 = time.time()
log = lambda *x: print("[%6.1fs]" % (time.time() - T0), *x, flush=True)
OUT = "results"
RNG = np.random.default_rng(a.seed)
RNG_DS = np.random.default_rng(20261009)      # 抽稀专用（与集抽样解耦，保证可复现）


# ------------------------------------------------------------------ obs 读取
def obs_col(g, col):
    """从任意 h5ad DataFrame 组（obs / var / raw.var）读一列，兼容 Dataset 与分类 Group"""
    if col not in g:
        return None
    d = g[col]
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


f = h5py.File(a.path, "r")
log("打开", a.path)

cand = []
if "X" in f:
    cand.append(("X", f["X"]))
if "raw" in f and "X" in f["raw"]:
    cand.append(("raw/X", f["raw"]["X"]))
if "layers" in f:
    for k in f["layers"].keys():
        cand.append(("layers/" + k, f["layers"][k]))
chosen = None
for nm, g in cand:
    try:
        d = g["data"][:200000]
        fr = float(np.mean(d == np.floor(d)))
        log("  候选矩阵 %-10s 整数占比=%.3f" % (nm, fr))
        if fr > 0.999 and chosen is None:
            chosen = (nm, g)
    except Exception:
        pass
if chosen is None:
    chosen = cand[0]
mname, grp = chosen
n_all, n_genes = grp.attrs["shape"]
n_all, n_genes = int(n_all), int(n_genes)
log("矩阵 =", mname, "|", n_all, "细胞 ×", n_genes, "基因")

vgrp = f["raw"]["var"] if mname.startswith("raw/") else f["var"]
symbols = None
for col in ("gene_symbols", "feature_name", "feature_name_original"):
    if col in vgrp:
        arr = np.asarray(obs_col(vgrp, col), dtype=object)
        if len(arr) and not str(arr[0]).startswith("ENSG"):
            symbols = arr; break
if symbols is None:
    ic = vgrp.attrs.get("_index")
    symbols = (np.asarray(obs_col(vgrp, ic), dtype=object) if ic and ic in vgrp
               else np.array(["G%d" % i for i in range(n_genes)], dtype=object))
log("基因符号示例 =", list(symbols[:5]))

# ------------------------------------------------------------------ 细胞过滤
donor = obs_col(f["obs"], "donor_id")
donor = np.array(["d%d" % i for i in range(n_all)], dtype=object) if donor is None else donor
disease = obs_col(f["obs"], "disease")
disease = np.array(["unknown"] * n_all, dtype=object) if disease is None else disease
ct = obs_col(f["obs"], "cell_type")
ct = np.array(["unknown"] * n_all, dtype=object) if ct is None else ct
prim = obs_col(f["obs"], "is_primary_data")
prim = np.zeros(n_all, dtype=bool) if prim is None else np.asarray(prim).astype(bool)
infl = obs_col(f["obs"], "inflammation")
infl = np.array([""] * n_all, dtype=object) if infl is None else infl

dk = {p.strip().lower() for p in a.disease_keep.split(",") if p.strip()}
keep = np.array([str(d).lower() in dk for d in disease])
log("过滤后细胞 %d / %d" % (keep.sum(), n_all))

sel = np.where(keep)[0]
donor, disease, ct = donor[sel], disease[sel], ct[sel]
prim, infl = prim[sel], infl[sel]

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
    log("深度匹配：抽稀至 %d UMI；保留 %d / %d（%.1f%%）" % (
        a.downsample, okd.sum(), len(sel), 100 * okd.mean()))
    sel = sel[okd]
    donor, disease, ct = donor[okd], disease[okd], ct[okd]
    prim, infl = prim[okd], infl[okd]
n_sel = len(sel)

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
        Xc = sp.csr_matrix((RNG_DS.binomial(Xc.data.astype(np.int64), p[rows]).astype(np.float32),
                            Xc.indices, Xc.indptr.copy()), shape=Xc.shape)
        Xc.eliminate_zeros()
        lib = np.asarray(Xc.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    Xc = sp.diags((a.source / lib).astype(np.float32)) @ Xc
    Xc = Xc.tocsr()
    Xc.data = np.log1p(Xc.data)
    return Xc, lib


# ------------------------------------------------------------------ pass 1
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
log("深度中位 = %d" % np.median(nUMI))

# ------------------------------------------------------------------ 基因集
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
log("SGRG 面板：原始 %d → 可检出去重 %d（随机集大小 K=%d）" % (len(panel), K, K))

okg = np.where((gene_det > 0) & (mean_expr > 0))[0]
v = mean_expr[okg]
qs = np.quantile(v, np.linspace(0, 1, a.n_bins + 1))
qs[0] -= 1e-9
bins_all = np.clip(np.digitize(v, qs) - 1, 0, a.n_bins - 1)
bin_of = {int(g): int(b) for g, b in zip(okg, bins_all)}
pool = {b: okg[bins_all == b] for b in range(a.n_bins)}
pool = {b: arr for b, arr in pool.items() if len(arr) > 0}

sgrg_set = set(sgrg_idx.tolist())
base_genes = []
for b, arr in pool.items():
    cb = np.array([g for g in arr if g not in sgrg_set], dtype=np.int64)
    if len(cb) == 0:
        continue
    base_genes.append(RNG.choice(cb, size=min(a.ctrl_size, len(cb)), replace=False))
base_genes = np.concatenate(base_genes)
log("对照基线基因 = %d" % len(base_genes))

sgrg_bins = np.array([bin_of[g] for g in sgrg_idx if g in bin_of])
bin_counts = np.bincount(sgrg_bins, minlength=a.n_bins)

rand_sets = []
for j in range(a.n_rand):
    chosen = []
    for b in range(a.n_bins):
        k = int(bin_counts[b])
        if k == 0:
            continue
        cb = np.array([g for g in pool.get(b, np.array([], dtype=np.int64)) if g not in sgrg_set],
                      dtype=np.int64)
        if len(cb) == 0:
            continue
        chosen.append(RNG.choice(cb, size=min(k, len(cb)), replace=False))
    if chosen:
        rand_sets.append(np.concatenate(chosen))
log("随机基因集 %d 个（平均 %.1f 基因）" % (len(rand_sets), np.mean([len(s) for s in rand_sets])))

POS = ["HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
       "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_HYPOXIA",
       "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]
pos_sets = [(h.replace("HALLMARK_", "H_"), resolve(hall[h])) for h in POS if h in hall]


def ind_vec(idx):
    idx = np.asarray(idx, dtype=np.int64)
    idx = idx[(idx >= 0) & (idx < n_genes)]
    if len(idx) == 0:
        return None
    w = np.full(len(idx), 1.0 / len(idx), dtype=np.float32)
    return sp.csr_matrix((w, (idx, np.zeros(len(idx), dtype=np.int32))), shape=(n_genes, 1))


cols, names = [ind_vec(base_genes)], ["BASELINE"]
real_names = []
for nm, gi in [("SGRG", sgrg_idx), ("SGRG_cyto", resolve(cyto)), ("SGRG_nuc", resolve(nuc))]:
    gi = gi[gene_det[gi] > 0] if len(gi) else gi
    c = ind_vec(gi)
    if c is None:
        continue
    cols.append(c); names.append(nm); real_names.append(nm)
pos_names = []
for nm, gi in pos_sets:
    gi = gi[gene_det[gi] > 0]
    c = ind_vec(gi)
    if c is None:
        continue
    cols.append(c); names.append(nm); pos_names.append(nm)
rand_names = []
for j, gi in enumerate(rand_sets):
    c = ind_vec(gi)
    if c is None:
        continue
    cols.append(c); names.append("RAND_%03d" % j); rand_names.append("RAND_%03d" % j)
B = sp.hstack(cols, format="csr").tocsr()
log("指示矩阵 B = %d×%d（1 baseline + %d 真实 + %d 阳性 + %d 随机）" % (
    B.shape[0], B.shape[1], len(real_names), len(pos_names), len(rand_names)))
score_names = names[1:]

# ------------------------------------------------------------------ pass 2
log("pass 2：逐细胞打分")
SC = np.zeros((n_sel, len(score_names)), dtype=np.float32)
nGene = np.zeros(n_sel, dtype=np.int32)
for a0, b0, Xc in chunk_iter(a.chunk):
    Xn, lib = norm_cp10k(Xc)
    nGene[a0:b0] = np.asarray((Xc > 0).sum(axis=1)).ravel().astype(np.int32)
    S = np.asarray((Xn @ B).todense(), dtype=np.float32)
    SC[a0:b0, :] = S[:, 1:] - S[:, [0]]
    del Xc, Xn, S
log("打分完成")

df = pd.DataFrame(SC, columns=score_names)
df["donor"] = donor
df["disease"] = disease
df["cell_type"] = ct
df["primary"] = prim
df["inflammation"] = infl
df["nUMI"] = nUMI
del SC

df["_grp"] = np.where(df.disease.str.lower().str.contains(a.uc_match), "UC",
              np.where(df.disease.str.lower().str.contains(a.ctrl_match), "CTRL", "OTHER"))
df.to_pickle(os.path.join(OUT, "partD_%s_cells.pkl" % a.tag))
log("per-cell 表已存（%d × %d）" % df.shape)

# ================================================================== donor 级 pseudobulk
pb = df.groupby(["donor", "disease", "_grp", "primary", "cell_type"])[score_names + ["nUMI"]].mean()
pb["n_cells"] = df.groupby(["donor", "disease", "_grp", "primary", "cell_type"]).size()
pb = pb.reset_index()
pb.to_csv(os.path.join(OUT, "partD_%s_pseudobulk.tsv" % a.tag), sep="\t",
          index=False, float_format="%.5g")


def run_contrast(sub, tagA, groupcol, min_donors, min_cells, n_perm, order=("UC", "CTRL")):
    """sub: donor×celltype×组别 的 pseudobulk 表。
    显式固定组顺序 order=(病例, 对照)，保证 diff = 病例 − 对照 的符号稳定。"""
    sub = sub[sub.n_cells >= min_cells]
    valid = sub.groupby(["cell_type", groupcol]).size().unstack(fill_value=0)
    need = [g for g in order if g in set(sub[groupcol].unique())]
    ok_ct = [c for c in valid.index if all(valid.loc[c].get(g, 0) >= min_donors for g in need)]
    rows = []
    RT = np.random.default_rng(a.seed + 313)
    for c in ok_ct:
        s = sub[sub.cell_type == c]
        g0 = s[s[groupcol] == need[0]]      # 病例
        g1 = s[s[groupcol] == need[1]]      # 对照
        vals = np.concatenate([g0[score_names].values, g1[score_names].values])
        n0, n = len(g0), len(s)
        perm = RT.permuted(np.tile(np.arange(n), (n_perm, 1)), axis=1)
        Ai, Bi = perm[:, :n0], perm[:, n0:]
        for k, nm in enumerate(score_names):
            x, y = g0[nm].values.astype(np.float64), g1[nm].values.astype(np.float64)
            try:
                pm = mannwhitneyu(x, y, alternative="two-sided").pvalue
            except ValueError:
                pm = np.nan
            obs = abs(x.mean() - y.mean())
            mA = vals[Ai, k].mean(axis=1); mB = vals[Bi, k].mean(axis=1)
            pp = (np.sum(np.abs(mA - mB) >= obs - 1e-12) + 1) / (n_perm + 1)
            rows.append({"contrast": tagA, "cell_type": c, "geneset": nm,
                         "n_a": len(x), "n_b": len(y),
                         "mean_from": need[0], "mean_against": need[1],
                         "mean_a": x.mean(), "mean_b": y.mean(),
                         "sd_a": x.std(ddof=1), "sd_b": y.std(ddof=1),
                         "diff": x.mean() - y.mean(),
                         "cohen_d": cohen_d(x, y), "p": pm, "p_perm": pp})
    t = pd.DataFrame(rows)
    for cp, cq in [("p", "q_BH"), ("p_perm", "q_perm_BH")]:
        t[cq] = np.nan
        for nm, idx in t.groupby("geneset").groups.items():
            t.loc[idx, cq] = bh(t.loc[idx, cp].values)
    return t, ok_ct


# ---- A. naive：UC vs 全部 normal（复现上一轮口径）
pbA = pb[(pb._grp.isin(["UC", "CTRL"]))]
tA, ctA = run_contrast(pbA, "A_naive_all_normal", "_grp", 3, a.min_cells, a.n_perm)
# ---- B. primary：UC vs 本研究自己的 normal
pbB = pb[(pb._grp.isin(["UC", "CTRL"])) & ((pb._grp == "UC") | pb.primary)]
tB, ctB = run_contrast(pbB, "B_primary_only", "_grp", a.min_donors, a.min_cells, a.n_perm)
T = pd.concat([tA, tB], ignore_index=True)
T.to_csv(os.path.join(OUT, "partD_%s_contrast_tests.tsv" % a.tag), sep="\t",
         index=False, float_format="%.5g")
log("A 口径纳入 %d 细胞类型 | B 口径纳入 %d 细胞类型" % (len(ctA), len(ctB)))

# ---- C. within-donor：炎症 vs 非炎症（同 donor 配对）
rows = []
RT = np.random.default_rng(a.seed + 777)
okc = []
for c, s in df[df.inflammation != ""].groupby("cell_type"):
    piv = s.pivot_table(index="donor", columns="inflammation", values="SGRG", aggfunc="mean")
    cnt = s.pivot_table(index="donor", columns="inflammation", values="SGRG", aggfunc="size")
    if not {"inflamed", "noninflamed"}.issubset(piv.columns):
        continue
    ok = (cnt[["inflamed", "noninflamed"]] >= a.min_cells).all(axis=1)
    if ok.sum() < 3:
        continue
    okc.append((c, int(ok.sum())))
    sub = s[s.donor.isin(piv.index[ok])]
    din = sub[sub.inflammation == "inflamed"].groupby("donor")[score_names].mean()
    dno = sub[sub.inflammation == "noninflamed"].groupby("donor")[score_names].mean()
    j = din.join(dno, lsuffix="_infl", rsuffix="_non")
    nd = len(j)
    for nm in score_names:
        d = (j[nm + "_infl"] - j[nm + "_non"]).values.astype(np.float64)
        try:
            pw = wilcoxon(d).pvalue
        except ValueError:
            pw = np.nan
        # 符号翻转置换
        sign = RT.choice([-1, 1], size=(2000, nd))
        null = np.abs((sign * d).mean(axis=1))
        pp = (np.sum(null >= abs(d.mean()) - 1e-12) + 1) / 2001
        rows.append({"contrast": "C_within_donor_inflamed", "cell_type": c, "geneset": nm,
                     "n_donors": nd, "mean_diff": d.mean(),
                     "mean_inflamed": j[nm + "_infl"].mean(), "mean_noninflamed": j[nm + "_non"].mean(),
                     "p_wilcoxon": pw, "p_signflip": pp,
                     "n_pos": int((d > 0).sum())})
C = pd.DataFrame(rows)
if len(C):
    C["q_BH"] = np.nan
    for nm, idx in C.groupby("geneset").groups.items():
        C.loc[idx, "q_BH"] = bh(C.loc[idx, "p_signflip"].values)
    C.to_csv(os.path.join(OUT, "partD_%s_inflamed_tests.tsv" % a.tag), sep="\t",
             index=False, float_format="%.5g")
log("C 口径（同 donor 炎症配对）纳入细胞类型：%s" % okc)

# ================================================================== 汇总 & 校准
def summ(tbl, kindcol="geneset"):
    r = []
    for nm, d in tbl.groupby(kindcol):
        r.append({"geneset": nm,
                  "kind": ("random" if nm.startswith("RAND_") else
                           ("positive" if nm.startswith("H_") else "real")),
                  "n_ct": len(d),
                  "n_p05": int((d.p < 0.05).sum()),
                  "n_q05": int((d.q_BH < 0.05).sum()),
                  "n_p05_perm": int((d.p_perm < 0.05).sum()) if "p_perm" in d else np.nan,
                  "min_p": float(np.nanmin(d.p)) if "p" in d else np.nan,
                  "n_up": int((d["diff"] > 0).sum()) if "diff" in d else np.nan,
                  "median_abs_d": float(np.nanmedian(np.abs(d["diff"]))) if "diff" in d else np.nan})
    return pd.DataFrame(r)


Ssum = []
for tagname, tbl in [("A_naive", tA), ("B_primary", tB)]:
    s = summ(tbl); s["contrast"] = tagname; Ssum.append(s)
if len(C):
    cc = C.copy(); cc["p"] = cc.p_wilcoxon
    cc["q_BH"] = cc.q_BH; cc["diff"] = cc.mean_diff
    s = summ(cc); s["contrast"] = "C_inflamed"; Ssum.append(s)
S_all = pd.concat(Ssum, ignore_index=True)
S_all.to_csv(os.path.join(OUT, "partD_%s_geneset_summary.tsv" % a.tag), sep="\t",
             index=False, float_format="%.5g")


def report(tag, S, tbl):
    R = S[(S.contrast == tag) & (S.kind == "random")]
    Re = S[(S.contrast == tag) & (S.kind == "real")]
    P = S[(S.contrast == tag) & (S.kind == "positive")]
    n_ct = int(S[S.contrast == tag].n_ct.iloc[0]) if len(S[S.contrast == tag]) else 0
    print("\n" + "=" * 84)
    print("口径 %s（%d 个细胞类型）" % (tag, n_ct))
    print("=" * 84)
    print("真实面板：")
    print(Re.to_string(index=False, float_format=lambda v: "%.4f" % v))
    print("阳性对照（炎症通路，应升高）：")
    print(P.to_string(index=False, float_format=lambda v: "%.4f" % v))
    print("随机基因集零分布（N=%d）：显著数中位 %.0f（期望 %.1f，膨胀 %.1f×）| |Δ| 中位 %.4f" % (
        len(R), R.n_p05.median(), n_ct * 0.05,
        R.n_p05.median() / max(n_ct * 0.05, 1e-9), R.median_abs_d.median()))
    for _, r in Re.iterrows():
        p05 = np.sort(R.n_p05.values)
        ep = (np.sum(p05 >= r.n_p05) + 1) / (len(p05) + 1)
        print("  %-11s 显著 %2d（经验 p=%.3f）| 最小 p=%.3g | 效应量分位 %.0f%%" % (
            r.geneset, r.n_p05, ep, r.min_p, 100 * np.mean(R.median_abs_d <= r.median_abs_d)))


report("A_naive", S_all, tA)
report("B_primary", S_all, tB)
if len(C):
    report("C_inflamed", S_all, C)

out = {"tag": a.tag, "label": a.label, "path": a.path, "matrix": mname,
       "downsample_to": a.downsample, "n_cells_used": int(n_sel),
       "n_rand": len(rand_names), "sgrg_K": int(K),
       "cell_types_A_naive": ctA, "cell_types_B_primary": ctB,
       "cell_types_C_inflamed": [c for c, _ in okc],
       "n_donors_by_disease": {str(k): int(v) for k, v in
                               df.groupby(["disease", "_grp"]).donor.nunique().items()},
       "geneset_summary": json.loads(S_all.to_json(orient="records"))}
json.dump(out, open(os.path.join(OUT, "partD_%s_summary.json" % a.tag), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False, default=str)
log("完成 %.1f s" % (time.time() - T0))
f.close()
