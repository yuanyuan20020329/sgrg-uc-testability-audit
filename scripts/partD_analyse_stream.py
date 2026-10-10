# -*- coding: utf-8 -*-
"""Part D 分块流式单细胞分析（针对超大 h5ad，如 ds4 Tearle 320k 细胞）

与 partD_analyse_generic.py 输出同构（summary / tests / pseudobulk / composition），
但用 h5py 按行分块读取，两遍扫描：pass1 求每基因均值表达与每细胞深度，pass2 打分。

用法:
  python partD_analyse_stream.py --path X.h5ad --tag ds4_tearle --label "Tearle2026" \
      --tissue-keep "colon" --chunk 12000
"""
import argparse
import json
import os
import sys
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from partD_common import (control_matrix, bh, perm_test_paired_diff, cohen_d)

ap = argparse.ArgumentParser()
ap.add_argument("--path", required=True)
ap.add_argument("--tag", required=True)
ap.add_argument("--label", default="")
ap.add_argument("--tissue-keep", default="")
ap.add_argument("--disease-keep", default="ulcerative colitis,Crohn disease,normal,inflammatory bowel disease,primary sclerosing cholangitis")
ap.add_argument("--uc-match", default="ulcer")        # UC 分组匹配（子串，小写）
ap.add_argument("--ctrl-match", default="normal")     # 对照组匹配
ap.add_argument("--matrix", default="auto")           # auto | X | raw/X
ap.add_argument("--chunk", type=int, default=12000)
ap.add_argument("--min-cells", type=int, default=20)
ap.add_argument("--min-donors", type=int, default=3)
ap.add_argument("--write-cells", type=int, default=150000)
ap.add_argument("--source", type=float, default=1e4)
ap.add_argument("--downsample", type=int, default=0,
                help=">0 时把所有细胞抽稀到该 UMI 数（深度匹配），丢弃不足者")
ap.add_argument("--seed", type=int, default=17)
a = ap.parse_args()

T0 = time.time()
log = lambda *x: print("[%6.1fs]" % (time.time() - T0), *x, flush=True)
OUT = "results"
NP_PERM = {"SGRG": 200, "SGRG_cyto": 100, "SGRG_nuc": 100}


# ------------------------------------------------------------------ obs 读取
def obs_col(f, col):
    if col not in f["obs"]:
        return None
    d = f["obs"][col]
    if isinstance(d, h5py.Group):
        cats = [x.decode() if isinstance(x, bytes) else str(x) for x in d["categories"][:]]
        codes = d["codes"][:]
        out = np.empty(len(codes), dtype=object)
        cats = np.array(cats, dtype=object)
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
log("观测数 (n_obs, n_var) =", f["X"].attrs.get("shape") if "X" in f else "?")

# 1) 选矩阵（要求 csr + 整数 counts）
def enc(g):
    return g.attrs.get("encoding-type", g.attrs.get("encoding_type", ""))


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
        log("  候选矩阵 %-10s enc=%-10s 整数占比=%.3f" % (name, enc(g), fr))
        if fr > 0.999 and chosen is None:
            chosen = (name, g)
    if chosen is None and cand:
        chosen = cand[0]

mname, grp = chosen
log("使用矩阵 =", mname, "| 状态 = counts（整数）")
n_all, n_genes = grp.attrs["shape"]
n_all, n_genes = int(n_all), int(n_genes)

# 2) var 基因符号
def read_col(g, col):
    """读 h5ad dataframe 的一列：可能是 Dataset，也可能是分类 Group"""
    d = g[col]
    if isinstance(d, h5py.Group):
        cats = [x.decode() if isinstance(x, bytes) else x for x in d["categories"][:]]
        cats = np.array(cats, dtype=object)
        codes = d["codes"][:]
        out = np.empty(len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        out[~ok] = ""
        return out
    a = d[:]
    if a.dtype.kind in "OS":
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a], dtype=object)
    return a


def var_symbols(vgrp):
    for col in ("gene_symbols", "feature_name", "feature_name_original"):
        if col in vgrp:
            arr = np.asarray(read_col(vgrp, col), dtype=object)
            if len(arr) and not str(arr[0]).startswith("ENSG"):
                return arr
    ic = vgrp.attrs.get("_index")
    if ic and ic in vgrp:
        return np.asarray(read_col(vgrp, ic), dtype=object)
    return np.array(["G%d" % i for i in range(n_genes)], dtype=object)

vgrp = f["raw"]["var"] if mname.startswith("raw/") else f["var"]
symbols = var_symbols(vgrp)
log("基因数 =", n_genes, "| 示例 =", list(symbols[:5]))

# 3) obs 过滤
donor = obs_col(f, "donor_id")
donor = np.array(["d%d" % i for i in range(n_all)], dtype=object) if donor is None else donor
disease = obs_col(f, "disease")
disease = np.array(["unknown"] * n_all, dtype=object) if disease is None else disease
ct = obs_col(f, "cell_type")
ct = np.array(["unknown"] * n_all, dtype=object) if ct is None else ct
tissue = obs_col(f, "tissue")
tissue = np.array([""] * n_all, dtype=object) if tissue is None else tissue

log("疾病分布：", dict(pd.Series(disease).value_counts()))
log("组织分布（前 12）：", dict(pd.Series(tissue).value_counts().head(12)))

keep = np.ones(n_all, dtype=bool)
if a.tissue_keep:
    pats = [p.strip().lower() for p in a.tissue_keep.split(",") if p.strip()]
    keep &= np.array([any(p in str(t).lower() for p in pats) for t in tissue])
if a.disease_keep:
    dk = {p.strip().lower() for p in a.disease_keep.split(",") if p.strip()}
    keep &= np.array([str(d).lower() in dk for d in disease])
log("过滤后细胞 %d / %d" % (keep.sum(), n_all))
if keep.sum() < 200:
    log("细胞过少，退出"); sys.exit(1)

sel = np.where(keep)[0]
donor, disease, ct, tissue = donor[sel], disease[sel], ct[sel], tissue[sel]

# ---- 深度匹配（可选）：先算每细胞总 counts，丢掉低于阈值的细胞 ----
RNG = np.random.default_rng(a.seed)
indptr_full = grp["indptr"][:]


def _totals(idx):
    """快速求每细胞总 counts（只读 data，不读 indices）"""
    t = np.zeros(len(idx), dtype=np.int64)
    ip = indptr_full
    for i0 in range(0, len(idx), 200000):
        block = idx[i0:i0 + 200000]
        s = int(ip[block[0]]); e = int(ip[block[-1] + 1])
        off = np.append(ip[block] - s, e - s)      # 长度 n+1
        seg = np.diff(off)                          # 长度 n
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
log("过滤后疾病分布：", dict(pd.Series(disease).value_counts()))
log("过滤后组织分布：", dict(pd.Series(tissue).value_counts().head(10)))

# 4) 分块读取器
indptr_full = grp["indptr"][:]
row_start = indptr_full[sel]
row_end = indptr_full[sel + 1]
nnz_tot = int(row_end[-1] - row_start[0])
log("选取细胞 nnz = %d（平均 %.0f 基因/细胞）" % (nnz_tot, nnz_tot / n_sel))


def chunk_iter(chunk):
    a0 = 0
    while a0 < n_sel:
        b0 = min(a0 + chunk, n_sel)
        s0 = int(row_start[a0]); s1 = int(row_end[b0 - 1])
        data = grp["data"][s0:s1].astype(np.float32)
        ind = grp["indices"][s0:s1].astype(np.int32)
        lip = np.concatenate([[0], (row_end[a0:b0] - s0)]).astype(np.int64)
        Xc = sp.csr_matrix((data, ind, lip), shape=(b0 - a0, n_genes))
        yield a0, b0, Xc
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
    if (a0 // a.chunk) % 5 == 0:
        log("  pass1 %d / %d" % (b0, n_sel))
mean_expr = (gene_sum / n_sel).astype(np.float32)
log("深度(中位) = %d" % np.median(nUMI))

# ------------------------------------------------------------------ 基因集
panel = [l.split("\t")[0] for l in
         open("results/sgrg_gene_set_clean.tsv", encoding="utf-8").read().splitlines()[1:] if l.strip()]
go = json.load(open("results/partD_go_sgrg_subsets.json", encoding="utf-8"))
cyto = [g for g in panel if g in set(go["cytoplasmic_stress_granule"])]
nuc = [g for g in panel if g in set(go["nuclear_stress_granule"])]
sym2idx = {}
for i, s in enumerate(symbols):
    sym2idx.setdefault(str(s), i)
sym2idx = {s: i for i, s in enumerate(symbols)}


def resolve(genes):
    return np.array(sorted({sym2idx[g] for g in genes if g in sym2idx}), dtype=np.int64)


HALL_USE = ["HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
            "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
            "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_HYPOXIA",
            "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
            "HALLMARK_APOPTOSIS", "HALLMARK_MTORC1_SIGNALING"]


def load_gmt(path):
    sets = {}
    for line in open(path, encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            sets[p[0]] = [g for g in p[2:] if g]
    return sets


hall = load_gmt("data/raw/genesets/h.all.v2023.2.Hs.symbols.gmt")

sets = []          # (输出名, gene_idx, n_perm)
for name, genes in [("SGRG", panel), ("SGRG_cyto", cyto), ("SGRG_nuc", nuc)]:
    gi = resolve(genes)
    gi = gi[gene_det[gi] > 0] if len(gi) else gi
    sets.append((name, gi, NP_PERM[name]))
    log("基因集 %-11s 命中 %3d / %d（可检出 %d）" % (name, len(gi), len(genes), len(gi)))
for h in HALL_USE:
    if h in hall:
        gi = resolve(hall[h])
        gi = gi[gene_det[gi] > 0]
        sets.append(("H_" + h.replace("HALLMARK_", ""), gi, 60))

# 预先构建零分布控制矩阵（分基因集）
Vmap = {}
for name, gi, np_perm in sets:
    mask = np.zeros(n_genes, dtype=bool)
    mask[gi] = True
    Vmap[name] = control_matrix(mean_expr, mask, n_bins=25, ctrl_size=50, n_perm=np_perm, seed=7)
log("控制矩阵构建完成")

# ------------------------------------------------------------------ pass 2
log("pass 2：逐细胞打分")
cols = {name: np.zeros(n_sel, dtype=np.float32) for name, _, _ in sets}
for name in [n for n, _, _ in sets]:
    if not name.startswith("H_"):
        cols[name + "_z"] = np.zeros(n_sel, dtype=np.float32)

nGene = np.zeros(n_sel, dtype=np.int32)
for a0, b0, Xc in chunk_iter(a.chunk):
    Xn, lib = norm_cp10k(Xc)
    nGene[a0:b0] = np.asarray((Xc > 0).sum(axis=1)).ravel().astype(np.int32)
    for name, gi, np_perm in sets:
        V = Vmap[name]
        obs = np.asarray(Xn[:, gi].mean(axis=1)).ravel()
        null = np.asarray((Xn @ V).todense(), dtype=np.float32)
        nm = null.mean(axis=1); nsd = null.std(axis=1, ddof=1)
        sc = obs - nm
        sd = np.where(nsd > 0, nsd, np.nan)
        cols[name][a0:b0] = sc.astype(np.float32)
        if name + "_z" in cols:
            cols[name + "_z"][a0:b0] = (sc / sd).astype(np.float32)
    del Xc, Xn
    if (a0 // a.chunk) % 5 == 0:
        log("  pass2 %d / %d" % (b0, n_sel))
log("打分完成")

df = pd.DataFrame(cols)
df["donor"] = donor; df["disease"] = disease; df["cell_type"] = ct; df["tissue"] = tissue
df["nUMI"] = nUMI; df["log_nUMI"] = np.log10(np.maximum(nUMI, 1)); df["nGene"] = nGene

for name in ["SGRG", "SGRG_cyto", "SGRG_nuc"]:
    log("打分 %-11s mean=%.4f z中位=%.2f" % (name, df[name].mean(), np.median(df[name + "_z"])))

if n_sel <= a.write_cells:
    df.to_csv(os.path.join(OUT, "partD_%s_cell_scores.tsv" % a.tag), sep="\t",
              index=False, float_format="%.5g")
    log("写出 per-cell 打分表")

comp = df.groupby(["disease", "cell_type"]).size().unstack(fill_value=0)
comp.to_csv(os.path.join(OUT, "partD_%s_composition.tsv" % a.tag), sep="\t")
log("细胞组成：\n" + comp.to_string())

# ------------------------------------------------------------------ pseudobulk
SC = [c for c in df.columns if c.startswith(("SGRG", "H_")) and not c.endswith("_z")]
pb = df.groupby(["donor", "disease", "cell_type"])[SC + ["nUMI", "log_nUMI"]].mean()
pb["n_cells"] = df.groupby(["donor", "disease", "cell_type"]).size()
pb = pb.reset_index()
pb.to_csv(os.path.join(OUT, "partD_%s_pseudobulk_score.tsv" % a.tag), sep="\t",
          index=False, float_format="%.5g")

pb_f = pb[pb.n_cells >= a.min_cells].copy()
pb_f["_grp"] = np.where(pb_f.disease.str.lower().str.contains(a.uc_match), "UC",
                        np.where(pb_f.disease.str.lower().str.contains(a.ctrl_match), "CTRL", "OTHER"))
valid_ct = pb_f[pb_f._grp != "OTHER"].groupby(["cell_type", "_grp"]).size().unstack(fill_value=0)
ok_ct = [c for c in valid_ct.index
         if valid_ct.loc[c].get("UC", 0) >= a.min_donors and valid_ct.loc[c].get("CTRL", 0) >= a.min_donors]
log("纳入检验的细胞类型（UC/CTRL 各 >=%d donor）：%s" % (a.min_donors, ok_ct))

rows = []
for c in ok_ct:
    s = pb_f[(pb_f.cell_type == c) & (pb_f._grp != "OTHER")].copy()
    b = np.polyfit(s["log_nUMI"].values, s["SGRG"].values, 1)
    s["SGRG_adj"] = s["SGRG"].values - np.polyval(b, s["log_nUMI"].values)
    uc, hc = s[s._grp == "UC"], s[s._grp == "CTRL"]
    for sc in ["SGRG", "SGRG_cyto", "SGRG_nuc", "SGRG_adj"]:
        d, p, _ = perm_test_paired_diff(uc[sc].values, hc[sc].values, n_perm=20000, seed=1)
        rows.append({"cell_type": c, "score": sc, "n_uc": len(uc), "n_hc": len(hc),
                     "mean_uc": uc[sc].mean(), "mean_hc": hc[sc].mean(),
                     "diff": d, "cohen_d": cohen_d(uc[sc].values, hc[sc].values), "p_perm": p,
                     "median_nUMI_uc": uc["nUMI"].median(), "median_nUMI_hc": hc["nUMI"].median()})
tt = pd.DataFrame(rows)
if len(tt):
    tt["q_BH"] = np.nan
    for sc, idx in tt.groupby("score").groups.items():
        tt.loc[idx, "q_BH"] = bh(tt.loc[idx, "p_perm"].values)
    tt.to_csv(os.path.join(OUT, "partD_%s_uc_vs_hc_tests.tsv" % a.tag), sep="\t",
              index=False, float_format="%.4g")
    print("\n===== %s : UC vs CTRL（pseudobulk，donor 级）=====" % a.tag)
    print(tt.to_string(index=False))

print("\n===== 基线：SGRG 评分按细胞类型（全部细胞）=====")
base = df.groupby("cell_type")["SGRG"].agg(["mean", "median", "count"]).sort_values("mean", ascending=False)
print(base.to_string())

print("\n===== 深度混杂：SGRG vs log10(nUMI) =====")
for c in sorted(df.cell_type.unique()):
    s = df[df.cell_type == c]
    if len(s) < 50:
        continue
    r, p = spearmanr(s["SGRG"], s["log_nUMI"])
    print("  %-50s rho=%+.3f p=%.1e n=%d" % (c, r, p, len(s)))

summary = {"tag": a.tag, "label": a.label, "path": a.path,
           "n_cells_used": int(n_sel), "n_genes": int(n_genes), "matrix": mname,
           "downsample_to": a.downsample,
           "tissue_keep": a.tissue_keep, "disease_keep": a.disease_keep,
           "min_cells": a.min_cells, "min_donors": a.min_donors,
           "cell_types_tested": ok_ct,
           "n_donors_by_disease": {str(k): int(v) for k, v in
                                   df.groupby("disease")["donor"].nunique().items()},
           "baseline_by_celltype": json.loads(base.reset_index().to_json(orient="records")),
           "composition": json.loads(comp.to_json(orient="index"))}
json.dump(summary, open(os.path.join(OUT, "partD_%s_summary.json" % a.tag), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False, default=str)
log("完成 %.1f s" % (time.time() - T0))
f.close()
