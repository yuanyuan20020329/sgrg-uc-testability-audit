# -*- coding: utf-8 -*-
"""Part D / 数据集1：逐基因 pseudobulk + 深度混杂诊断"""
import json
import os
import sys
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from partD_common import load_h5ad, normalize_log1p, is_integer_counts, resolve, bh, perm_test_paired_diff, cohen_d

T0 = time.time()
DS = "ds1"
MIN_CELLS = 20
log = lambda *a: print("[%6.1fs]" % (time.time() - T0), *a, flush=True)

X, obs, var = load_h5ad("data/raw/sc/CxG_IBD_HC_10x.h5ad")
sym2idx = {s: i for i, s in enumerate(var["symbol"].values)}
panel = [l.split("\t")[0] for l in
         open("results/sgrg_gene_set_clean.tsv", encoding="utf-8").read().splitlines()[1:] if l.strip()]
Xn = normalize_log1p(X, 1e4)
log("载入完成", X.shape)

donor = obs["donor_id"].astype(str).values
disease = obs["disease"].astype(str).values
ct = obs["cell_type"].astype(str).values

# ---------- per-donor×cell_type 组均值（用指示矩阵一次算完） ----------
key = np.array(["%s|%s" % (d, c) for d, c in zip(donor, ct)])
ukey, kidx = np.unique(key, return_inverse=True)
G = sp.csr_matrix((np.ones(len(key), np.float32), (np.arange(len(key)), kidx)),
                  shape=(len(key), len(ukey)))
gsize = np.asarray(G.sum(axis=0)).ravel()
Gm = sp.diags(1.0 / gsize) @ G.T          # (n_group, n_cells)
sums = np.asarray((Gm @ Xn).todense())     # (n_group, n_genes)
log("组均值矩阵", sums.shape)

meta = pd.DataFrame({"key": ukey})
meta["donor"] = [k.split("|")[0] for k in meta["key"]]
meta["cell_type"] = [k.split("|")[1] for k in meta["key"]]
dmap = {d: dd for d, dd in zip(donor, disease)}
meta["disease"] = meta["donor"].map(dmap)
meta["n_cells"] = gsize
keep = meta["n_cells"] >= MIN_CELLS
log("保留组数 %d / %d" % (keep.sum(), len(meta)))

# ---------- 逐基因 pseudobulk 检验 ----------
gene_idx = resolve(panel, sym2idx)
gene_idx = gene_idx[sorted(range(len(gene_idx)), key=lambda i: panel.index(var["symbol"].values[gene_idx[i]]))]
gene_sym = [var["symbol"].values[i] for i in gene_idx]
gnames = [sym2idx[g] for g in gene_sym]
gene_sym = [var["symbol"].values[i] for i in gnames]

rows = []
for c in sorted(meta.loc[keep, "cell_type"].unique()):
    m = keep & (meta["cell_type"] == c)
    sub = meta[m]
    mat = sums[m.values][:, gnames]
    for gi, g in enumerate(gene_sym):
        v = mat[:, gi]
        uc = v[(sub["disease"] == "ulcerative colitis").values]
        hc = v[(sub["disease"] == "normal").values]
        cd = v[(sub["disease"] == "Crohn disease").values]
        rec = {"cell_type": c, "gene": g, "mean_expr": v.mean(),
               "n_uc": len(uc), "n_hc": len(hc), "n_cd": len(cd),
               "uc": uc.mean(), "hc": hc.mean(), "cd": cd.mean()}
        for lab, a, b in [("UC_vs_normal", uc, hc), ("CD_vs_normal", cd, hc)]:
            if len(a) >= 2 and len(b) >= 2:
                d_, p_, _ = perm_test_paired_diff(a, b, n_perm=20000, seed=3)
                rec["diff_" + lab] = d_
                rec["d_" + lab] = cohen_d(a, b)
                rec["p_" + lab] = p_
            else:
                rec["diff_" + lab] = rec["d_" + lab] = rec["p_" + lab] = np.nan
        rows.append(rec)
pg = pd.DataFrame(rows)
for lab in ["UC_vs_normal", "CD_vs_normal"]:
    pg["q_" + lab] = np.nan
    for c, idx in pg.groupby("cell_type").groups.items():
        pg.loc[idx, "q_" + lab] = bh(pg.loc[idx, "p_" + lab].values)
pg.to_csv("results/partD_%s_sgrg_gene_pseudobulk.tsv" % DS, sep="\t", index=False, float_format="%.4g")
log("逐基因结果写出", pg.shape)

print("\n===== UC vs normal：按 |Cohen's d| 排序的 SGRG 基因（各细胞类型）=====")
for c in sorted(pg["cell_type"].unique()):
    s = pg[pg.cell_type == c].dropna(subset=["d_UC_vs_normal"])
    s = s.reindex(s["d_UC_vs_normal"].abs().sort_values(ascending=False).index).head(6)
    print("\n-- %s" % c)
    print(s[["gene", "hc", "uc", "diff_UC_vs_normal", "d_UC_vs_normal", "p_UC_vs_normal",
             "q_UC_vs_normal"]].to_string(index=False))

# ---------- 深度混杂诊断 ----------
cell = pd.read_csv("results/partD_%s_cell_scores.tsv" % DS, sep="\t")
print("\n===== 深度混杂诊断：SGRG 评分 vs log10(nUMI) =====")
diag = []
for c in sorted(cell["cell_type"].unique()):
    s = cell[cell.cell_type == c]
    r, p = spearmanr(s["SGRG"], s["log_nUMI"])
    r2, p2 = spearmanr(s["SGRG_z"], s["log_nUMI"])
    diag.append({"cell_type": c, "n": len(s), "rho_score_nUMI": r, "rho_z_nUMI": r2,
                 "median_nUMI_UC": s[s.disease == "ulcerative colitis"]["nUMI"].median(),
                 "median_nUMI_normal": s[s.disease == "normal"]["nUMI"].median()})
dd = pd.DataFrame(diag)
print(dd.to_string(index=False))

# 组均值层面
pb = pd.read_csv("results/partD_%s_pseudobulk_score.tsv" % DS, sep="\t")
pb = pb[pb.n_cells >= MIN_CELLS]
print("\n--- pseudobulk 层面 ---")
for c in sorted(pb["cell_type"].unique()):
    s = pb[pb.cell_type == c]
    r, p = spearmanr(s["SGRG"], s["log_nUMI"])
    print("  %-45s rho=%+.3f p=%.3f n=%d" % (c, r, p, len(s)))

dd.to_csv("results/partD_%s_depth_diagnostics.tsv" % DS, sep="\t", index=False, float_format="%.4g")

# ---------- 深度校正后的 UC vs normal（组均值对 log_nUMI 回归取残差） ----------
print("\n===== 深度校正后 UC vs normal（残差法）=====")
res_rows = []
for c in sorted(pb["cell_type"].unique()):
    s = pb[pb.cell_type == c].copy()
    if len(s) < 6:
        continue
    x = s["log_nUMI"].values
    y = s["SGRG"].values
    beta = np.polyfit(x, y, 1)
    res = y - np.polyval(beta, x)
    s["SGRG_resid"] = res
    uc = s[s.disease == "ulcerative colitis"]["SGRG_resid"].values
    hc = s[s.disease == "normal"]["SGRG_resid"].values
    cd_arr = s.loc[s["disease"] == "Crohn disease", "SGRG_resid"].values
    d1, p1, _ = perm_test_paired_diff(uc, hc, n_perm=20000, seed=5)
    d2, p2, _ = perm_test_paired_diff(cd_arr, hc, n_perm=20000, seed=5)
    res_rows.append({"cell_type": c, "n_uc": len(uc), "n_hc": len(hc),
                     "resid_uc": uc.mean(), "resid_hc": hc.mean(), "diff_UC": d1, "p_UC": p1,
                     "diff_CD": d2, "p_CD": p2})
    print("  %-44s UC-HC=%+.4f p=%.3f | CD-HC=%+.4f p=%.3f" % (c, d1, p1, d2, p2))
pd.DataFrame(res_rows).to_csv("results/partD_%s_depth_adjusted_tests.tsv" % DS, sep="\t",
                              index=False, float_format="%.4g")
log("完成")
