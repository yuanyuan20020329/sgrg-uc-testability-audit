# -*- coding: utf-8 -*-
"""Part D 通用单细胞分析：SGRG 面板/子集打分 + pseudobulk 检验 + 深度混杂诊断

用法:
  python partD_analyse_generic.py --path X.h5ad --tag ds2 \
     --tissue-keep "colon,mucosa,rectum" --min-cells 20
"""
import argparse
import json
import os
import sys
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from partD_common import (load_h5ad, normalize_log1p, is_integer_counts, score_gene_set,
                          load_gmt, resolve, bh, perm_test_paired_diff, cohen_d)

ap = argparse.ArgumentParser()
ap.add_argument("--path", required=True)
ap.add_argument("--tag", required=True)
ap.add_argument("--tissue-keep", default="")
ap.add_argument("--disease-keep", default="ulcerative colitis,Crohn disease,normal,inflammatory bowel disease")
ap.add_argument("--min-cells", type=int, default=20)
ap.add_argument("--min-donors", type=int, default=3)
ap.add_argument("--label", default="")
ap.add_argument("--write-cells", type=int, default=150000)
a = ap.parse_args()

T0 = time.time()
log = lambda *x: print("[%6.1fs]" % (time.time() - T0), *x, flush=True)
OUT = "results"

log("读取", a.path)
X, obs, var, state = load_h5ad(a.path)
log("shape", X.shape, "矩阵状态 =", state)
sym2idx = {s: i for i, s in enumerate(var["symbol"].values)}

# ---------------- 过滤 ----------------
donor = obs["donor_id"].astype(str).values
disease = obs["disease"].astype(str).values
ct = obs["cell_type"].astype(str).values
tissue = obs["tissue"].astype(str).values if "tissue" in obs.columns else np.array([""] * X.shape[0])

keep = np.ones(X.shape[0], dtype=bool)
if a.tissue_keep:
    pats = [p.strip().lower() for p in a.tissue_keep.split(",") if p.strip()]
    keep &= np.array([any(p in t.lower() for p in pats) for t in tissue])
if a.disease_keep:
    dk = {p.strip() for p in a.disease_keep.split(",") if p.strip()}
    keep &= np.array([d in dk for d in disease])
log("过滤后细胞 %d / %d；组织=[%s] 疾病=[%s]" % (
    keep.sum(), len(keep),
    ", ".join(sorted(set(tissue[keep]))[:6]), ", ".join(sorted(set(disease[keep])))))

if keep.sum() < 200:
    log("细胞过少，退出"); sys.exit(1)
X = X[keep].tocsr()
donor, disease, ct, tissue = donor[keep], disease[keep], ct[keep], tissue[keep]
log("过滤后 shape", X.shape)

if state == "counts":
    nUMI = np.asarray(X.sum(axis=1)).ravel().astype(np.int32)
    Xn = normalize_log1p(X, 1e4)
else:
    nUMI = np.asarray(np.expm1(X).sum(axis=1)).ravel().astype(np.int32)  # 由 log 标准化矩阵反推深度
    Xn = X
log("深度(中位) = %d" % np.median(nUMI))
mean_expr = np.asarray(Xn.mean(axis=0)).ravel().astype(np.float32)

# ---------------- 打分 ----------------
panel = [l.split("\t")[0] for l in
         open("results/sgrg_gene_set_clean.tsv", encoding="utf-8").read().splitlines()[1:] if l.strip()]
go = json.load(open("results/partD_go_sgrg_subsets.json", encoding="utf-8"))
cyto = [g for g in panel if g in set(go["cytoplasmic_stress_granule"])]
nuc = [g for g in panel if g in set(go["nuclear_stress_granule"])]

scores = {}
for name, genes, np_ in [("SGRG", panel, 200), ("SGRG_cyto", cyto, 100), ("SGRG_nuc", nuc, 100)]:
    gi = resolve(genes, sym2idx)
    s, z, _, _ = score_gene_set(Xn, mean_expr, gi, n_perm=np_, seed=7)
    scores[name] = s; scores[name + "_z"] = z
    log("打分 %-11s n=%3d mean=%.4f z中位=%.2f" % (name, len(gi), s.mean(), np.median(z)))

HALL_USE = ["HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
            "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
            "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_HYPOXIA",
            "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
            "HALLMARK_APOPTOSIS", "HALLMARK_MTORC1_SIGNALING"]
hall = load_gmt("data/raw/genesets/h.all.v2023.2.Hs.symbols.gmt")
for h in HALL_USE:
    gi = resolve(hall[h], sym2idx)
    s, _, _, _ = score_gene_set(Xn, mean_expr, gi, n_perm=60, seed=11)
    scores["H_" + h.replace("HALLMARK_", "")] = s

df = pd.DataFrame(scores)
df["donor"] = donor; df["disease"] = disease; df["cell_type"] = ct; df["tissue"] = tissue
df["nUMI"] = nUMI; df["log_nUMI"] = np.log10(nUMI)
df["nGene"] = np.asarray((X > 0).sum(axis=1)).ravel()
if X.shape[0] <= a.write_cells:
    df.to_csv(os.path.join(OUT, "partD_%s_cell_scores.tsv" % a.tag), sep="\t",
              index=False, float_format="%.5g")
    log("写出 per-cell 打分表")

comp = df.groupby(["disease", "cell_type"]).size().unstack(fill_value=0)
comp.to_csv(os.path.join(OUT, "partD_%s_composition.tsv" % a.tag), sep="\t")
log("细胞组成：\n" + comp.to_string())

# ---------------- pseudobulk ----------------
SC = [c for c in df.columns if c.startswith(("SGRG", "H_"))]
pb = df.groupby(["donor", "disease", "cell_type"])[SC + ["nUMI", "log_nUMI"]].mean()
pb["n_cells"] = df.groupby(["donor", "disease", "cell_type"]).size()
pb = pb.reset_index()
pb.to_csv(os.path.join(OUT, "partD_%s_pseudobulk_score.tsv" % a.tag), sep="\t",
          index=False, float_format="%.5g")

pb_f = pb[pb.n_cells >= a.min_cells].copy()
# 组内样本数门槛
valid_ct = pb_f.groupby(["cell_type", "disease"]).size().unstack(fill_value=0)
ok_ct = [c for c in valid_ct.index
         if all(valid_ct.loc[c].get(d, 0) >= a.min_donors for d in ["ulcerative colitis", "normal"])]
log("纳入检验的细胞类型（UC/HC 各 >=%d donor）：%s" % (a.min_donors, ok_ct))

rows = []
for c in ok_ct:
    s = pb_f[pb_f.cell_type == c].copy()
    # 深度校正（组内对 log_nUMI 回归取残差）
    b = np.polyfit(s["log_nUMI"].values, s["SGRG"].values, 1)
    s["SGRG_adj"] = s["SGRG"].values - np.polyval(b, s["log_nUMI"].values)
    uc, hc = s[s.disease == "ulcerative colitis"], s[s.disease == "normal"]
    for sc in ["SGRG", "SGRG_cyto", "SGRG_nuc", "SGRG_adj"]:
        d, p, _ = perm_test_paired_diff(uc[sc].values, hc[sc].values, n_perm=20000, seed=1)
        rows.append({"cell_type": c, "score": sc, "n_uc": len(uc), "n_hc": len(hc),
                     "mean_uc": uc[sc].mean(), "mean_hc": hc[sc].mean(),
                     "diff": d, "cohen_d": cohen_d(uc[sc].values, hc[sc].values), "p_perm": p,
                     "median_nUMI_uc": uc["nUMI"].median(), "median_nUMI_hc": hc["nUMI"].median()})
tt = pd.DataFrame(rows)
tt["q_BH"] = np.nan
for sc, idx in tt.groupby("score").groups.items():
    tt.loc[idx, "q_BH"] = bh(tt.loc[idx, "p_perm"].values)
tt.to_csv(os.path.join(OUT, "partD_%s_uc_vs_hc_tests.tsv" % a.tag), sep="\t",
          index=False, float_format="%.4g")

print("\n===== %s : UC vs normal（pseudobulk，donor 级）=====" % a.tag)
print(tt[["cell_type", "score", "n_uc", "n_hc", "mean_uc", "mean_hc", "diff", "cohen_d",
          "p_perm", "q_BH", "median_nUMI_uc", "median_nUMI_hc"]].to_string(index=False))

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
           "n_cells_used": int(X.shape[0]), "n_genes": int(X.shape[1]),
           "tissue_keep": a.tissue_keep, "disease_keep": a.disease_keep,
           "min_cells": a.min_cells, "min_donors": a.min_donors,
           "cell_types_tested": ok_ct,
           "n_donors_by_disease": df.groupby("disease")["donor"].nunique().to_dict(),
           "baseline_by_celltype": json.loads(base.reset_index().to_json(orient="records")),
           "composition": json.loads(comp.to_json(orient="index"))}
json.dump(summary, open(os.path.join(OUT, "partD_%s_summary.json" % a.tag), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False, default=str)
log("完成 %.1f s" % (time.time() - T0))
