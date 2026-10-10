# -*- coding: utf-8 -*-
"""Part D / 数据集1（Garrido-Trigo 2023, Nat Commun; CxG）：
   应激颗粒面板在 UC 黏膜中的细胞类型定位与炎症状态关联
   分析：per-cell 打分 -> pseudobulk（donor×cell_type）-> 置换检验 / 逐基因
"""
import json
import os
import sys
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from partD_common import (load_h5ad, normalize_log1p, is_integer_counts, score_gene_set,
                          load_gmt, resolve, bh, perm_test_paired_diff, cohen_d)

T0 = time.time()
DS = "ds1"
PATH = "data/raw/sc/CxG_IBD_HC_10x.h5ad"
OUT = "results"


def log(*a):
    print("[%6.1fs]" % (time.time() - T0), *a, flush=True)


# ------------------------------------------------------------------ 读入
log("读取", PATH)
X, obs, var = load_h5ad(PATH)
log("shape =", X.shape, " 整数占比 = %.3f" % is_integer_counts(X))
sym2idx = {s: i for i, s in enumerate(var["symbol"].values)}

panel = [l.split("\t")[0] for l in
         open("results/sgrg_gene_set_clean.tsv", encoding="utf-8").read().splitlines()[1:] if l.strip()]
go = json.load(open("results/partD_go_sgrg_subsets.json", encoding="utf-8"))
cyto_genes = [g for g in panel if g in set(go["cytoplasmic_stress_granule"])]
nuc_genes = [g for g in panel if g in set(go["nuclear_stress_granule"])]
log("面板 %d / 胞质 %d / 核 %d" % (len(panel), len(cyto_genes), len(nuc_genes)))

nUMI = np.asarray(X.sum(axis=1)).ravel().astype(np.int32)
nGene = np.asarray((X > 0).sum(axis=1)).ravel().astype(np.int32)

log("归一化 CP10K + log1p")
Xn = normalize_log1p(X, 1e4)
mean_expr = np.asarray(Xn.mean(axis=0)).ravel().astype(np.float32)

# ------------------------------------------------------------------ 打分
scores = {}
for name, genes, nperm in [("SGRG", panel, 200), ("SGRG_cyto", cyto_genes, 100),
                           ("SGRG_nuc", nuc_genes, 100)]:
    gi = resolve(genes, sym2idx)
    s, z, nm, ns = score_gene_set(Xn, mean_expr, gi, n_perm=nperm, seed=7)
    scores[name] = s
    scores[name + "_z"] = z
    log("打分 %-11s n_genes=%3d  mean=%.4f  z(中位)=%.3f" % (name, len(gi), s.mean(), np.median(z)))

hall = load_gmt("data/raw/genesets/h.all.v2023.2.Hs.symbols.gmt")
HALL_USE = ["HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
            "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
            "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_HYPOXIA",
            "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
            "HALLMARK_APOPTOSIS", "HALLMARK_MTORC1_SIGNALING",
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
for h in HALL_USE:
    gi = resolve(hall[h], sym2idx)
    s, z, _, _ = score_gene_set(Xn, mean_expr, gi, n_perm=100, seed=11)
    scores["H_" + h.replace("HALLMARK_", "")] = s
    log("打分 %-45s n=%4d" % (h, len(gi)))

df = pd.DataFrame(scores)
df["donor"] = obs["donor_id"].astype(str).values
df["disease"] = obs["disease"].astype(str).values
df["cell_type"] = obs["cell_type"].astype(str).values
df["tissue"] = obs["tissue"].astype(str).values
df["nUMI"] = nUMI
df["nGene"] = nGene
df["log_nUMI"] = np.log10(nUMI)
df.to_csv(os.path.join(OUT, "partD_%s_cell_scores.tsv" % DS), sep="\t", index=False, float_format="%.5g")
log("写出 per-cell 打分表", df.shape)

# ------------------------------------------------------------------ 细胞组成
comp = df.groupby(["disease", "cell_type"]).size().unstack(fill_value=0)
comp.to_csv(os.path.join(OUT, "partD_%s_composition.tsv" % DS), sep="\t")
donor_ct = df.groupby(["donor", "disease", "cell_type"]).size().unstack(fill_value=0)
log("细胞组成（疾病×细胞类型）：\n" + comp.to_string())

# ------------------------------------------------------------------ pseudobulk
SCORE_COLS = [c for c in df.columns if c.startswith(("SGRG", "H_"))]
pb = df.groupby(["donor", "disease", "cell_type"])[SCORE_COLS + ["nUMI", "log_nUMI"]].mean()
pb["n_cells"] = df.groupby(["donor", "disease", "cell_type"]).size()
pb = pb.reset_index()
pb.to_csv(os.path.join(OUT, "partD_%s_pseudobulk_score.tsv" % DS), sep="\t",
          index=False, float_format="%.5g")
log("pseudobulk 行数 =", len(pb), " (donor×cell_type)")

MIN_CELLS = 20  # 每 donor×cell_type 至少 20 个细胞才纳入
pb_f = pb[pb["n_cells"] >= MIN_CELLS].copy()

test_rows = []
for ct in sorted(pb_f["cell_type"].unique()):
    sub = pb_f[pb_f["cell_type"] == ct]
    uc = sub[sub["disease"] == "ulcerative colitis"]
    hc = sub[sub["disease"] == "normal"]
    cd = sub[sub["disease"] == "Crohn disease"]
    for sc in SCORE_COLS:
        for g1, g2, lab in [(uc, hc, "UC_vs_normal"), (cd, hc, "CD_vs_normal"), (uc, cd, "UC_vs_CD")]:
            if len(g1) < 2 or len(g2) < 2:
                continue
            d, p, _ = perm_test_paired_diff(g1[sc].values, g2[sc].values, n_perm=20000, seed=1)
            test_rows.append({"cell_type": ct, "score": sc, "contrast": lab,
                              "n1": len(g1), "n2": len(g2),
                              "mean1": g1[sc].mean(), "mean2": g2[sc].mean(),
                              "diff": d, "cohen_d": cohen_d(g1[sc].values, g2[sc].values), "p_perm": p})
tt = pd.DataFrame(test_rows)
tt["q_BH_within_score"] = np.nan
for (sc, lab), idx in tt.groupby(["score", "contrast"]).groups.items():
    tt.loc[idx, "q_BH_within_score"] = bh(tt.loc[idx, "p_perm"].values)
tt.to_csv(os.path.join(OUT, "partD_%s_pseudobulk_tests.tsv" % DS), sep="\t",
          index=False, float_format="%.4g")

print("\n===== UC vs normal：SGRG 面板评分（pseudobulk, donor 级）=====")
show = tt[(tt.contrast == "UC_vs_normal") & (tt.score.isin(["SGRG", "SGRG_cyto", "SGRG_nuc"]))]
print(show[["cell_type", "score", "n1", "n2", "mean1", "mean2", "diff", "cohen_d", "p_perm"]]
      .to_string(index=False))
print("\n===== CD vs normal =====")
show2 = tt[(tt.contrast == "CD_vs_normal") & (tt.score.isin(["SGRG", "SGRG_cyto", "SGRG_nuc"]))]
print(show2[["cell_type", "score", "n1", "n2", "mean1", "mean2", "diff", "cohen_d", "p_perm"]]
      .to_string(index=False))

print("\n===== 基线：SGRG 评分按细胞类型（全部细胞）=====")
base = df.groupby("cell_type")["SGRG"].agg(["mean", "median", "count"]).sort_values("mean", ascending=False)
print(base.to_string())
base_z = df.groupby("cell_type")["SGRG_z"].agg(["mean", "median"]).sort_values("mean", ascending=False)
print(base_z.to_string())

print("\n===== 炎症通路评分 vs SGRG（pseudobulk, 相关）=====")
cor_rows = []
for ct in sorted(pb_f["cell_type"].unique()):
    sub = pb_f[pb_f["cell_type"] == ct]
    if len(sub) < 6:
        continue
    for c in SCORE_COLS:
        if not c.startswith("H_"):
            continue
        r = np.corrcoef(sub["SGRG"].values, sub[c].values)[0, 1]
        cor_rows.append({"cell_type": ct, "signature": c, "pearson_r": r, "n": len(sub)})
cr = pd.DataFrame(cor_rows)
cr.to_csv(os.path.join(OUT, "partD_%s_sgrg_vs_hallmark.tsv" % DS), sep="\t",
          index=False, float_format="%.4g")
print(cr.pivot(index="signature", columns="cell_type", values="pearson_r").round(3).to_string())

json.dump({"dataset": DS, "path": PATH, "n_cells": int(X.shape[0]), "n_genes": int(X.shape[1]),
           "panel": len(panel), "cyto": len(cyto_genes), "nuc": len(nuc_genes),
           "min_cells": MIN_CELLS,
           "baseline_by_celltype": json.loads(base.reset_index().to_json(orient="records")),
           "composition": json.loads(comp.to_json(orient="index"))},
          open(os.path.join(OUT, "partD_%s_summary.json" % DS), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
log("完成，总耗时 %.1f s" % (time.time() - T0))
