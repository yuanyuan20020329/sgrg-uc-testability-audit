# -*- coding: utf-8 -*-
"""Part D / 数据集1：细胞类型内「炎症状态 × 应激颗粒」耦合分析（控制测序深度）
   问题：在哪种炎症状态下 SGRG 程序被激活？
"""
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DS = "ds1"
MIN_CELLS = 200
cell = pd.read_csv("results/partD_%s_cell_scores.tsv" % DS, sep="\t")
HL = [c for c in cell.columns if c.startswith("H_")]
print("hallmark 列：", HL)

# 炎症复合评分（各 hallmark 先按全细胞 z 标准化再平均）
INFLAM = ["H_TNFA_SIGNALING_VIA_NFKB", "H_INFLAMMATORY_RESPONSE",
          "H_INTERFERON_GAMMA_RESPONSE", "H_INTERFERON_ALPHA_RESPONSE",
          "H_IL6_JAK_STAT3_SIGNALING"]
z = cell[INFLAM].apply(lambda s: (s - s.mean()) / s.std())
cell["INFLAM"] = z.mean(axis=1)

rows, decile_rows = [], []
for c in sorted(cell.cell_type.unique()):
    s = cell[cell.cell_type == c].copy()
    if len(s) < MIN_CELLS:
        continue
    # 控制深度：对 log_nUMI 取残差
    b = np.polyfit(s["log_nUMI"].values, s["SGRG"].values, 1)
    s["SGRG_res"] = s["SGRG"].values - np.polyval(b, s["log_nUMI"].values)
    b2 = np.polyfit(s["log_nUMI"].values, s["INFLAM"].values, 1)
    s["INFLAM_res"] = s["INFLAM"].values - np.polyval(b2, s["log_nUMI"].values)

    for dis in ["normal", "ulcerative colitis"]:
        g = s[s.disease == dis]
        if len(g) < 100:
            continue
        r_raw, _ = spearmanr(g["SGRG"], g["INFLAM"])
        r_adj, p_adj = spearmanr(g["SGRG_res"], g["INFLAM_res"])
        rows.append({"cell_type": c, "disease": dis, "n": len(g),
                     "rho_raw": r_raw, "rho_depth_adj": r_adj, "p_depth_adj": p_adj})
    # 分十分位（用全组 quintile 更稳）
    s["INFLAM_q"] = pd.qcut(s["INFLAM_res"], 5, labels=[1, 2, 3, 4, 5])
    for q, g in s.groupby("INFLAM_q", observed=True):
        decile_rows.append({"cell_type": c, "quintile": int(q), "n": len(g),
                            "sgrg_mean": g["SGRG_res"].mean(),
                            "sgrg_sem": g["SGRG_res"].std() / np.sqrt(len(g)),
                            "inflam_mean": g["INFLAM_res"].mean()})

t1 = pd.DataFrame(rows)
t2 = pd.DataFrame(decile_rows)
t1.to_csv("results/partD_%s_state_coupling.tsv" % DS, sep="\t", index=False, float_format="%.4g")
t2.to_csv("results/partD_%s_state_quintile.tsv" % DS, sep="\t", index=False, float_format="%.5g")

print("\n===== 细胞类型内：SGRG（深度校正）vs 炎症复合评分（深度校正）的 Spearman =====")
print(t1.pivot(index="cell_type", columns="disease", values="rho_depth_adj").round(3).to_string())
print("\n--- 原始（未校正深度）---")
print(t1.pivot(index="cell_type", columns="disease", values="rho_raw").round(3).to_string())

print("\n===== 炎症复合评分五分位 × 深度校正 SGRG 评分 =====")
pv = t2.pivot(index="cell_type", columns="quintile", values="sgrg_mean")
print(pv.round(4).to_string())
print("\n(对比：五分位区间跨度——Q5-Q1)")
delta = (pv[5] - pv[1]).rename("Q5_minus_Q1").reset_index()
print(delta.round(4).to_string(index=False))

json.dump({"coupling": json.loads(t1.to_json(orient="records")),
           "quintile_trend": json.loads(delta.to_json(orient="records"))},
          open("results/partD_%s_state_summary.json" % DS, "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
print("\n[saved] results/partD_%s_state_coupling.tsv / state_quintile.tsv / state_summary.json" % DS)
