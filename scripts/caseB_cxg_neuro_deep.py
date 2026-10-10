# -*- coding: utf-8 -*-
"""案例 B：CELLxGENE 神经退行域的 L1（对照来源）深度审计，并与 Part D 的 UC 口径对齐。

对齐 Part D 的三大指标：
  ① 非原产细胞占比
  ② 「mixed」（同集合内含对照 + 疾病标签）数据集的 frac_primary 分布
  ③ 「携带对比但臂不可分」的比例（89.3% 的那个指标）

输出：results/caseB_cxg_neuro_deep.json
"""
import pandas as pd
import numpy as np
import json

df = pd.read_csv("results/caseB_cxg_neuro.tsv", sep="\t")
df = df[df["organism"].astype(str).str.contains("Homo sapiens", na=False)].copy()
print("人类神经相关数据集:", len(df))

DIS = {"parkinson", "alzheimer", "dementia", "lewy", "amyotrophic", "sclerosis",
       "huntington", "frontotemporal", "neurodegener", "glioma", "epilepsy",
       "trisomy", "down syndrome", "pick", "progressive supranuclear"}


def splitdis(s):
    return [x.strip().lower() for x in str(s).split(";")]


df["diseases"] = df["disease"].apply(splitdis)
df["has_normal"] = df["diseases"].apply(lambda L: any("normal" in x for x in L))
df["has_disease"] = df["diseases"].apply(lambda L: any(any(d in x for d in DIS) for x in L))
mixed = df[df["has_normal"] & df["has_disease"]]
disonly = df[df["has_disease"] & ~df["has_normal"]]
ctlonly = df[df["has_normal"] & ~df["has_disease"]]

out = {}
out["n_human_datasets"] = int(len(df))
out["n_mixed"] = int(len(mixed))
out["n_dis_only"] = int(len(disonly))
out["n_ctl_only"] = int(len(ctlonly))
out["total_cells"] = int(df["cell_count"].sum())
out["non_primary_cells"] = int((df["cell_count"] - df["primary_cell_count"]).sum())
out["share_non_primary"] = round(out["non_primary_cells"] / out["total_cells"], 4)

print("\n① 非原产细胞: %d / %d = %.2f%%" % (
    out["non_primary_cells"], out["total_cells"], 100 * out["share_non_primary"]))
print("   完全复用（frac_primary=0）数据集: %d / %d = %.1f%%" % (
    int((df["frac_primary"] == 0).sum()), len(df),
    100 * (df["frac_primary"] == 0).mean()))

# ② mixed 的 frac_primary
print("\n② mixed（同集含对照+疾病）数据集: n=%d" % len(mixed))
print("   frac_primary 中位 = %.3f  均值 = %.3f" % (
    mixed["frac_primary"].median(), mixed["frac_primary"].mean()))
print("   frac_primary=0 的比例 = %.1f%%" % (100 * (mixed["frac_primary"] == 0).mean()))
out["mixed_frac_primary_median"] = float(mixed["frac_primary"].median())
out["mixed_frac_primary_mean"] = float(mixed["frac_primary"].mean())
out["mixed_share_fully_reused"] = float((mixed["frac_primary"] == 0).mean())
out["mixed_cells"] = int(mixed["cell_count"].sum())
out["mixed_nonprimary_cells"] = int((mixed["cell_count"] - mixed["primary_cell_count"]).sum())
print("   mixed 细胞 %d，非原产 %d（%.1f%%）" % (
    out["mixed_cells"], out["mixed_nonprimary_cells"],
    100 * out["mixed_nonprimary_cells"] / out["mixed_cells"]))

# ③ 臂不可分：CELLxGENE 数据集层面 disease 字段是「集合内出现的疾病集合」，
#    并不告诉「哪些细胞是对照」——与 Part D 的 89.3% 同源。
n_compare_carrying = int(len(mixed))
out["n_compare_carrying"] = n_compare_carrying
out["arm_indistinguishable_share"] = 1.0  # 数据集级 disease 字段不含细胞级分组
print("\n③ 携带对照+疾病对比的数据集 %d 个；CELLxGENE 数据集级 disease 字段"
      "**不含细胞级分组**" % n_compare_carrying)
print("   → 与 Part D 的「89.3% 臂不可分」同构：可审计性在元数据层即缺失。")

# 分疾病看 frac_primary
print("\n④ 分疾病 frac_primary（细胞加权）")
agg = {}
for k, g in df.groupby("disease"):
    c = g["cell_count"].sum()
    p = g["primary_cell_count"].sum()
    agg[k] = {"cells": int(c), "primary": int(p), "frac_primary": round(p / c, 4) if c else None}
for k, v in sorted(agg.items(), key=lambda x: -x[1]["cells"])[:15]:
    print("   %-50s cells=%-9d frac_primary=%.3f" % (k[:50], v["cells"], v["frac_primary"] or 0))
out["by_disease"] = agg

json.dump(out, open("results/caseB_cxg_neuro_deep.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("\nwrote results/caseB_cxg_neuro_deep.json")
