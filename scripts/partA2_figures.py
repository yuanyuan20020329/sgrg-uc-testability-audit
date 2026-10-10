# -*- coding: utf-8 -*-
"""A2 图件：公共 UC 队列池的批次质控标注 × 下游数据侧后果（SCI 版式）"""
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.stats import spearmanr, mannwhitneyu

rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "svg.fonttype": "none",
})
MM = 1 / 25.4
ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT = ROOT + "/results"

df = pd.read_csv(OUT + "/partA2_cohorts.tsv", sep="\t")
A = json.load(open(OUT + "/partA2_analysis.json", encoding="utf-8"))
df["log2n"] = np.log2(pd.to_numeric(df["n_samples"], errors="coerce"))

ORDER = ["批次伪影", "信息不足", "无问题", "已校正"]
EN = {"批次伪影": "Batch artifact", "信息不足": "Insufficient batch info",
      "无问题": "No batch problem", "已校正": "Batch corrected"}
ENW = {"批次伪影": "Batch\nartifact", "信息不足": "Insufficient\nbatch info",
       "无问题": "No batch\nproblem", "已校正": "Batch\ncorrected"}
ENSH = {"批次伪影": "Artifact", "信息不足": "Insuff.\ninfo",
        "无问题": "No\nproblem", "已校正": "Corrected"}
COL = {"批次伪影": "#C44E52", "信息不足": "#DD8452", "无问题": "#4C72B0", "已校正": "#55A868"}


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig("%s/%s.%s" % (OUT, name, ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[fig]", name)


def star(p):
    if p is None or not np.isfinite(p):
        return "n.s."
    return "***" if p < 1e-3 else ("**" if p < 1e-2 else ("*" if p < 5e-2 else "n.s."))


def boxplot(ax, groups, values, colors, jitter=True, width=0.55):
    data = [values[g].dropna().values for g in groups]
    bp = ax.boxplot(data, widths=width, patch_artist=True, showfliers=False,
                    medianprops=dict(color="white", lw=1.1),
                    whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
    for patch, g in zip(bp["boxes"], groups):
        patch.set_facecolor(colors[g])
        patch.set_edgecolor("#333333")
        patch.set_linewidth(0.6)
        patch.set_alpha(0.85)
    if jitter:
        rng = np.random.default_rng(3)
        for i, g in enumerate(groups):
            v = values[g].dropna().values
            ax.scatter(np.full(len(v), i + 1) + rng.uniform(-0.16, 0.16, len(v)), v,
                       s=5, color="#2b2b2b", alpha=0.55, linewidths=0, zorder=3)
    return bp


# ============================================================ F1 地景
fig, axes = plt.subplots(1, 3, figsize=(174 * MM, 48 * MM))
cnt = df["annot_status"].value_counts().reindex(ORDER).fillna(0)
ax = axes[0]
y = np.arange(len(ORDER))[::-1]
ax.barh(y, cnt.values, 0.62, color=[COL[g] for g in ORDER], edgecolor="#333", linewidth=0.5)
for yi, g in zip(y, ORDER):
    n = int(cnt[g])
    ax.text(n + 1.2, yi, "%d (%.0f%%)" % (n, 100 * n / len(df)), va="center", fontsize=6.2)
ax.set_yticks(y)
ax.set_yticklabels([EN[g] for g in ORDER], fontsize=6.5)
ax.set_xlabel("Number of cohorts (n = %d)" % len(df))
ax.set_xlim(0, cnt.max() * 1.34)
ax.set_title("a  Pool composition", loc="left", fontsize=7.5, fontweight="bold")

ax = axes[1]
vals = {g: df.loc[df["annot_status"] == g, "geeq_public"] for g in ORDER}
boxplot(ax, ORDER, vals, COL)
ax.set_ylabel("GEEQ public quality score")
ax.set_xticks(range(1, len(ORDER) + 1))
ax.set_xticklabels([ENSH[g] for g in ORDER], fontsize=6)
ax.set_ylim(-0.05, 1.08)
ax.axhline(0.5, color="#999", lw=0.6, ls=":")
ax.set_title("b  Data-derived QC score", loc="left", fontsize=7.5, fontweight="bold")

ax = axes[2]
frac = df.groupby("annot_status")["batch_confound"].apply(
    lambda s: s.fillna("").str.len().gt(0).mean()).reindex(ORDER).fillna(0)
ax.bar(range(len(ORDER)), frac.values * 100, 0.62, color=[COL[g] for g in ORDER],
       edgecolor="#333", linewidth=0.5)
for i, g in enumerate(ORDER):
    ax.text(i, frac[g] * 100 + 2, "%.0f%%" % (frac[g] * 100), ha="center", fontsize=6.2)
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels([ENSH[g] for g in ORDER], fontsize=6)
ax.set_ylabel("Cohorts with batch-confound\nannotation (%)")
ax.set_ylim(0, 108)
ax.set_title("c  Confounding flagged", loc="left", fontsize=7.5, fontweight="bold")
fig.tight_layout(w_pad=1.6)
save(fig, "fig_A2_1_pool_landscape")

# ============================================================ F2 核心：PC1 批次可分性
sv = df.dropna(subset=["auc_pc1_batch", "auc_pc1_primary"]).copy()
sv["normalized_batch_auc"] = sv["auc_pc1_batch"] - sv["auc_pc1_primary"]
fig, axes = plt.subplots(1, 2, figsize=(120 * MM, 50 * MM))
ax = axes[0]
vals = {g: sv.loc[sv["annot_status"] == g, "auc_pc1_batch"] for g in ORDER}
boxplot(ax, ORDER, vals, COL)
ax.axhline(0.5, color="#999", lw=0.6, ls=":")
ax.set_ylabel("AUC of PC1 for batch factor")
ax.set_xticks(range(1, len(ORDER) + 1))
ax.set_xticklabels([ENW[g] for g in ORDER], fontsize=6)
ax.set_ylim(0.42, 1.04)
for i, g in enumerate(ORDER):
    ax.text(i + 1, 1.005, "n=%d" % vals[g].notna().sum(), ha="center", fontsize=5.8, color="#555")
# p 值
def pw(a, b):
    x = sv.loc[sv["annot_status"] == a, "auc_pc1_batch"].dropna()
    y2 = sv.loc[sv["annot_status"] == b, "auc_pc1_batch"].dropna()
    if len(x) < 3 or len(y2) < 3:
        return None
    return float(mannwhitneyu(x, y2).pvalue)
p_ac = pw("批次伪影", "已校正")
p_an = pw("批次伪影", "无问题")
ax.plot([1, 1, 4, 4], [0.955, 0.975, 0.975, 0.955], color="#333", lw=0.7)
ax.text(2.5, 0.978, "%s  (p = %.1e)" % (star(p_ac), p_ac), ha="center", fontsize=6)
ax.plot([1, 1, 3, 3], [0.915, 0.932, 0.932, 0.915], color="#333", lw=0.7)
ax.text(2, 0.935, star(p_an), ha="center", fontsize=6)
ax.set_title("a  Batch structure on the leading axis", loc="left", fontsize=7.5, fontweight="bold")

ax = axes[1]
vals2 = {g: sv.loc[sv["annot_status"] == g, "normalized_batch_auc"] for g in ORDER}
boxplot(ax, ORDER, vals2, COL)
ax.axhline(0, color="#333", lw=0.8, ls="--")
ax.set_ylabel("AUC$_{batch}$ − AUC$_{biology}$  (PC1)")
ax.set_xticks(range(1, len(ORDER) + 1))
ax.set_xticklabels([ENW[g] for g in ORDER], fontsize=6)
ax.set_title("b  Batch minus biology on PC1", loc="left", fontsize=7.5, fontweight="bold")
for i, g in enumerate(ORDER):
    ax.text(i + 1, 0.66, "n=%d" % vals2[g].notna().sum(), ha="center", fontsize=5.8, color="#555")
ax.set_ylim(-0.72, 0.72)
fig.tight_layout(w_pad=1.8)
save(fig, "fig_A2_2_batch_on_pc1")

# ============================================================ F3 生物学可分性（阴性结果）
fig, ax = plt.subplots(figsize=(88 * MM, 48 * MM))
vals3 = {g: sv.loc[sv["annot_status"] == g, "auc_pc1_primary"] for g in ORDER}
boxplot(ax, ORDER, vals3, COL)
ax.axhline(0.5, color="#999", lw=0.6, ls=":")
ax.set_ylabel("AUC of PC1 for the biological factor")
ax.set_xticks(range(1, len(ORDER) + 1))
ax.set_xticklabels([ENW[g] for g in ORDER], fontsize=6)
ax.set_ylim(0.42, 1.04)
m = A["models_binary"].get("auc_pc1_primary")
if m and m.get("F"):
    ax.set_title("Biological signal on PC1 is preserved\n(Kruskal/F$_{perm}$ p = %.2f)" % m["F"]["p_perm"],
                 loc="left", fontsize=7.5, fontweight="bold")
fig.tight_layout()
save(fig, "fig_A2_3_biology_preserved")

# ============================================================ F4 不一致表 + 阈值敏感性
fig, axes = plt.subplots(1, 2, figsize=(132 * MM, 52 * MM))
ax = axes[0]
ts = pd.read_csv(OUT + "/partA2_threshold_sensitivity.tsv", sep="\t")
ax.plot(ts["thr"], ts["sens"], "-o", ms=3, lw=1.1, color="#4C72B0", label="Sensitivity")
ax.plot(ts["thr"], ts["spec"], "-s", ms=3, lw=1.1, color="#55A868", label="Specificity")
ax.plot(ts["thr"], ts["ppv"], "-^", ms=3, lw=1.1, color="#C44E52", label="PPV")
ax.set_xlabel("Threshold: AUC$_{batch}$ for calling batch-dominated")
ax.set_ylabel("Proportion")
ax.set_ylim(0, 1.05)
ax.legend(fontsize=6)
ax.set_title("a  Annotation vs data-derived outcome", loc="left", fontsize=7.5, fontweight="bold")

ax = axes[1]
dd = A["discordance"]["table"]
mat = np.array([[dd["both_pos"], dd["annot_only"]], [dd["svd_only"], dd["both_neg"]]], dtype=float)
im = ax.imshow(mat, cmap="Blues", vmin=0, vmax=mat.max() * 1.1)
for i in range(2):
    for j in range(2):
        ax.text(j, i, "%d" % mat[i, j], ha="center", va="center", fontsize=8,
                color="white" if mat[i, j] > mat.max() * 0.6 else "#222")
ax.set_xticks([0, 1]); ax.set_xticklabels(["Data: yes", "Data: no"], fontsize=6.5)
ax.set_yticks([0, 1]); ax.set_yticklabels(["Annot.: yes", "Annot.: no"], fontsize=6.5)
ax.set_title("b  Discordance (threshold 0.75)", loc="left", fontsize=7.5, fontweight="bold")
for s in ax.spines.values():
    s.set_visible(False)
ax.text(1.02, -0.62, "n = %d" % A["discordance"]["n"], fontsize=6, color="#555")
fig.tight_layout(w_pad=1.8)
save(fig, "fig_A2_4_discordance")

# ============================================================ F5 GEEQ 交叉验证
gq = df.dropna(subset=["geeq_public", "auc_pc1_batch"])
fig, ax = plt.subplots(figsize=(80 * MM, 52 * MM))
for g in ORDER:
    sub = gq[gq["annot_status"] == g]
    ax.scatter(sub["geeq_public"], sub["auc_pc1_batch"], s=11, color=COL[g],
               alpha=0.8, linewidths=0.3, edgecolors="white", label=EN[g])
rho, p = spearmanr(gq["geeq_public"], gq["auc_pc1_batch"])
z = np.polyfit(gq["geeq_public"], gq["auc_pc1_batch"], 1)
xs = np.linspace(gq["geeq_public"].min(), gq["geeq_public"].max(), 20)
ax.plot(xs, np.polyval(z, xs), color="#333", lw=0.9, ls="--")
ax.axhline(0.5, color="#999", lw=0.6, ls=":")
ax.set_xlabel("GEEQ public quality score")
ax.set_ylabel("AUC of PC1 for batch factor")
ax.set_title("$\\rho$ = %.2f (p = %.3g), n = %d" % (rho, p, len(gq)),
             loc="left", fontsize=7.5, fontweight="bold")
ax.legend(fontsize=5.6, loc="upper right")
fig.tight_layout()
save(fig, "fig_A2_5_geeq_crosscheck")

# ============================================================ F6 UC 专属子集
uc = pd.read_csv(OUT + "/partA2_uc_subset.tsv", sep="\t")
uc["d_auc"] = uc["auc_pc1_disease"].fillna(uc["auc_pc1_primary"])
fig, axes = plt.subplots(1, 2, figsize=(132 * MM, 50 * MM))
ax = axes[0]
order2 = [g for g in ORDER if (uc["annot_status"] == g).sum() >= 2]
vals4 = {g: uc.loc[uc["annot_status"] == g, "d_auc"] for g in order2}
boxplot(ax, order2, vals4, COL)
ax.axhline(0.5, color="#999", lw=0.6, ls=":")
ax.set_ylabel("AUC of PC1 for UC vs control")
ax.set_xticks(range(1, len(order2) + 1))
ax.set_xticklabels([ENW[g] for g in order2], fontsize=6)
ax.set_ylim(0.4, 1.05)
for i, g in enumerate(order2):
    ax.text(i + 1, 1.02, "n=%d" % vals4[g].notna().sum(), ha="center", fontsize=5.8, color="#555")
ax.set_title("a  UC-vs-control cohorts (arms $\\geq$ 4)", loc="left", fontsize=7.5, fontweight="bold")

ax = axes[1]
ax.scatter(uc["n_samples"], uc["d_auc"], s=14, c=[COL[g] for g in uc["annot_status"]],
           edgecolors="white", linewidths=0.4)
ax.axhline(0.5, color="#999", lw=0.6, ls=":")
ax.set_xscale("log")
import matplotlib.ticker as mticker
ax.xaxis.set_major_locator(mticker.LogLocator(base=10, numticks=5))
ax.xaxis.set_major_formatter(mticker.FuncFormatter(
    lambda v, _: ("10$^{%d}$" % round(np.log10(v))) if abs(np.log10(v) - round(np.log10(v))) < 0.01 else ""))
ax.xaxis.set_minor_locator(mticker.NullLocator())
ax.set_xlabel("Cohort size (bioassays, log$_{10}$)")
ax.set_ylabel("AUC of PC1 for UC vs control")
ax.set_title("b  Only %d of %d cohorts resolve UC vs control" % (len(uc), len(df)),
             loc="left", fontsize=7.5, fontweight="bold")
fig.tight_layout(w_pad=1.8)
save(fig, "fig_A2_6_uc_subset")

print("\n图件完成。核心数字：")
print("  PC1 批次可分性 伪影组中位 %.3f vs 无问题 %.3f vs 已校正 %.3f" % (
    sv.loc[sv.annot_status == "批次伪影", "auc_pc1_batch"].median(),
    sv.loc[sv.annot_status == "无问题", "auc_pc1_batch"].median(),
    sv.loc[sv.annot_status == "已校正", "auc_pc1_batch"].median()))
print("  Fisher p (伪影 vs 已校正) = %.3g" % p_ac)
