# -*- coding: utf-8 -*-
"""Part D 图件（SCI 版式：mm 单位、Arial、300 dpi、无冗余边框）"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.stats import spearmanr

rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "svg.fonttype": "none",
})
MM = 1 / 25.4
DS = "ds1"
C = {"normal": "#4C72B0", "ulcerative colitis": "#C44E52", "Crohn disease": "#DD8452"}
LAB = {"normal": "Normal", "ulcerative colitis": "UC", "Crohn disease": "CD"}

cell = pd.read_csv("results/partD_%s_cell_scores.tsv" % DS, sep="\t")
pb = pd.read_csv("results/partD_%s_pseudobulk_score.tsv" % DS, sep="\t")
pb = pb[pb.n_cells >= 20]
tests = pd.read_csv("results/partD_%s_pseudobulk_tests.tsv" % DS, sep="\t")
adj = pd.read_csv("results/partD_%s_depth_adjusted_tests.tsv" % DS, sep="\t")
pg = pd.read_csv("results/partD_%s_sgrg_gene_pseudobulk.tsv" % DS, sep="\t")
q5 = pd.read_csv("results/partD_%s_state_quintile.tsv" % DS, sep="\t")
ct_order = ["colon epithelial cell", "myeloid cell", "T cell of anorectum",
            "plasma cell", "stromal cell of lamina propria of colon"]
SHORT = {"colon epithelial cell": "Epithelial", "myeloid cell": "Myeloid",
         "T cell of anorectum": "T", "plasma cell": "Plasma",
         "stromal cell of lamina propria of colon": "Stromal"}


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig("results/%s.%s" % (name, ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[fig]", name)


# ---------- F1 细胞组成 ----------
fig, ax = plt.subplots(figsize=(88 * MM, 45 * MM))
ctab = cell.groupby(["disease", "cell_type"]).size().unstack(fill_value=0)
ctab = ctab.reindex(["normal", "ulcerative colitis", "Crohn disease"])
frac = ctab.div(ctab.sum(axis=1), axis=0) * 100
colors = ["#8DA0CB", "#FC8D62", "#66C2A5", "#E78AC3", "#A6D854"]
bottom = np.zeros(len(frac))
x = np.arange(len(frac))
for i, c in enumerate(ct_order):
    v = frac[c].values
    ax.bar(x, v, 0.62, bottom=bottom, color=colors[i], label=SHORT[c],
           edgecolor="white", linewidth=0.6)
    for xi, (b, vi) in enumerate(zip(bottom, v)):
        if vi > 6:
            ax.text(xi, b + vi / 2, "%.0f%%" % vi, ha="center", va="center",
                    fontsize=6, color="white")
    bottom += v
ax.set_xticks(x); ax.set_xticklabels([LAB[d] for d in frac.index])
ax.set_ylabel("Cell composition (%)")
ax.set_ylim(0, 100)
ax.legend(ncol=5, bbox_to_anchor=(0.5, 1.16), loc="upper center", fontsize=6,
          handlelength=1.0, columnspacing=1.0, handletextpad=0.4)
ax.text(0.99, -0.30, "n = 46,700 cells / 18 donors (GSE)", transform=ax.transAxes,
        ha="right", fontsize=6, color="#555555")
save(fig, "fig_D1_composition")

# ---------- F2 pseudobulk SGRG 评分 ----------
pb["SGRG_adj"] = np.nan
for c in pb.cell_type.unique():
    m = pb.cell_type == c
    b = np.polyfit(pb.loc[m, "log_nUMI"].values, pb.loc[m, "SGRG"].values, 1)
    pb.loc[m, "SGRG_adj"] = pb.loc[m, "SGRG"].values - np.polyval(b, pb.loc[m, "log_nUMI"].values)

fig, axes = plt.subplots(1, 2, figsize=(180 * MM, 48 * MM), sharey=True)
dis_order = ["normal", "ulcerative colitis", "Crohn disease"]
for ax, (col, ttl) in zip(axes, [("SGRG", "Observed panel score"),
                                 ("SGRG_adj", "Depth-adjusted (residual)")]):
    for i, c in enumerate(ct_order):
        for j, d in enumerate(dis_order):
            v = pb[(pb.cell_type == c) & (pb.disease == d)][col].values
            if len(v) == 0:
                continue
            xpos = i + (j - 1) * 0.24
            ax.scatter(np.full(len(v), xpos) + np.random.default_rng(j).normal(0, 0.018, len(v)),
                       v, s=7, color=C[d], alpha=0.9, linewidths=0)
            ax.plot([xpos - 0.07, xpos + 0.07], [np.mean(v)] * 2, color=C[d], lw=1.4)
    for i, c in enumerate(ct_order):
        t = tests[(tests.cell_type == c) & (tests.score == "SGRG") & (tests.contrast == "UC_vs_normal")]
        a2 = adj[adj.cell_type == c]
        if col == "SGRG":
            y = pb["SGRG"].max() * 1.02
            pv = t.iloc[0]["p_perm"] if len(t) else 1
        else:
            y = pb["SGRG_adj"].max() * 1.02
            pv = a2.iloc[0]["p_UC"] if len(a2) else 1
        if pv < 0.05:
            ax.text(i, y, "*", ha="center", fontsize=9)
    ax.set_xticks(range(len(ct_order)))
    ax.set_xticklabels([SHORT[c] for c in ct_order], rotation=25, ha="right")
    ax.set_title(ttl, fontsize=7)
axes[0].set_ylabel("SGRG panel score")
h = [plt.Line2D([], [], marker="o", ls="", color=C[d], markersize=4, label=LAB[d]) for d in dis_order]
axes[1].legend(handles=h, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.20), fontsize=6)
axes[0].text(-0.16, 1.13, "n = 6 donors / group", transform=axes[0].transAxes, fontsize=6, color="#555555")
save(fig, "fig_D2_pseudobulk_score")

# ---------- F3 深度混杂 ----------
fig, axes = plt.subplots(1, len(ct_order), figsize=(180 * MM, 42 * MM))
for ax, c in zip(axes, ct_order):
    s = cell[cell.cell_type == c]
    r, p = spearmanr(s["SGRG"], s["log_nUMI"])
    ax.hexbin(s["log_nUMI"], s["SGRG"], gridsize=28, cmap="Greys", mincnt=1, linewidths=0)
    idx = np.argsort(s["log_nUMI"].values)
    xx = s["log_nUMI"].values[idx]
    ax.plot(np.unique(xx), np.polyval(np.polyfit(s["log_nUMI"], s["SGRG"], 1), np.unique(xx)),
            color="#C44E52", lw=1.2)
    ax.set_title("%s\nrho = %+.2f" % (SHORT[c], r), fontsize=6.5)
    ax.set_xlabel("log10(UMI)")
axes[0].set_ylabel("SGRG panel score")
save(fig, "fig_D3_depth_confounding")

# ---------- F4 炎症五分位 ----------
fig, ax = plt.subplots(figsize=(85 * MM, 52 * MM))
cmap = plt.get_cmap("viridis")
for i, c in enumerate(ct_order):
    s = q5[q5.cell_type == c].sort_values("quintile")
    ax.errorbar(s["quintile"], s["sgrg_mean"], yerr=s["sgrg_sem"], marker="o", ms=3,
                lw=1.2, capsize=1.6, color=cmap(i / (len(ct_order) - 1) * 0.9), label=SHORT[c])
ax.axhline(0, color="#999999", lw=0.6, ls="--")
ax.set_xticks([1, 2, 3, 4, 5])
ax.set_xlabel("Inflammation-score quintile (depth-adjusted)")
ax.set_ylabel("Depth-adjusted SGRG score")
ax.legend(fontsize=6, ncol=2)
save(fig, "fig_D4_inflammation_quintile")

# ---------- F5 逐基因 Cohen's d 热图 ----------
dcol = "d_UC_vs_normal"
piv = pg.pivot(index="gene", columns="cell_type", values=dcol)
piv = piv[ct_order]
# 每个细胞类型内按 |d| 取前 12
sel = sorted(set(sum([[piv[c].abs().sort_values(ascending=False).index[i] for i in range(min(10, piv[c].notna().sum()))]
                      for c in ct_order], [])))
# 按最大 |d| 排序
mat = piv.loc[sel]
mat = mat.loc[mat.abs().max(axis=1).sort_values(ascending=False).index]
fig, ax = plt.subplots(figsize=(75 * MM, 100 * MM))
v = np.nanmax(np.abs(mat.values))
im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
ax.set_xticks(range(len(ct_order)))
ax.set_xticklabels([SHORT[c] for c in ct_order], rotation=30, ha="right")
ax.set_yticks(range(len(mat)))
ax.set_yticklabels(mat.index, fontsize=5.5)
ax.set_xticks(np.arange(-.5, len(ct_order), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(mat), 1), minor=True)
ax.grid(which="minor", color="white", linewidth=0.5)
ax.tick_params(which="minor", length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
cb.set_label("Cohen's d (UC vs Normal)", fontsize=6.5)
cb.ax.tick_params(labelsize=6)
ax.set_title("SGRG genes, pseudobulk (donor level)", fontsize=7)
save(fig, "fig_D5_gene_effect_heatmap")

print("\n完成 5 张图")
