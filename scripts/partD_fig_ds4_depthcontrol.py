# -*- coding: utf-8 -*-
"""图 D7｜ds4（Tearle，320,979 细胞）三种深度处理口径的对照
显示：观测 → 残差校正 → 深度匹配抽稀，效应量一致向 0 收缩，无一通过 FDR。
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False,
})
MM = 1 / 25.4
R = "results"

SHORT = {
    "CD4-positive, alpha-beta T cell": "CD4 T", "CD8-positive, alpha-beta T cell": "CD8 T",
    "CD8-positive, alpha-beta memory T cell": "CD8 T (mem)", "T follicular helper cell": "Tfh",
    "T-helper 17 cell": "Th17", "capillary endothelial cell": "Capillary EC",
    "conventional dendritic cell": "cDC", "endothelial cell of artery": "Artery EC",
    "epithelial cell of large intestine": "Epithelial", "fibroblast": "Fibroblast",
    "gamma-delta T cell": "gdT", "group 3 innate lymphoid cell": "ILC3",
    "intestinal crypt stem cell of large intestine": "Crypt stem",
    "intestine goblet cell": "Goblet", "macrophage": "Macrophage", "mast cell": "Mast",
    "memory B cell": "Memory B", "myofibroblast cell": "Myofibroblast",
    "naive B cell": "Naive B", "pericyte": "Pericyte", "plasma cell": "Plasma",
    "regulatory T cell": "Treg", "tuft cell of colon": "Tuft",
    "type I enteric glial cell": "Glial", "vein endothelial cell": "Vein EC",
}

obs = pd.read_csv(f"{R}/partD_ds4_tearle_uc_vs_hc_tests.tsv", sep="\t")
rar = pd.read_csv(f"{R}/partD_ds4_tearle_ds1k_uc_vs_hc_tests.tsv", sep="\t")


def get(df, key):
    d = df[df.score == key].copy()
    return d.set_index("cell_type")


O = get(obs, "SGRG")
A = get(obs, "SGRG_adj")
D = get(rar, "SGRG")

cts = [c for c in O.index if c in A.index and c in D.index]
step = (O.loc[cts, "diff"] - D.loc[cts, "diff"]).abs().sort_values(ascending=False)
cts = list(step.index)          # 收缩最大的排在上面
y = np.arange(len(cts))[::-1]

fig, axes = plt.subplots(1, 3, figsize=(180 * MM, max(58, 3.0 * len(cts) + 26) * MM),
                         sharey=True)
specs = [
    ("Observed", O, "raw", "#C44E52", "p_perm"),
    ("Residual-adjusted", A, "SGRG_adj", "#DD8452", "p_perm"),
    ("Depth-matched (rarefied)", D, "SGRG", "#4C72B0", "p_perm"),
]
lim = float(np.nanmax(np.abs(np.concatenate([O.loc[cts, "diff"].values,
                                             A.loc[cts, "diff"].values,
                                             D.loc[cts, "diff"].values])))) * 1.35
lim = max(lim, 0.06)

for ax, (ttl, df, key, col, pcol) in zip(axes, specs):
    ax.axvline(0, color="#888888", lw=0.8, ls=(0, (3, 2)), zorder=1)
    for i, c in enumerate(cts):
        v = float(df.loc[c, "diff"])
        p = float(df.loc[c, pcol])
        sig = p < 0.05
        ax.plot([0, v], [y[i], y[i]], color=col, lw=1.3, alpha=0.85,
                solid_capstyle="round", zorder=2)
        ax.scatter([v], [y[i]], s=26 if sig else 17, color=col,
                   edgecolor="black" if sig else "white",
                   linewidths=0.7 if sig else 0.5, zorder=3)
    ax.set_title(ttl, fontsize=7.2, pad=4)
    ax.set_xlabel("Δ SGRG score  (UC − normal)", fontsize=6.8)
    ax.tick_params(axis="x", labelsize=6)
    ax.set_xlim(-lim, lim)

axes[0].set_yticks(y)
axes[0].set_yticklabels([SHORT.get(c, c) for c in cts], fontsize=6)

fig.text(0.5, 0.008,
         "Filled dot + black outline = p < 0.05 (permutation test, donor-level pseudobulk).  "
         "ds4 = Tearle 2026, 8 UC vs 7 normal donors.",
         ha="center", fontsize=5.6, color="#555555")
fig.tight_layout(rect=(0, 0.045, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(f"{R}/fig_D7_ds4_depth_control.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("[fig] fig_D7_ds4_depth_control  (%d cell types)" % len(cts))
print(pd.DataFrame({"obs": O.loc[cts, "diff"], "obs_p": O.loc[cts, "p_perm"],
                    "adj": A.loc[cts, "diff"], "adj_p": A.loc[cts, "p_perm"],
                    "rar": D.loc[cts, "diff"], "rar_p": D.loc[cts, "p_perm"]}
                   ).to_string(float_format=lambda v: "%.4f" % v))
