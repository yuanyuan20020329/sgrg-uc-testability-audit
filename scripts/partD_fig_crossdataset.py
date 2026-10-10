# -*- coding: utf-8 -*-
"""图 D6｜跨数据集效应量森林图：原始面板评分 vs 深度校正后（UC − 正常/对照）
这是 Part D 的核心结论图：名义效应在深度校正后一致向 0 收缩。
"""
import os
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
    "T cell of anorectum": "T cell", "colon epithelial cell": "Epithelial",
    "myeloid cell": "Myeloid", "plasma cell": "Plasma",
    "stromal cell of lamina propria of colon": "Stromal",
    "monocyte": "Monocyte", "colonocyte": "Colonocyte",
    "intestine goblet cell": "Goblet", "transit amplifying cell": "TA cell",
    "B cell": "B cell", "T cell": "T cell", "epithelial cell": "Epithelial",
    "macrophage": "Macrophage", "endothelial cell": "Endothelial",
    "fibroblast": "Fibroblast", "NK cell": "NK cell", "mast cell": "Mast cell",
}
DSCOL = {"ds1": "#4C72B0", "ds2": "#DD8452", "ds3": "#55A868", "ds4": "#C44E52"}
DSLAB = {"ds1": "ds1 (discovery)", "ds2": "ds2 myeloid", "ds3": "ds3 epithelium",
         "ds4": "ds4 Tearle (validation)"}


def rd(p):
    return pd.read_csv(p, sep="\t") if os.path.exists(p) else None


rows = []      # (ds, cell_type, raw_diff, adj_diff, p_raw, p_adj, n_uc, n_ctrl)

# ds1：raw 用 pseudobulk_tests（UC_vs_normal），adj 用 depth_adjusted_tests
t1 = rd(f"{R}/partD_ds1_pseudobulk_tests.tsv")
a1 = rd(f"{R}/partD_ds1_depth_adjusted_tests.tsv")
if t1 is not None and a1 is not None:
    t1u = t1[(t1.score == "SGRG") & (t1.contrast == "UC_vs_normal")]
    for _, r in t1u.iterrows():
        ar = a1[a1.cell_type == r["cell_type"]]
        rows.append(("ds1", r["cell_type"], r["diff"], float(ar["diff_UC"].iloc[0]) if len(ar) else np.nan,
                     r["p_perm"], float(ar["p_UC"].iloc[0]) if len(ar) else np.nan,
                     int(r["n1"]), int(r["n2"])))

# ds2/ds3/ds4：同一张 tests 表内含 SGRG 与 SGRG_adj
for ds, f in [("ds2", "partD_ds2_myeloid_uc_vs_hc_tests.tsv"),
              ("ds3", "partD_ds3_largeintestine_uc_vs_hc_tests.tsv"),
              ("ds4", "partD_ds4_tearle_uc_vs_hc_tests.tsv")]:
    t = rd(f"{R}/{f}")
    if t is None:
        continue
    # ds4 的「深度控制」用抽稀到 1000 UMI 的独立运行（而非残差法），避免共线性伪影
    tctl = rd(f"{R}/partD_ds4_tearle_ds1k_uc_vs_hc_tests.tsv") if ds == "ds4" else None
    for c in t.cell_type.unique():
        s = t[(t.cell_type == c) & (t.score == "SGRG")]
        if tctl is not None:
            j = tctl[(tctl.cell_type == c) & (tctl.score == "SGRG")]
        else:
            j = t[(t.cell_type == c) & (t.score == "SGRG_adj")]
        if not len(s) or not len(j):
            continue
        rows.append((ds, c, float(s["diff"].iloc[0]), float(j["diff"].iloc[0]),
                     float(s["p_perm"].iloc[0]), float(j["p_perm"].iloc[0]),
                     int(s["n_uc"].iloc[0]), int(s["n_hc"].iloc[0])))

df = pd.DataFrame(rows, columns=["ds", "cell_type", "raw", "adj", "p_raw", "p_adj", "n_uc", "n_ctrl"])
df["lab"] = df.cell_type.map(lambda c: SHORT.get(c, c))
df["y"] = np.arange(len(df))[::-1]

fig, axes = plt.subplots(1, 2, figsize=(180 * MM, max(48, 5.2 * len(df) + 22) * MM), sharey=True)
titles = ["Observed panel score", "Depth-controlled"]
cols = ["raw", "adj"]
pcols = ["p_raw", "p_adj"]

for ax, c, pc, ttl in zip(axes, cols, pcols, titles):
    ax.axvline(0, color="#888888", lw=0.8, ls=(0, (3, 2)), zorder=1)
    for _, r in df.iterrows():
        col = DSCOL[r["ds"]]
        sig = r[pc] < 0.05
        ax.plot([0, r[c]], [r["y"], r["y"]], color=col, lw=1.3, alpha=0.85,
                solid_capstyle="round", zorder=2)
        ax.scatter([r[c]], [r["y"]], s=26 if sig else 18, color=col,
                   edgecolor="black" if sig else "white", linewidths=0.7 if sig else 0.5,
                   zorder=3)
        ax.text(r[c] + (0.0035 if r[c] >= 0 else -0.0035), r["y"],
                "n=%d/%d" % (r["n_uc"], r["n_ctrl"]), va="center",
                ha="left" if r[c] >= 0 else "right", fontsize=5.2, color="#555555")
    ax.set_title(ttl, fontsize=7.2, pad=4)
    ax.set_xlabel("Δ SGRG score  (UC − control)", fontsize=6.8)
    ax.tick_params(axis="x", labelsize=6)
    lim = max(0.035, float(np.nanmax(np.abs(df[["raw", "adj"]].values))) * 1.45)
    ax.set_xlim(-lim, lim)

axes[0].set_yticks(df["y"].values)
axes[0].set_yticklabels(["%s · %s" % (r["lab"], r["ds"]) for _, r in df.iterrows()], fontsize=6)

handles = [plt.Line2D([], [], color=DSCOL[d], marker="o", ls="-", lw=1.3, ms=4, label=DSLAB[d])
           for d in ["ds1", "ds2", "ds3", "ds4"] if d in set(df["ds"])]
axes[1].legend(handles=handles, fontsize=6, loc="lower right")
fig.text(0.5, 0.005,
         "Filled dot + black outline = p < 0.05 (permutation, donor-level pseudobulk);  n = UC / control donors.  "
         "Depth-controlled: ds1/ds2/ds3 = regression residual on log10(nUMI); ds4 = rarefaction to 1000 UMI per cell.",
         ha="center", fontsize=5.6, color="#555555")
fig.tight_layout(rect=(0, 0.035, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(f"{R}/fig_D6_cross_dataset_forest.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("[fig] fig_D6_cross_dataset_forest")
print(df[["ds", "lab", "raw", "adj", "p_raw", "p_adj", "n_uc", "n_ctrl"]].to_string(index=False))
