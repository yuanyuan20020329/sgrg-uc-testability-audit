# -*- coding: utf-8 -*-
"""图 D10：对照组的「来源组成」如何决定了结论
(A) ds4 细胞级深度分布：UC / 本研究对照 / 外部复用「正常」
(B) ds4 donor 级中位深度（按来源着色）
(C) 四套数据的「对照同源性」：对照细胞中来自含病例批次/研究的比例
(D) SGRG 效应量：朴素口径 vs 同源口径（若结果存在）
"""
import json
import os

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

MM = 1 / 25.4
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 6.2,
    "axes.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
R = "results"
C_UC = "#d1495b"
C_CTRLP = "#2e6f9e"
C_EXT = "#9aa5b1"
C_SGRG = "#b5442e"
C_RAND = "#8d99ae"


def rd(p, **kw):
    p = os.path.join(R, p) if not os.path.isabs(p) else p
    return pd.read_csv(p, sep="\t", **kw) if os.path.exists(p) else None


# ---------------------------------------------------------------- ds4 细胞级元数据
def obs_col(g, col):
    d = g[col]
    if isinstance(d, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x)
                         for x in d["categories"][:]], dtype=object)
        c = d["codes"][:]
        o = np.empty(len(c), dtype=object); ok = c >= 0
        o[ok] = cats[c[ok]]; o[~ok] = ""
        return o
    arr = d[:]
    if arr.dtype.kind in "OS":
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr], dtype=object)
    return arr


f = h5py.File("data/raw/sc/CxG_Tearle2026_PSC_UC.h5ad", "r")
obs = f["obs"]
d4 = pd.DataFrame({"disease": obs_col(obs, "disease"),
                   "donor": obs_col(obs, "donor_id"),
                   "primary": np.asarray(obs_col(obs, "is_primary_data")).astype(bool),
                   "batch": obs_col(obs, "batch_condition")})
d4["nUMI"] = np.load("data/raw/sc/ds4_cell_nUMI.npy")
f.close()

d4["grp4"] = np.where(d4.disease.str.contains("ulcer"), "UC (primary)",
            np.where(d4.disease.str.contains("primary sclerosing"), "PSC (primary)",
            np.where(d4.primary, "Control (primary)", "Control (reused, external)")))
CATC = {"UC (primary)": C_UC, "PSC (primary)": "#7a5195",
        "Control (primary)": C_CTRLP, "Control (reused, external)": C_EXT}

# ---------------------------------------------------------------- 图
fig = plt.figure(figsize=(180 * MM, 104 * MM))
gs = GridSpec(2, 3, figure=fig, hspace=0.62, wspace=0.46,
              left=0.195, right=0.985, top=0.90, bottom=0.115)

# (A) 细胞级深度分布
ax = fig.add_subplot(gs[0, 0])
order = ["UC (primary)", "PSC (primary)", "Control (primary)", "Control (reused, external)"]
for i, k in enumerate(order):
    v = np.log10(d4.loc[d4.grp4 == k, "nUMI"].clip(lower=1).values)
    if not len(v):
        continue
    parts = ax.violinplot(v, positions=[i], widths=0.82, showextrema=False)
    for b in parts["bodies"]:
        b.set_facecolor(CATC[k]); b.set_alpha(0.75); b.set_edgecolor("none")
    ax.plot([i], [np.median(v)], "o", color="white", ms=2.0, zorder=5)
    ax.plot([i], [np.median(v)], "o", mfc="none", mec="black", ms=3.6, mew=0.6, zorder=6)
ax.set_xticks(range(len(order)))
ax.set_xticklabels(["UC\n(primary)", "PSC\n(primary)", "Control\n(primary)",
                    "Control\n(reused)"], fontsize=5.4)
ax.set_ylabel("log$_{10}$(UMI per cell)", fontsize=6.0)
ax.set_title("(A)  ds4 cell-level depth", fontsize=6.2, pad=3)
for i, k in enumerate(order):
    n = int((d4.grp4 == k).sum())
    ax.text(i, ax.get_ylim()[1] * 0.995, "n=%s" % (("%.1fk" % (n / 1000)) if n >= 1000 else n),
            ha="center", va="top", fontsize=4.8, color="#444444")

# (B) donor 级中位深度
ax = fig.add_subplot(gs[0, 1])
g = d4.groupby(["grp4", "donor", "disease", "primary"]).nUMI.median().reset_index()
g = g.sort_values("nUMI").reset_index(drop=True)
for i, r in g.iterrows():
    ax.scatter(i, r.nUMI, s=13, color=CATC[r.grp4], zorder=4,
               edgecolor="black", linewidth=0.35)
ax.set_yscale("log")
ax.set_ylabel("median UMI per donor", fontsize=6.0)
ax.set_xticks([])
ax.set_xlabel("donors (sorted)", fontsize=6.0)
ax.set_title("(B)  ds4 donor-level depth", fontsize=6.2, pad=3)
ucm = float(np.median(d4.loc[d4.grp4 == "UC (primary)", "nUMI"]))
hcp = float(np.median(d4.loc[d4.grp4 == "Control (primary)", "nUMI"]))
hcx = float(np.median(d4.loc[d4.grp4 == "Control (reused, external)", "nUMI"]))
for y, lab, col, xa in [(ucm, "UC %.0f" % ucm, C_UC, 0.02),
                        (hcp, "own ctrl %.0f" % hcp, C_CTRLP, 0.02),
                        (hcx, "reused %.0f" % hcx, C_EXT, 0.02)]:
    ax.axhline(y, color=col, lw=0.6, ls="--", zorder=1)
    ax.text(xa, y, " " + lab, transform=ax.get_yaxis_transform(),
            fontsize=4.6, color=col, va="bottom", bbox=dict(fc="white", ec="none", pad=0.4))
ax.text(0.97, 0.06, "own-ctrl fold = %.2f×\nnaive fold = %.2f×" % (ucm / hcp, ucm / hcx),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=5.0, color="#333333")

# (C) 四套数据的对照同源性
ax = fig.add_subplot(gs[0, 2])
S = rd("partD_batch_audit_summary.tsv")
names, vals, labs = [], [], []
order_ds = ["ds1_garrido", "ds2_myeloid", "ds3_largeintestine", "ds4_tearle"]
LAB = {"ds1_garrido": "ds1\nGarrido\nTrigo", "ds2_myeloid": "ds2\nOliver\nmyeloid",
       "ds3_largeintestine": "ds3\nOliver\nlarge int.", "ds4_tearle": "ds4\nTearle\nPSC/UC"}
which = {"ds1_garrido": "frac_ctrl_cells_samebatch_assay",
         "ds2_myeloid": "frac_ctrl_cells_samebatch_study",
         "ds3_largeintestine": "frac_ctrl_cells_samebatch_study",
         "ds4_tearle": "frac_ctrl_cells_samebatch_batch_condition"}
for ds in order_ds:
    r = S[S.dataset == ds]
    if not len(r):
        continue
    v = float(r[which[ds]].iloc[0])
    names.append(LAB[ds]); vals.append(v)
    labs.append("same %s" % which[ds].split("_")[-1].replace("batch_condition", "batch"))
bars = ax.bar(range(len(vals)), vals, color=["#4c9f70" if v > 0.6 else "#c1443c" for v in vals],
              edgecolor="black", linewidth=0.4, width=0.66)
for i, (v, l) in enumerate(zip(vals, labs)):
    ax.text(i, v + 0.03, "%.0f%%" % (100 * v), ha="center", fontsize=5.2)
    ax.text(i, 0.02, l, ha="center", fontsize=4.2, color="white", rotation=90, va="bottom")
ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, fontsize=5.0)
ax.set_ylim(0, 1.16)
ax.set_ylabel("control cells from a\ncase-containing batch", fontsize=5.6)
ax.set_title("(C)  Control-group homology", fontsize=6.2, pad=3)
ax.axhline(1.0, color="#4c9f70", lw=0.6, ls=":")
ax.axhline(0.5, color="#888888", lw=0.6, ls=":")

# (D) SGRG 效应量：朴素 vs 同源
ax = fig.add_subplot(gs[1, :])
T = rd("partD_ds4_prov_contrast_tests.tsv")
if T is not None:
    a = T[(T.contrast == "A_naive_all_normal") & (T.geneset == "SGRG")].set_index("cell_type")
    b = T[(T.contrast == "B_primary_only") & (T.geneset == "SGRG")].set_index("cell_type")
    ct = sorted(set(a.index) | set(b.index),
                key=lambda c: -(abs(b["diff"].get(c, 0)) if c in b.index else abs(a["diff"].get(c, 0))))
    ct = ct[:20][::-1]
    y = np.arange(len(ct))
    for i, c in enumerate(ct):
        va = a.loc[c, "diff"] if c in a.index else np.nan
        vb = b.loc[c, "diff"] if c in b.index else np.nan
        if np.isfinite(va):
            ax.plot([0, va], [i + 0.16, i + 0.16], color=C_EXT, lw=0.7, zorder=2)
            ax.plot([va], [i + 0.16], "o", color=C_EXT, ms=3.4, zorder=3)
        if np.isfinite(vb):
            ax.plot([0, vb], [i - 0.16, i - 0.16], color=C_SGRG, lw=0.7, zorder=2)
            ax.plot([vb], [i - 0.16], "D", color=C_SGRG, ms=2.8, zorder=3)
    ax.axvline(0, color="black", lw=0.5, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels([c[:42] for c in ct], fontsize=4.5)
    ax.set_xlabel("$\\Delta$ SGRG score  (case − control)", fontsize=6.0)
    ax.set_title("(D)  Same data, two control definitions:  grey = UC vs all 'normal' (raw, 94% reused);  "
                 "red = UC vs the study's own controls (n=2 donors)", fontsize=5.8, pad=3)
    ax.legend(handles=[plt.Line2D([], [], color=C_EXT, marker="o", ls="-", ms=3.4,
                                  lw=0.7, label="naive (all 'normal')"),
                       plt.Line2D([], [], color=C_SGRG, marker="D", ls="-", ms=2.8,
                                  lw=0.7, label="own controls only")],
              fontsize=5.0, loc="lower right", frameon=False)
else:
    ax.text(0.5, 0.5, "contrast results not found", ha="center", va="center",
            transform=ax.transAxes, fontsize=6)
    ax.set_axis_off()

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(R, "fig_D10_provenance_confounding.%s" % ext), dpi=600)
print("saved fig_D10_provenance_confounding.png / .pdf")
