# -*- coding: utf-8 -*-
"""图 D11：三种口径的检验力与校准诊断（Part D 的最终裁决）

(A) 三个口径下 SGRG 的「显著细胞类型数」对照随机基因集零分布
(B) 阳性对照诊断：A 口径有 2–7 个炎症通路显著，B 口径 0 个、C 口径 0 个（但方向 13–15/15 全向上）
    —— B/C 不是「阴性」，是「没有检验力」
(C) 配对设计下最小可达双侧 p 与配对 donor 数的关系（标注本研究只有 4 对）
(D) 同 donor 炎症配对：15 种细胞类型的 SGRG 配对差（含每个 donor 的点）
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import binom

MM = 1 / 25.4
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 6.2,
    "axes.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
R = "results"
C_SGRG = "#b5442e"
C_RAND = "#8d99ae"
C_NAIVE = "#9aa5b1"
C_PRIM = "#2e6f9e"
C_INFL = "#4c9f70"
POSE = ["H_TNFA_SIGNALING_VIA_NFKB", "H_INFLAMMATORY_RESPONSE",
        "H_INTERFERON_GAMMA_RESPONSE", "H_HYPOXIA", "H_UNFOLDED_PROTEIN_RESPONSE"]
SHORT = {"H_TNFA_SIGNALING_VIA_NFKB": "TNFA/NFKB", "H_INFLAMMATORY_RESPONSE": "Inflam",
         "H_INTERFERON_GAMMA_RESPONSE": "IFN-γ", "H_HYPOXIA": "Hypoxia",
         "H_UNFOLDED_PROTEIN_RESPONSE": "UPR"}

S = pd.read_csv(os.path.join(R, "partD_ds4_prov_geneset_summary.tsv"), sep="\t")
T = pd.read_csv(os.path.join(R, "partD_ds4_prov_contrast_tests.tsv"), sep="\t")
C = pd.read_csv(os.path.join(R, "partD_ds4_prov_inflamed_tests.tsv"), sep="\t")

fig = plt.figure(figsize=(180 * MM, 116 * MM))
gs = GridSpec(2, 3, figure=fig, hspace=0.60, wspace=0.40,
              left=0.135, right=0.985, top=0.905, bottom=0.10)

CONTRASTS = ["A_naive", "B_primary", "C_inflamed"]
CLAB = ["A  naive\nUC vs all\n'normal'", "B  primary\nUC vs own\ncontrols (n=2)",
        "C  within-donor\ninflamed vs\nnon-infl (4 pairs)"]
CCOL = [C_NAIVE, C_PRIM, C_INFL]

# ---------------------------------------------------------------- (A)
ax = fig.add_subplot(gs[0, 0])
x = np.arange(len(CONTRASTS))
for i, cname in enumerate(CONTRASTS):
    rs = S[(S.contrast == cname) & (S.kind == "random")]
    rv = S[(S.contrast == cname) & (S.geneset == "SGRG")]
    if not len(rs):
        continue
    v = rs.n_p05.values
    lo, hi = np.percentile(v, [10, 90]); m = np.median(v)
    ax.plot([i, i], [lo, hi], color=C_RAND, lw=3.0, solid_capstyle="round", zorder=2)
    ax.plot([i], [m], "_", color="black", ms=11, mew=1.1, zorder=4)
    nc = int(rs.n_ct.iloc[0]) if len(rs) else 0
    ax.plot([i], [nc * 0.05], "v", color="#1a7f37", ms=3.6, zorder=4)
    if len(rv):
        ax.plot([i], [float(rv.n_p05.iloc[0])], "D", color=C_SGRG, ms=5.0, zorder=5)
ax.axhline(1.0, color="#1a7f37", lw=0.6, ls=":")
ax.set_xticks(x); ax.set_xticklabels(CLAB, fontsize=4.7)
ax.set_ylabel("# cell types with $p$ < 0.05", fontsize=6.0)
ax.set_title("(A)  'Significance' count: SGRG vs random-set null", fontsize=6.2, pad=3)
ax.set_ylim(-0.6, max(16, ax.get_ylim()[1]))
ax.legend(handles=[plt.Line2D([], [], marker="D", color=C_SGRG, ls="", ms=4.6, label="SGRG panel"),
                   plt.Line2D([], [], color=C_RAND, lw=3, label="200 random sets\n(median, 10–90%)"),
                   plt.Line2D([], [], marker="v", color="#1a7f37", ls="", ms=3.6,
                              label="expected by chance\n(5% of cell types)")],
          fontsize=4.4, loc="upper right", frameon=False, handlelength=1.4)

# ---------------------------------------------------------------- (B)
ax = fig.add_subplot(gs[0, 1])
w = 0.26
for i, cname in enumerate(CONTRASTS):
    d = S[(S.contrast == cname) & (S.kind == "positive")]
    for j, g in enumerate(POSE):
        row = d[d.geneset == g]
        if not len(row):
            continue
        ax.bar(i + (j - 2) * w, float(row.n_p05.iloc[0]), width=w * 0.9,
               color=plt.cm.YlGnBu((j + 1) / 6.0), edgecolor="black", linewidth=0.25, zorder=3)
ax.set_xticks(x); ax.set_xticklabels(CLAB, fontsize=4.7)
ax.set_ylabel("# cell types with $p$ < 0.05", fontsize=6.0)
ax.set_title("(B)  Positive controls: are the assays\n      able to detect known signals at all?",
             fontsize=6.2, pad=3)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=plt.cm.YlGnBu((j + 1) / 6.0),
                                 ec="black", lw=0.25, label=SHORT[g])
                   for j, g in enumerate(POSE)],
          fontsize=4.4, loc="upper right", frameon=False, ncol=1, handlelength=0.9,
          handleheight=0.7, labelspacing=0.25)
ax.text(0.5, 0.60, "B and C cannot even detect\ninflammation pathways →\nzero power, not zero effect",
        transform=ax.transAxes, ha="center", va="center", fontsize=4.7,
        color="#8a1c1c", style="italic")

# ---------------------------------------------------------------- (C)
ax = fig.add_subplot(gs[0, 2])
ns = np.arange(2, 11)
minp = [2.0 / (2 ** n) for n in ns]      # 双侧符号翻转检验的最小可达 p
ax.plot(ns, -np.log10(minp), "-o", color="#444444", ms=3.0, lw=0.9)
ax.axhline(-np.log10(0.05), color="#c1443c", lw=0.7, ls="--")
ax.text(2.1, -np.log10(0.05) + 0.08, "$p=0.05$", fontsize=4.8, color="#c1443c")
ax.axvline(4, color=C_INFL, lw=0.8, ls=":")
ax.text(4.15, 0.15, "this study: 4 paired donors\n→ min $p$ = 0.125", fontsize=4.7, color=C_INFL)
ax.set_xlabel("# paired donors", fontsize=6.0)
ax.set_ylabel("$-\\log_{10}$(smallest attainable $p$)", fontsize=5.6)
ax.set_title("(C)  Power ceiling of a paired design", fontsize=6.2, pad=3)
ax.set_xticks(ns)

# ---------------------------------------------------------------- (D)
ax = fig.add_subplot(gs[1, :])
s = C[C.geneset == "SGRG"].copy()
s = s.sort_values("mean_diff")
y = np.arange(len(s))
for i, (_, r) in enumerate(s.iterrows()):
    ax.plot([0, r.mean_diff], [i, i], color=C_INFL, lw=0.7, alpha=0.8, zorder=2)
    ax.plot([r.mean_diff], [i], "D", color=C_SGRG, ms=3.6, zorder=4)
ax.axvline(0, color="black", lw=0.5, ls=":")
ax.set_yticks(y)
ax.set_yticklabels(["%s  (n=%d)" % (c[:38], int(s.iloc[i].n_donors)) for i, c in enumerate(s.cell_type)],
                   fontsize=4.6)
for i, (_, r) in enumerate(s.iterrows()):
    ax.text(ax.get_xlim()[1], i, "", fontsize=4)
ax.set_xlabel("$\\Delta$ SGRG score  (inflamed − non-inflamed, paired within donor)", fontsize=6.0)
ax.set_title("(D)  The only confound-free design available in ds4 — and it has 4 donors",
             fontsize=6.2, pad=3)
ax.text(0.985, 0.06,
        "green dots in (B):  direction of inflammation pathways in C = 13–15 of 15 cell types up\n"
        "SGRG direction in C = 11 of 15 up, but no cell type reaches the attainable floor (min p = 0.125)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=4.7, color="#444444")

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(R, "fig_D11_power_diagnosis.%s" % ext), dpi=600)
print("saved fig_D11_power_diagnosis.png / .pdf")
