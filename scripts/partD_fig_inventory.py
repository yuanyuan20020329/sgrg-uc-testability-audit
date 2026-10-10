# -*- coding: utf-8 -*-
"""图 AD4｜单细胞数据集级技术轴清单（A2 队列级宽表的等构对应物）

产出
----
  fig_AD4_technical_inventory.png/pdf

面板
----
  A 深度轴的两个口径（donor 级中位 vs 细胞级中位）—— 说明「UC 测得更深」不是稳健事实
  B 主轴占据度 × 深度-疾病共线（象限图）—— 技术轴「占据主轴」与「与分组共线」是两件事
  C ds4 的四道干预：名义 → 抽稀 → 仅自有对照 → 同 donor 炎症配对

注意：图内文字一律英文（图件用 Arial，中文会缺字形）。
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import NullFormatter

rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False,
})

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT = ROOT + "/results"

C_DONOR = "#4C72B0"
C_CELL = "#DD8452"
C_PRIM = "#55A868"
C_TECH = "#8172B2"
C_NOM = "#C44E52"

T = pd.read_csv(OUT + "/partD_technical_inventory.tsv", sep="\t").set_index("dataset")


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig("%s/%s.%s" % (OUT, name, ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[fig]", name)


fig = plt.figure(figsize=(7.1, 2.9))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.0], wspace=0.42,
                      left=0.075, right=0.985, top=0.86, bottom=0.17)

# ------------------------------------------------------------------ A 深度轴
ax = fig.add_subplot(gs[0, 0])
ds = list(T.index)
x = np.arange(len(ds))
w = 0.36
v1 = T.loc[ds, "depth_fold_uc_over_ctrl"].values          # donor 级
v2 = T.loc[ds, "depth_fold_celllevel"].values             # 细胞级
ax.bar(x - w / 2, v1, w, color=C_DONOR, label="donor-level median")
ax.bar(x + w / 2, v2, w, color=C_CELL, label="cell-level median")
for i, (a1, a2) in enumerate(zip(v1, v2)):
    ax.text(i - w / 2, a1 * 1.08, "%.2f" % a1, ha="center", fontsize=5.0, color=C_DONOR)
    ax.text(i + w / 2, a2 * 1.08, "%.2f" % a2, ha="center", fontsize=5.0, color=C_CELL)
for i, d in enumerate(ds):                                # ds4 第三口径
    vp = T.loc[d, "depth_fold_celllevel_primaryctrl"]
    if np.isfinite(vp):
        ax.bar(i + w / 2 + 0.20, vp, w * 0.5, color=C_PRIM)
        ax.text(i + w / 2 + 0.20, vp * 1.12, "%.2f" % vp, ha="center", va="bottom",
                fontsize=5.4, color=C_PRIM)
ax.axhline(1.0, color="#888", lw=0.7, ls="--")
ax.set_yscale("log")
ax.set_ylim(0.4, 9.0)
ax.set_yticks([0.5, 1, 2, 4, 8])
ax.set_yticklabels(["0.5", "1.0", "2.0", "4.0", "8.0"])
ax.yaxis.set_minor_formatter(NullFormatter())
ax.set_xticks(x)
ax.set_xticklabels(["%s\nAUC %.2f" % (d, T.loc[d, "depth_separability"]) for d in ds],
                   fontsize=6.0, linespacing=1.1)
ax.set_ylabel("median UMI per cell: UC / control")
ax.set_title("A  Depth axis is not a robust fact", fontsize=7.2, loc="left", pad=4)
ax.legend(fontsize=5.2, loc="upper left", handlelength=1.0, borderpad=0.2,
          labelspacing=0.25)
ax.text(3.42, 5.6, "green: cell-level,\nprimary ctrl only", ha="right", va="center",
        fontsize=5.0, color=C_PRIM, linespacing=1.15)

# ------------------------------------------------------------------ B 象限图
ax = fig.add_subplot(gs[0, 1])
for d in ds:
    r = T.loc[d]
    size = 14 + 210 * (r["n_ct_p05_perm"] / max(r["n_ct_tested"], 1))
    ax.scatter(r["depth_separability"], r["pc_tech_share"], s=size, color=C_TECH,
               alpha=0.75, lw=0.6, edgecolor="w", zorder=3)
    ax.annotate(d, (r["depth_separability"], r["pc_tech_share"]),
                textcoords="offset points", xytext=(0, 9), ha="center", fontsize=6.2)
    ax.plot([r["depth_separability"]] * 2,
            [r["pc_tech_share_cw"], r["pc_tech_share"]],
            color="#bbb", lw=0.7, zorder=1)
ax.axvline(0.5, color="#bbb", lw=0.6, ls=":")
ax.set_xlim(0.48, 1.0)
ax.set_ylim(0, 0.62)
ax.set_xlabel("depth-disease separability (AUC)")
ax.set_ylabel("main-axis technical share  PC1-3")
ax.set_title("B  Occupancy vs confounding", fontsize=7.2, loc="left", pad=4)
ax.text(0.985, 0.03, "marker size ~ nominal\nSGRG hits (p<0.05)\ngrey line: depth\n definition B",
        ha="right", va="bottom", fontsize=5.0, color="#666", linespacing=1.25)

# ------------------------------------------------------------------ C ds4 干预
ax = fig.add_subplot(gs[0, 2])
labels = ["naive\n(all controls)", "rarefied\n1,000 UMI", "primary\ncontrols only",
          "within-donor\ninflamed"]
num = [T.loc["ds4", "n_ct_p05_perm"], T.loc["ds4", "n_ct_p05_after_rarefy"],
       T.loc["ds4", "n_ct_p05_after_primary_ctrl"], T.loc["ds4", "n_ct_p05_after_within_donor"]]
den = [T.loc["ds4", "n_ct_tested"], 25, T.loc["ds4", "n_ct_in_primary_ctrl"],
       T.loc["ds4", "n_ct_in_within_donor"]]
frac = [n / d for n, d in zip(num, den)]
cols = [C_NOM, C_TECH, C_PRIM, C_DONOR]
yy = np.arange(len(labels))[::-1].astype(float)
ax.barh(yy, frac, 0.6, color=cols)
for y, n, d, f in zip(yy, num, den, frac):
    ax.text(f + 0.012, y, "%d/%d" % (n, d), va="center", fontsize=6.0)
ax.axvline(0.05, color="#888", lw=0.8, ls="--")
ax.text(0.058, yy.max() + 0.56, "expected 5%", fontsize=5.0, color="#666", va="center")
ax.set_yticks(yy)
ax.set_yticklabels(labels, fontsize=5.6, linespacing=1.05)
ax.set_ylim(-0.75, len(labels) - 0.35)
ax.set_xlim(0, 0.72)
ax.set_xlabel("cell types with SGRG p < 0.05")
ax.set_title("C  ds4: four interventions", fontsize=7.2, loc="left", pad=4)

save(fig, "fig_AD4_technical_inventory")
