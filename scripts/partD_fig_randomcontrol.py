# -*- coding: utf-8 -*-
"""Part D 图 D8：随机基因集负对照
——SGRG 的效应量/显著数是否落在「同等大小的随机基因集」零分布之内。
"""
import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 6.4,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
R = "results"
MM = 1 / 25.4

C_RAND = "#9aa3ae"
C_SGRG = "#c0392b"
C_POS = "#1f6fb4"


def rd(p):
    return pd.read_csv(p, sep="\t") if os.path.exists(p) else None


def jload(p, d=None):
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else d


def panel_hist(ax, vals, xlabel, sgrg=None, pos=None, bins=22, title=""):
    ax.hist(vals, bins=bins, color=C_RAND, edgecolor="white", linewidth=0.3)
    if sgrg is not None:
        ax.axvline(sgrg, color=C_SGRG, lw=1.2, zorder=5)
    for i, pv in enumerate(pos or []):
        ax.axvline(pv, color=C_POS, lw=0.9, ls="--", zorder=4, alpha=0.85 + 0.15 * (i == 0))
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Random gene sets (n)", fontsize=6.0)
    if title:
        ax.set_title(title, fontsize=6.4, pad=3)


def main():
    S = jload(f"{R}/partD_ds4_tearle_rand_summary.json")
    if S is None:
        print("缺 ds4_tearle_rand_summary.json"); return
    T = rd(f"{R}/partD_ds4_tearle_rand_geneset_tests.tsv")
    rnd = T[T.geneset.str.startswith("RAND_")].copy()

    # 每个随机集的汇总
    g = rnd.groupby("geneset")
    rs = pd.DataFrame({
        "n_p05": g["p"].apply(lambda v: (v < 0.05).sum()),
        "n_q05": g["q_BH"].apply(lambda v: (v < 0.05).sum()),
        "n_p05_perm": g["p_perm"].apply(lambda v: (v < 0.05).sum()),
        "min_p": g["p"].min(),
        "median_abs_d": g["diff"].apply(lambda v: np.median(np.abs(v))),
        "n_up": g["diff"].apply(lambda v: (v > 0).sum()),
    })
    real = pd.DataFrame(S["real"]).set_index("geneset")
    pos = pd.DataFrame(S["positive"]).set_index("geneset")
    n_ct = int(S["n_ct_tested"])
    exp_fp = n_ct * 0.05

    fig = plt.figure(figsize=(180 * MM, 118 * MM))
    gs = fig.add_gridspec(2, 3, hspace=0.52, wspace=0.34)

    # (A) Mann-Whitney 口径
    ax = fig.add_subplot(gs[0, 0])
    panel_hist(ax, rs.n_p05.values,
               "Cell types with $p$ < 0.05  (of %d)" % n_ct,
               sgrg=float(real.loc["SGRG", "n_p05"]),
               pos=[float(pos.loc[i, "n_p05"]) for i in pos.index],
               bins=np.arange(-0.5, max(rs.n_p05.max(), float(real.loc["SGRG", "n_p05"])) + 1.5, 1),
               title="(A)  Significant cell types — Mann–Whitney $U$")
    ax.axvline(exp_fp, color="#1a7f37", lw=0.9, ls="-.")
    ax.text(0.04, 0.96, "median %d = %.1f× null" % (
        np.median(rs.n_p05), np.median(rs.n_p05) / exp_fp),
        transform=ax.transAxes, va="top", fontsize=5.2)

    # (B) 置换检验口径
    ax = fig.add_subplot(gs[0, 1])
    panel_hist(ax, rs.n_p05_perm.values,
               "Cell types with $p_{\\rm perm}$ < 0.05  (of %d)" % n_ct,
               sgrg=float(real.loc["SGRG", "n_p05_perm"]),
               pos=[float(pos.loc[i, "n_p05_perm"]) for i in pos.index],
               bins=np.arange(-0.5, max(rs.n_p05_perm.max(),
                                        float(real.loc["SGRG", "n_p05_perm"])) + 1.5, 1),
               title="(B)  Significant cell types — permutation test")
    ax.axvline(exp_fp, color="#1a7f37", lw=0.9, ls="-.")
    ax.text(0.04, 0.96, "median %d = %.1f× null" % (
        np.median(rs.n_p05_perm), np.median(rs.n_p05_perm) / exp_fp),
        transform=ax.transAxes, va="top", fontsize=5.2)

    # (C) 效应量
    ax = fig.add_subplot(gs[0, 2])
    panel_hist(ax, rs.median_abs_d.values,
               "Median |$\\Delta$| effect size",
               sgrg=float(real.loc["SGRG", "median_abs_d"]),
               pos=[float(pos.loc[i, "median_abs_d"]) for i in pos.index],
               title="(C)  Effect magnitude")
    ax.text(0.04, 0.96, "SGRG at %dth pct of random" % int(
        100 * np.mean(rs.median_abs_d.values <=
                      float(real.loc["SGRG", "median_abs_d"]))),
        transform=ax.transAxes, va="top", fontsize=5.2)

    # (D) 逐细胞类型：随机集 Δ 的分布 vs SGRG Δ
    ax = fig.add_subplot(gs[1, :])
    cts = S["cell_types_tested"]
    sgrg_diff = T[T.geneset == "SGRG"].set_index("cell_type")["diff"]
    order = sorted(cts, key=lambda c: -abs(float(sgrg_diff.get(c, 0.0))))
    for i, c in enumerate(order):
        v = rnd[rnd.cell_type == c]["diff"].values
        if len(v) == 0:
            continue
        lo, hi = np.percentile(v, [5, 95])
        ax.plot([lo, hi], [i, i], color=C_RAND, lw=1.6, solid_capstyle="round")
        ax.plot([np.median(v)], [i], "o", color=C_RAND, ms=2.4, zorder=3)
        ax.plot([float(sgrg_diff.get(c, 0.0))], [i], "D", color=C_SGRG, ms=3.0, zorder=5)
    ax.axvline(0, color="black", lw=0.5, ls=":")
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([c if len(c) <= 40 else c[:38] + "…" for c in order], fontsize=5.3)
    ax.set_ylim(-1, len(order))
    ax.set_xlabel("UC $-$ control difference in panel score (donor-level pseudobulk)")
    ax.set_title("(D)  Per-cell-type effect: SGRG (red diamond) vs random gene sets "
                 "(grey = median and 5–95%% range of %d random sets)" % len(rs), fontsize=6.2, pad=3)
    ax.legend(handles=[Line2D([], [], color=C_RAND, lw=1.6, label="Random gene sets (n=%d)" % len(rs)),
                       Line2D([], [], color=C_SGRG, marker="D", ls="", ms=3.0, label="SGRG panel")],
              loc="lower right", fontsize=5.4, frameon=False, handlelength=1.6)

    for ax_ in fig.axes:
        for s in ("top", "right"):
            ax_.spines[s].set_visible(False)

    fig.text(0.5, -0.012,
             "Random gene sets are expression-stratified matched controls (identical size, identical "
             "mean-expression bin composition, SGRG genes excluded) and pass through exactly the same pipeline.",
             ha="center", fontsize=5.4, color="#555555")
    fig.savefig(f"{R}/fig_D8_random_control.png", dpi=300)
    fig.savefig(f"{R}/fig_D8_random_control.pdf")
    print("wrote fig_D8_random_control.{png,pdf}")

    # 打印关键数字
    print("\n随机集 n_p05: 中位 %.0f  5-95%% [%.0f, %.0f]  max %d" % (
        np.median(rs.n_p05), np.percentile(rs.n_p05, 5), np.percentile(rs.n_p05, 95), rs.n_p05.max()))
    print("SGRG n_p05 = %s ; 阳性对照 = %s" % (
        real.loc["SGRG", "n_p05"], {i: pos.loc[i, "n_p05"] for i in pos.index}))


if __name__ == "__main__":
    main()
