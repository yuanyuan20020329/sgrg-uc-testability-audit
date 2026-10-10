# -*- coding: utf-8 -*-
"""A2 + Part D 合并稿图件：跨分辨率剂量-反应 & 干预-响应

产出
----
  fig_AD1_dose_response.png/pdf          跨分辨率剂量-反应（并排，不合并）
  fig_AD2_intervention_response.png/pdf  干预-响应 2×2（方法学主论点）

注意：图内文字一律英文（本项目图件用 Arial，中文会缺字形）。
"""
import json
import os

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

C_BULK = "#4C72B0"
C_SC = "#C44E52"
C_TECH = "#8172B2"
C_NOM = "#DD8452"


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig("%s/%s.%s" % (OUT, name, ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[fig]", name)


def sp(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 6:
        return np.nan, np.nan, int(m.sum())
    r, p = spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def pfmt(p):
    if not np.isfinite(p):
        return "n.s."
    return "p = %.1e" % p if p < 0.01 else "p = %.3f" % p


# ==================================================================== 读数据
df = pd.read_csv(OUT + "/partA2_cohorts.tsv", sep="\t")
df["d_bio"] = pd.to_numeric(df["auc_pc1_batch"], errors="coerce") - \
              pd.to_numeric(df["auc_pc1_primary"], errors="coerce")
coll = pd.read_csv(OUT + "/partD_ds4_depth_bio_collinearity.tsv", sep="\t")

# 深度比：来自原流水线检验表（描述性统计，与框架无关）
raw_t = pd.read_csv(OUT + "/partD_ds4_tearle_uc_vs_hc_tests.tsv", sep="\t")
ds1_t = pd.read_csv(OUT + "/partD_ds4_tearle_ds1k_uc_vs_hc_tests.tsv", sep="\t")

# 显著性：一律用「随机对照框架」内的 SGRG（与 200 个随机集同管线、同基线，可与零分布直接比）
raw_sg = pd.read_csv(OUT + "/partD_ds4_tearle_rand_geneset_tests.tsv", sep="\t")
raw_sg = raw_sg[raw_sg.geneset == "SGRG"]
ds1_sg = pd.read_csv(OUT + "/partD_ds4_tearle_rand_ds1k_geneset_tests.tsv", sep="\t")
ds1_sg = ds1_sg[ds1_sg.geneset == "SGRG"]


def n_sig(t, score="SGRG"):
    """(置换 p<0.05, FDR q<0.05, 检验的细胞类型数) —— 随机对照框架口径"""
    d = t
    return int((d.p_perm < 0.05).sum()), int((d.q_perm_BH < 0.05).sum()), int(len(d))


def depth_ratio(t):
    d = t[t.score == "SGRG"]
    r = d["median_nUMI_uc"].values / d["median_nUMI_hc"].values
    return float(np.nanmedian(r))


n_raw, n_raw_q, n_raw_tot = n_sig(raw_sg)
n_ds1, n_ds1_q, n_ds1_tot = n_sig(ds1_sg)
dr_raw, dr_ds1 = depth_ratio(raw_t), depth_ratio(ds1_t)

A2 = {g: df[df.annot_status == g] for g in ["批次伪影", "无问题", "已校正", "信息不足"]}


def _med(v, col):
    return float(pd.to_numeric(v[col], errors="coerce").median())


def _n(v, col):
    return int(pd.to_numeric(v[col], errors="coerce").notna().sum())


bulk_auc = {g: _med(v, "auc_pc1_batch") for g, v in A2.items()}
bulk_nom = {g: float(v["d_bio"].median()) for g, v in A2.items()}
bulk_n = {g: _n(v, "auc_pc1_batch") for g, v in A2.items()}

# ==================================================================== 图 1：剂量-反应
fig, axes = plt.subplots(1, 2, figsize=(174 * MM, 62 * MM))

# ---- a) A2（bulk）：平台质量分 × (批次占据 − 生物占据)
ax = axes[0]
x = pd.to_numeric(df["geeq_pub_batch_confound"], errors="coerce").values
y = df["d_bio"].values
cls = df["annot_status"].values
cmap = {"批次伪影": "#C44E52", "信息不足": "#DD8452", "无问题": "#4C72B0", "已校正": "#55A868"}
for g in ["信息不足", "无问题", "已校正", "批次伪影"]:
    m = cls == g
    ax.scatter(x[m], y[m], s=11, c=cmap[g], edgecolors="white", linewidths=0.3,
               label={"批次伪影": "Batch artifact", "信息不足": "Insufficient info",
                      "无问题": "No problem", "已校正": "Batch corrected"}[g], zorder=3)
m = np.isfinite(x) & np.isfinite(y)
z = np.polyfit(x[m], y[m], 1)
xx = np.linspace(np.nanmin(x[m]), np.nanmax(x[m]), 50)
ax.plot(xx, np.polyval(z, xx), color="0.25", lw=0.9, ls="--", zorder=2)
rho, p, n = sp(x.astype(float), y.astype(float))
ax.axhline(0, color="0.75", lw=0.6, zorder=1)
ax.set_xlabel("Gemma GEEQ batch score\n(higher = cleaner)", fontsize=7)
ax.set_ylabel("AUC$_{batch}$ − AUC$_{biology}$  (PC1)", fontsize=7)
ax.set_title("a  Bulk cohorts (n = %d)\n$\\rho$ = %+.2f, %s" % (n, rho, pfmt(p)),
             fontsize=7, loc="left")
ax.legend(fontsize=5.6, loc="lower right", handletextpad=0.3, borderpad=0.2)

# ---- b) Part D（单细胞）：深度-疾病共线强度 × 效应量收缩幅度
ax = axes[1]
x2 = coll["r2_depth_disease"].values
y2 = coll["shrink_abs"].values
sc = ax.scatter(x2, y2, s=13, c=x2, cmap="viridis", edgecolors="white", linewidths=0.3, zorder=3)
m2 = np.isfinite(x2) & np.isfinite(y2)
z2 = np.polyfit(x2[m2], y2[m2], 1)
xx2 = np.linspace(np.nanmin(x2[m2]), np.nanmax(x2[m2]), 50)
ax.plot(xx2, np.polyval(z2, xx2), color="0.25", lw=0.9, ls="--", zorder=2)
thr = float(np.median(x2[m2]))
ax.axvline(thr, color="0.8", lw=0.6, ls=":", zorder=1)
rho2, p2, n2 = sp(x2, y2)
ax.set_xlabel("R²(log nUMI ~ disease)  (donor level)", fontsize=7)
ax.set_ylabel("Shrinkage of effect size after downsampling", fontsize=7)
ax.set_title("b  Single-cell, ds4 (n = %d cell types)\n$\\rho$ = %+.2f, %s" % (n2, rho2, pfmt(p2)),
             fontsize=7, loc="left")
ax.text(0.98, 0.04, "high-collinearity cell types:\nuncorrected $p$<0.05  %d/%d → %d/%d"
        % (n_raw, n_raw_tot, n_ds1, n_ds1_tot),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8, color="0.3")

fig.suptitle("Same qualitative law at two resolutions — metrics differ, panels are NOT poolable",
             fontsize=7.2, x=0.01, ha="left", y=1.045)
fig.text(0.01, -0.10,
         "a  Bulk: platform batch-corrected cohorts have a cleaner technical axis, and the "
         "biological contrast reclaims PC1 (negative slope).\n"
         "b  Single-cell: the more the donor-level sequencing depth is collinear with disease, "
         "the more the nominal effect shrinks after depth matching.",
         fontsize=6.0, color="0.25", ha="left", va="top")
save(fig, "fig_AD1_dose_response")

# ==================================================================== 图 2：干预-响应
fig, axes = plt.subplots(2, 2, figsize=(174 * MM, 96 * MM))
XL, XR = 0, 1
ARROW = dict(arrowstyle="-|>", color="0.45", lw=0.8, shrinkA=2, shrinkB=2)


def slopegraph(ax, v0, v1, l0, l1, col, ylab, title, fmt="%.3f", better=""):
    ax.plot([XL, XR], [v0, v1], color=col, lw=1.3, zorder=3,
            marker="o", ms=4.5, mfc=col, mec="white", mew=0.6)
    ax.annotate("", xy=(XR, v1), xytext=(XL, v0), arrowprops=ARROW, zorder=2)
    ax.scatter([XL, XR], [v0, v1], s=20, c=col, zorder=4, edgecolors="white", linewidths=0.5)
    ax.text(XL - 0.06, v0, fmt % v0, ha="right", va="center", fontsize=6.2, color=col)
    ax.text(XR + 0.06, v1, fmt % v1, ha="left", va="center", fontsize=6.2, color=col)
    ax.set_xlim(-0.42, 1.42)
    pad = 0.22 * max(abs(v1 - v0), 1e-9)
    lo, hi = min(v0, v1) - pad, max(v0, v1) + pad
    if lo < 0 < hi or (min(v0, v1) < 0 < max(v0, v1)):
        lo, hi = min(v0, v1, 0) - pad, max(v0, v1, 0) + pad
        ax.axhline(0, color="0.85", lw=0.6, zorder=1)
    ax.set_ylim(lo, hi)
    ax.set_xticks([XL, XR])
    ax.set_xticklabels([l0, l1], fontsize=6.2)
    ax.set_ylabel(ylab, fontsize=6.6)
    ax.set_title(title, fontsize=6.8, loc="left")
    if better:
        ax.text(0.5, 0.97, better, transform=ax.transAxes, ha="center", va="top",
                fontsize=5.8, color="0.35")


# ---- a) bulk 技术轴
slopegraph(axes[0, 0], bulk_auc["批次伪影"], bulk_auc["已校正"],
           "Before\n(batch artifact, n=%d)" % bulk_n["批次伪影"],
           "After\n(batch corrected, n=%d)" % bulk_n["已校正"],
           C_TECH, "PC1 batch AUC", "a  Bulk · technical axis",
           better="↓ %.3f" % (bulk_auc["批次伪影"] - bulk_auc["已校正"]))
# ---- b) bulk 名义信号
slopegraph(axes[0, 1], bulk_nom["批次伪影"], bulk_nom["已校正"],
           "Before\n(n=%d)" % bulk_n["批次伪影"], "After\n(n=%d)" % bulk_n["已校正"],
           C_NOM, "AUC$_{batch}$ − AUC$_{biology}$", "b  Bulk · nominal signal",
           fmt="%+.3f", better="sign flips: batch stops dominating")
# ---- c) sc 技术轴
slopegraph(axes[1, 0], dr_raw, dr_ds1,
           "Raw\n(no control)", "Downsampled\nto 1,000 UMI",
           C_TECH, "Depth ratio UC / control (median)", "c  Single-cell · technical axis",
           fmt="%.2f", better="↓ %.2f" % (dr_raw - dr_ds1))
# ---- d) sc 名义信号
slopegraph(axes[1, 1], n_raw, n_ds1,
           "Raw\n(no control)", "Downsampled\nto 1,000 UMI",
           C_NOM, "Cell types, uncorrected $p$ < 0.05 (SGRG)", "d  Single-cell · nominal signal",
           fmt="%.0f", better="%d / %d  →  %d / %d\n(FDR $q$ < 0.05: %d → %d)"
           % (n_raw, n_raw_tot, n_ds1, n_ds1_tot, n_raw_q, n_ds1_q))

fig.suptitle("Two independent interventions remove the technical structure AND the nominal signal",
             fontsize=7.4, x=0.01, ha="left", y=1.035)
fig.text(0.01, -0.055,
         "Bulk (a,b): platform batch correction.  Single-cell (c,d): binomial downsampling to a common depth.\n"
         "Neither correction is a regression adjustment — both are interventions on the data, "
         "which is the point: within-unit collinearity makes residual-style adjustment unidentifiable.",
         fontsize=6.0, color="0.25", ha="left", va="top")
fig.tight_layout()
save(fig, "fig_AD2_intervention_response")

# ==================================================================== 数字落盘
summ = {
    "bulk_auc_before_artifact": bulk_auc["批次伪影"], "bulk_auc_after_corrected": bulk_auc["已校正"],
    "bulk_auc_noproblem": bulk_auc["无问题"],
    "bulk_nom_before": bulk_nom["批次伪影"], "bulk_nom_after": bulk_nom["已校正"],
    "bulk_nom_noproblem": bulk_nom["无问题"],
    "n_artifact": bulk_n["批次伪影"], "n_corrected": bulk_n["已校正"], "n_noproblem": bulk_n["无问题"],
    "sc_depth_ratio_raw": dr_raw, "sc_depth_ratio_ds1k": dr_ds1,
    "sc_nsig_raw": n_raw, "sc_nsig_q_raw": n_raw_q, "sc_nsig_raw_total": n_raw_tot,
    "sc_nsig_ds1k": n_ds1, "sc_nsig_q_ds1k": n_ds1_q, "sc_nsig_ds1k_total": n_ds1_tot,
    "dose_response_bulk": {"rho": rho, "p": p, "n": n},
    "dose_response_sc": {"rho": rho2, "p": p2, "n": n2},
}
json.dump(summ, open(os.path.join(OUT, "partAD_figure_numbers.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1, default=float)
print(json.dumps(summ, ensure_ascii=False, indent=1, default=float))

# ==================================================================== 图 3：随机基因集零分布
fig, axes = plt.subplots(1, 2, figsize=(174 * MM, 70 * MM), sharex=True, sharey=True)
ARMS = [("ds4_tearle_rand", "a  Raw (no depth control)", 26),
        ("ds4_tearle_rand_ds1k", "b  Downsampled to 1,000 UMI", 25)]
POS_COL = {"H_INFLAMMATORY_RESPONSE": "#DD8452", "H_TNFA_SIGNALING_VIA_NFKB": "#8172B2",
           "H_INTERFERON_GAMMA_RESPONSE": "#55A868"}
POS_LAB = {"H_INFLAMMATORY_RESPONSE": "Inflammatory response",
           "H_TNFA_SIGNALING_VIA_NFKB": "TNFa/NF-kB",
           "H_INTERFERON_GAMMA_RESPONSE": "IFN-gamma response"}

for ax, (tag, ttl, nct) in zip(axes, ARMS):
    s = pd.read_csv(OUT + "/partD_%s_geneset_summary.tsv" % tag, sep="\t")
    rnd = s[s.kind == "random"]
    real = s[s.kind == "real"]
    pos = s[s.kind == "positive"]
    exp_fp = nct * 0.05
    rng = np.random.default_rng(3)
    jx = rng.normal(0, 0.16, len(rnd))
    jy = rng.normal(0, 0.16, len(rnd))
    ax.scatter(rnd.n_p05_perm + jx, rnd.n_p05 + jy, s=7, c="0.78", edgecolors="none",
               label="Expression-matched\nrandom sets (n=%d)" % len(rnd), zorder=2)
    ax.axvline(exp_fp, color="0.6", lw=0.7, ls=":", zorder=1)
    ax.axhline(exp_fp, color="0.6", lw=0.7, ls=":", zorder=1)
    ax.text(0.02, 0.98, "dotted lines: expected under true null (%.1f cell types)" % exp_fp,
            transform=ax.transAxes, fontsize=5.4, color="0.45", ha="left", va="top")
    for _, r in real.iterrows():
        if r.geneset == "SGRG":
            ax.scatter(r.n_p05_perm, r.n_p05, marker="*", s=95, c="#C44E52",
                       edgecolors="white", linewidths=0.5, zorder=5, label="SGRG panel")
            ax.annotate("SGRG", (r.n_p05_perm, r.n_p05), textcoords="offset points",
                        xytext=(7, 5), fontsize=6.2, color="#C44E52", fontweight="bold")
        else:
            ax.scatter(r.n_p05_perm, r.n_p05, marker="o", s=16, c="#C44E52", alpha=0.55,
                       edgecolors="white", linewidths=0.4, zorder=4)
    for _, r in pos.iterrows():
        ax.scatter(r.n_p05_perm, r.n_p05, marker="D", s=22, c=POS_COL[r.geneset],
                   edgecolors="white", linewidths=0.4, zorder=5,
                   label="Hallmark (positive ctrl)" if r.geneset == "H_INFLAMMATORY_RESPONSE" else None)
    med = rnd.n_p05_perm.median()
    ax.axvline(med, color="0.5", lw=0.7, ls="--", zorder=1)
    ax.set_xlim(-1.5, nct * 0.82)
    ax.set_ylim(-1.5, nct * 0.62)
    ax.set_xlabel("Cell types, permutation $p$ < 0.05", fontsize=6.8)
    ax.set_title("%s\nrandom median: %d (permutation) / %d (Mann-Whitney)" % (
        ttl, med, rnd.n_p05.median()), fontsize=6.8, loc="left")
axes[0].set_ylabel("Cell types, Mann-Whitney $p$ < 0.05", fontsize=6.8)
axes[0].legend(fontsize=5.4, loc="lower right", handletextpad=0.3, borderpad=0.3)
fig.suptitle("Nominal signal carries no panel-specific information: SGRG sits inside the random null",
             fontsize=7.4, x=0.01, ha="left", y=1.04)
fig.text(0.01, -0.09,
         "Each grey point = one of 200 random gene sets matched to SGRG in size and mean-expression profile.\n"
         "Both test statistics are strongly inflated without depth control (expected false positives 1.3 / 1.2 cell types); "
         "downsampling restores calibration, but SGRG never leaves the random cloud — and neither do the Hallmark positive controls.",
         fontsize=6.0, color="0.25", ha="left", va="top")
save(fig, "fig_AD3_random_null")
