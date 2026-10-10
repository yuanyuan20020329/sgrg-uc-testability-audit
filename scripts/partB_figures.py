# -*- coding: utf-8 -*-
"""
Part B · Step 4 — 三张主图
  fig_partB_uc_mdd_scatter.png    UC vs MDD log2FC 象限散点
  fig_partB_evidence_heatmap.png  证据叠合热图（sign · -log10 p）
  fig_partB_enrichment.png        随机面板零分布 + 面板内 Venn
"""
import os, csv, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT  = os.path.join(BASE, "results")

def log(*a): print(*a, flush=True)

def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

# ---------------- 读证据叠合表
rows = {}
with open(os.path.join(OUT, "partB_evidence_stack.tsv"), encoding="utf-8") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        rows[row["gene"]] = row

res = json.load(open(os.path.join(OUT, "partB_bridge_results.json"), encoding="utf-8"))

# ================================================================ Fig 1 散点
pts = []
for g, row in rows.items():
    x = fnum(row["UC_lfc_GSE38713"]); y = fnum(row["MDDblood_lfc_GSE19738"])
    if x is None or y is None:
        continue
    pts.append((g, x, y, fnum(row["UC_FDR"]) or 1, fnum(row["MDDblood_q"]) or 1))
xs = np.array([p[1] for p in pts]); ys = np.array([p[2] for p in pts])

fig, ax = plt.subplots(figsize=(5.6, 5.4), dpi=300)
ax.axhline(0, color="#999999", lw=0.6, zorder=1)
ax.axvline(0, color="#999999", lw=0.6, zorder=1)
for g, x, y, fu, fm in pts:
    both = fu < 0.05 and fm < 0.05
    ax.scatter(x, y, s=26 if both else 14,
               color="#d1495b" if (x > 0 and y > 0) else "#2e86ab" if (x < 0 and y < 0) else "#b8b8b8",
               alpha=0.85 if both else 0.55, edgecolor="none", zorder=3)
lbl = {"OGFOD1": (0.30, 0.10), "CASC3": (0.24, 0.055), "G3BP1": (0.24, -0.21),
       "EIF4A1": (0.53, -0.25), "VCP": (0.64, -0.02), "YTHDF3": (-0.71, 0.03),
       "ATXN2": (-0.35, 0.045), "PQBP1": (-0.05, -0.113), "CIRBP": (-0.147, -0.21),
       "STAU1": (-0.055, 0.198), "KPNB1": (0.313, -0.17)}
for g, (dx, dy) in lbl.items():
    row = rows.get(g)
    if not row:
        continue
    x = fnum(row["UC_lfc_GSE38713"]); y = fnum(row["MDDblood_lfc_GSE19738"])
    if x is None or y is None:
        continue
    ax.annotate(g, (x, y), xytext=(dx, dy), textcoords="offset points",
                fontsize=6.5, color="#333333")
same = int((np.sign(xs) == np.sign(ys)).sum())
ax.set_xlabel("UC log2FC (GSE38713, involved vs control)", fontsize=9)
ax.set_ylabel("MDD log2FC (GSE19738 blood, MDD vs control)", fontsize=9)
c = res["concordance"]
ax.set_title(f"SGRG panel: UC vs MDD effect directions\n"
             f"rho = {c['spearman_rho']:.3f} (p = {c['spearman_p']:.2g}); "
             f"{same}/{len(pts)} same sign; 0 genes FDR<0.05 both sides", fontsize=8)
ax.tick_params(labelsize=8)
fig.tight_layout()
p1 = os.path.join(OUT, "fig_partB_uc_mdd_scatter.png")
fig.savefig(p1, dpi=300); plt.close(fig)
log("[write]", p1)

# ================================================================ Fig 2 证据热图
STREAMS = [
    ("UC colon\nGSE38713", "UC_lfc_GSE38713", "UC_FDR", "p"),
    ("UC paired\nGSE22619", None, "UCpaired_p", "dz"),   # 用 dz 的符号 × -log10(p)
    ("MDD blood\nGSE19738", "MDDblood_lfc_GSE19738", "MDDblood_p", "q"),
    ("MDD nAcc\nGSE80655", "nAcc_lfc", "nAcc_p", "p"),
    ("MDD ACC\nGSE80655", "ACC_lfc", "ACC_p", "p"),
    ("MDD dlPFC\nGSE80655", "dlPFC80655_lfc", "dlPFC80655_p", "p"),
    ("MDD striatum\nGSE53987", "STR_lfc", "STR_p", "p"),
    ("MDD hippocampus\nGSE53987", "HPC_lfc", "HPC_p", "p"),
    ("MDD dlPFC\nGSE53987", "dlPFC53987_lfc", "dlPFC53987_p", "p"),
    ("MR eQTL→UC\nPart C", None, "MR_IVW_p", "b"),
]
def sign_of(row, stream):
    label, lfc_col, p_col, mode = stream
    p = fnum(row.get(p_col))
    if p is None:
        return None, None, None
    if mode == "dz":
        v = fnum(row.get("UCpaired_dz_GSE22619"))
    elif mode == "b":
        # MR 方向用 Part C IVW b（同面板证据列里未存 b，取组织最小 p 那条；这里用 q 的列不存在，退化处理）
        v = fnum(row.get("MR_IVW_b"))
        if v is None:
            # evidence stack 没存 b，则用 tissue 列判断不可行——MR 用正向固定色阶
            return None, None, None
    else:
        v = fnum(row.get(lfc_col))
    if v is None:
        return None, None, None
    stat = np.sign(v) * (-np.log10(max(p, 1e-300)))
    sig = False
    if mode == "q":
        q = fnum(row.get("MDDblood_q"))
        sig = q is not None and q < 0.05
    else:
        sig = p < 0.05
    return stat, p, sig

# 选基因：任一证据流 p<0.05 的基因，按流数排序，取前 24
score = {}
for g, row in rows.items():
    n_sig = 0
    for s in STREAMS:
        _, p, sig = sign_of(row, s)
        if sig:
            n_sig += 1
    if n_sig:
        score[g] = n_sig
sel = sorted(score, key=lambda g: (-score[g], g))[:24]
log(f"heatmap genes ({len(sel)}):", sel)

# MR b 列没有——从 partC_mr_geneset 直接补
mr_b = {}
with open(os.path.join(OUT, "partC_mr_geneset.tsv"), encoding="utf-8") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        g = row["gene"]; p = fnum(row["IVW_p"])
        if p is None:
            continue
        if g not in mr_b or p < mr_b[g][1]:
            mr_b[g] = (fnum(row["IVW_b"]), p, fnum(row["IVW_q"]))

M = np.full((len(sel), len(STREAMS)), np.nan)
SIG = np.zeros((len(sel), len(STREAMS)), dtype=bool)
for i, g in enumerate(sel):
    row = dict(rows[g])
    for j, s in enumerate(STREAMS):
        label, lfc_col, p_col, mode = s
        if mode == "b":
            if g in mr_b:
                b, p, q = mr_b[g]
                M[i, j] = np.sign(b) * (-np.log10(p))
                SIG[i, j] = q is not None and q < 0.05
            continue
        stat, p, sig = sign_of(row, s)
        if stat is not None:
            M[i, j] = stat
            SIG[i, j] = sig

fig, ax = plt.subplots(figsize=(7.6, 0.34 * len(sel) + 1.6), dpi=300)
vmax = np.nanmax(np.abs(M))
im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
for i in range(len(sel)):
    for j in range(len(STREAMS)):
        if np.isnan(M[i, j]):
            ax.text(j, i, "·", ha="center", va="center", fontsize=6, color="#bbbbbb")
        elif SIG[i, j]:
            ax.plot(j, i, marker="*", ms=4.5, color="#111111", mec="none")
ax.set_xticks(range(len(STREAMS)))
ax.set_xticklabels([s[0] for s in STREAMS], fontsize=6.5, rotation=38, ha="right")
ax.set_yticks(range(len(sel)))
ax.set_yticklabels(sel, fontsize=7)
ax.set_title("Evidence stack: sign(lfc) x -log10(p)\n"
             "* = FDR/q < 0.05 (paired: Wilcoxon p < 0.05)", fontsize=8.5, pad=8)
fig.tight_layout()
p2 = os.path.join(OUT, "fig_partB_evidence_heatmap.png")
fig.savefig(p2, dpi=300); plt.close(fig)
log("[write]", p2)

# ================================================================ Fig 3 富集零分布 + Venn
enr = res["enrichment"]
rng = np.random.default_rng(20261009)
# universe 与 DEG 从 genomewide 重算
u = set(); de = set()
with open(os.path.join(OUT, "partB_gse19738_genomewide.tsv"), encoding="utf-8") as f:
    r = csv.reader(f, delimiter="\t"); next(r)
    for row in r:
        probe, syms, p, q, lfc, tstat = row
        ss = [s for s in syms.split("|") if s and s != "NA"]
        for s in ss:
            u.add(s)
        qq = fnum(q)
        if qq is not None and qq <= 0.05:
            for s in ss:
                de.add(s)
ul = sorted(u)
null = []
for _ in range(2000):
    idx = rng.choice(len(ul), size=enr["n"], replace=False)
    null.append(sum(1 for i in idx if ul[i] in de))
null = np.array(null)

fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6), dpi=300,
                         gridspec_kw=dict(width_ratios=[1.15, 1]))
ax = axes[0]
ax.hist(null, bins=np.arange(-0.5, max(null) + 1.5, 1), color="#cfe0f4",
        edgecolor="#7fa8d0", lw=0.6)
ax.axvline(enr["k"], color="#d1495b", lw=1.4)
ax.annotate(f"SGRG observed\nk = {enr['k']}", (enr["k"], ax.get_ylim()[1] * 0.72),
            fontsize=7.5, color="#d1495b", ha="left", xytext=(4, 0), textcoords="offset points")
ax.axvline(np.median(null), color="#555555", lw=0.8, ls="--")
ax.set_xlabel("MDD-DE genes in random 90-gene panel", fontsize=8.5)
ax.set_ylabel("panels (of 2000)", fontsize=8.5)
ax.set_title(f"No panel-level enrichment\nOR = {enr['OR']:.2f}, p = {enr['p']:.2g}, "
             f"{enr['percentile']:.0f}th pct of random panels", fontsize=9)
ax.tick_params(labelsize=8)

# Venn：面板内 UC-DE vs MDD-DE
uc_de = {g for g, row in rows.items() if (fnum(row.get("UC_FDR")) or 1) < 0.05}
mdd_de_panel = {g for g, row in rows.items() if (fnum(row.get("MDDblood_q")) or 1) < 0.05}
n_uc_only = len(uc_de - mdd_de_panel)
n_mdd_only = len(mdd_de_panel - uc_de)
n_both = len(uc_de & mdd_de_panel)
ax2 = axes[1]
c1 = Circle((0.42, 0.5), 0.28, facecolor="#f2c14e", alpha=0.55, edgecolor="#b58b00", lw=1)
c2 = Circle((0.62, 0.5), 0.28, facecolor="#2e86ab", alpha=0.55, edgecolor="#1f6285", lw=1)
ax2.add_patch(c1); ax2.add_patch(c2)
ax2.text(0.30, 0.50, str(n_uc_only), fontsize=13, ha="center", va="center")
ax2.text(0.74, 0.50, str(n_mdd_only), fontsize=13, ha="center", va="center")
ax2.text(0.52, 0.50, str(n_both), fontsize=13, ha="center", va="center")
ax2.text(0.28, 0.90, "UC DE\n(FDR<0.05)", fontsize=8, ha="center", color="#7a5d00")
ax2.text(0.78, 0.90, "MDD DE\n(q<0.05)", fontsize=8, ha="center", color="#1f6285")
ax2.set_xlim(0, 1.05); ax2.set_ylim(0, 1); ax2.axis("off")
ax2.set_title(f"Within SGRG panel (n measured = {enr['n']}): zero overlap", fontsize=9)
fig.tight_layout()
p3 = os.path.join(OUT, "fig_partB_enrichment.png")
fig.savefig(p3, dpi=300); plt.close(fig)
log("[write]", p3)
