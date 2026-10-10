# -*- coding: utf-8 -*-
"""Figure 1 -- the five-rung testability audit ladder (graphical abstract).

Rendered with matplotlib because headless Chrome/Edge are unavailable in this
sandbox (both binaries run but emit nothing at all, including --version).

Submission constraints honoured:
  * width exactly 17.0 cm (full page width)
  * every glyph >= 6 pt (min used = 6.0 pt)
  * font family Arial (DejaVu Sans fallback); no U+207B anywhere
Outputs: results/fig_F1_audit_ladder.{png,pdf}
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")

plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

W_CM = 17.0
H_CM = 12.85
L, R = 0.18, W_CM - 0.18

FS_TITLE, FS_SUB, FS_BH = 9.5, 6.8, 7.5
FS_Q, FS_K, FS_V, FS_FOOT, FS_TAG = 7.0, 6.1, 6.4, 6.3, 6.0

C_DARK, C_PILL = "#12395E", "#1F4E79"
C_BODY, C_MUT, C_LINE = "#22303C", "#5A6773", "#C9D3DC"
C_BAD_BG, C_BAD_ED, C_BAD_TX = "#FDECEA", "#E7A79E", "#8A1B12"
C_OK_BG, C_OK_ED, C_OK_TX = "#EAF6EC", "#A9D3B2", "#145C29"

fig = plt.figure(figsize=(W_CM / 2.54, H_CM / 2.54))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W_CM)
ax.set_ylim(0, H_CM)
ax.axis("off")


def box(x, y, w, h, fc, ec, r=0.09, lw=0.55, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0,rounding_size=%f" % r,
                                linewidth=lw, edgecolor=ec, facecolor=fc, zorder=z))


def txt(x, y, s, fs=FS_V, color=C_BODY, w="normal", ha="left", va="center", z=3):
    ax.text(x, y, s, fontsize=fs, color=color, fontweight=w,
            ha=ha, va=va, zorder=z, linespacing=1.30)


def arrow(x0, x1, y):
    ax.annotate("", xy=(x1, y), xytext=(x0, y),
                arrowprops=dict(arrowstyle="-|>", color="#7A8794",
                                lw=0.9, shrinkA=0, shrinkB=0), zorder=2)


y = H_CM

# ----------------------------------------------------------------- title
y -= 0.36
txt(L, y, "A five-rung testability audit ladder for gene-set \u00d7 disease studies", FS_TITLE, C_DARK, "bold")
y -= 0.30
sub = ("Each rung asks one question about whether the public data could, in principle, have returned a negative result. Rungs are ordered by the\n"
       "cost of the check; every rung is answerable from metadata that is already public.")
txt(L, y, sub, FS_SUB, C_MUT, va="top")
y -= 0.78

# ------------------------------------------------------- template band
txt(L, y, "The template being audited", FS_BH, C_DARK, "bold")
y -= 0.24
bh = 0.80
box(L, y - bh - 0.06, R - L, bh + 0.12, "#F5F7F9", "#DCE3E9", r=0.08, lw=0.5, z=0)
labels = [
    ("Curated gene set", "(GO term / published census)"),
    ("Intersect", "with disease DE genes"),
    ("Supervised learner", "(LASSO / SVM-RFE / RF)"),
    ("Reported panel", "\u201cAUC = 1.0\u201d"),
]
n = len(labels)
aw = 0.62
cw = ((R - L) - (n - 1) * aw) / n
for i, (a, b) in enumerate(labels):
    x0 = L + i * (cw + aw)
    box(x0, y - bh, cw, bh, "#FFFFFF", C_LINE, r=0.08, lw=0.6)
    txt(x0 + cw / 2, y - bh * 0.34, a, FS_Q, C_DARK, "bold", ha="center")
    txt(x0 + cw / 2, y - bh * 0.72, b, FS_K, C_MUT, ha="center")
    if i < n - 1:
        arrow(x0 + cw + 0.09, x0 + cw + aw - 0.09, y - bh / 2)
y -= bh
y -= 0.12
nh = 0.56
box(L, y - nh, R - L, nh, C_BAD_BG, C_BAD_ED, r=0.08, lw=0.5)
txt(L + 0.14, y - nh / 2,
    "The template shows that the panel IS associated with disease in one dataset. It never asks whether that dataset could have returned a\n"
    "negative \u2014 so a panel can be \u201cvalidated\u201d by arms that differ in ways unrelated to disease.",
    FS_K, "#6E1109", va="center")
y -= nh

# ------------------------------------------------------------- ladder
y -= 0.32
txt(L, y, "The ladder: five questions, five measured quantities, five failure signatures", FS_BH, C_DARK, "bold")
y -= 0.26

rows = [
    ("1", "PROVENANCE", "Where did the controls\ncome from?", "read one column of metadata",
     "fraction of control cells / samples primary to the study\n(is_primary_data, donor_id)",
     "controls largely re-used from\nexternal donors"),
    ("2", "NESTING", "Is the technical axis nested\nin the disease axis?", "metadata + one PCA",
     "Cram\u00e9r's V / \u03c7\u00b2 between batch and group;\nbatch-only classifier AUC",
     "batch-only AUC \u2248 or > panel AUC"),
    ("3", "TECHNICAL\nLOAD", "Is the score a function of\ntechnical load?", "score and depth per cell",
     "\u03c1(score, nUMI) within cell type; R\u00b2(depth ~ disease);\ncase/control depth ratio",
     "strong score\u2013depth correlation\nwith disease-aligned depth"),
    ("4", "CALIBRATION", "Is the null distribution valid\nfor this design?", "\u2265200 random gene sets",
     "random-gene-set null for the SAME test;\npermutation validity / exchangeability check",
     "more \u201csignificant\u201d units than\nthe null predicts"),
    ("5", "POWER\nFLOOR", "Could the design have\ndetected anything?", "known-positive control",
     "positive-control (known pathway) result;\nminimum attainable p (2/2\u207f for n pairs)",
     "known-positive controls all\nnon-significant"),
]

rh, gap = 0.98, 0.05
pill_w, q_w, sig_w = 1.60, 4.15, 3.72
rep_w = (R - L) - pill_w - q_w - sig_w

for i, (num, tag, q, sub2, rep, sig) in enumerate(rows):
    ry = y - (i + 1) * rh - i * gap
    box(L, ry, pill_w, rh, C_PILL, C_PILL, r=0.08, lw=0.5)
    two = "\n" in tag
    txt(L + pill_w / 2, ry + rh * (0.76 if two else 0.68), num, 8.4, "#FFFFFF", "bold", ha="center")
    txt(L + pill_w / 2, ry + rh * (0.30 if two else 0.26), tag, FS_TAG, "#D6E6F5", "bold", ha="center")

    box(L + pill_w, ry, R - L - pill_w, rh, "#FFFFFF", C_LINE, r=0.0, lw=0.55)

    qx = L + pill_w + 0.12
    txt(qx, ry + rh * (0.64 if "\n" in q else 0.62), q, FS_Q, C_DARK, "bold", va="center")
    txt(qx, ry + rh * 0.20, sub2, FS_K, C_MUT, va="center")

    ax.plot([L + pill_w + q_w] * 2, [ry + 0.03, ry + rh - 0.03], color="#E2E8ED", lw=0.5, zorder=2)
    rx = L + pill_w + q_w + 0.12
    txt(rx, ry + rh - 0.17, "Quantity to report", FS_K, "#7A8794", va="center")
    txt(rx, ry + rh - 0.34, rep, FS_V, C_BODY, va="top")

    sx = L + pill_w + q_w + rep_w
    box(sx, ry, sig_w, rh, C_BAD_BG, "#F2D6D2", r=0.0, lw=0.55)
    txt(sx + 0.12, ry + rh - 0.17, "Failure signature", FS_K, "#A6564D", va="center")
    txt(sx + 0.12, ry + rh - 0.34, sig, FS_V, C_BAD_TX, va="top")

y -= 5 * rh + 4 * gap

# ----------------------------------------------------------- outcomes
y -= 0.34
txt(L, y, "What the ladder returns", FS_BH, C_DARK, "bold")
y -= 0.24
oh = 1.32
ohw = (R - L - 0.22) / 2
box(L, y - oh, ohw, oh, C_OK_BG, C_OK_ED, r=0.08, lw=0.6)
txt(L + 0.13, y - 0.20, "Passes all five rungs", FS_Q, C_OK_TX, "bold")
txt(L + 0.13, y - 0.40,
    "The dataset can test the claim: report the calibrated result. In\n"
    "this project the one qualifying resource returned a weak,\n"
    "direction-consistent positive that emerged only after calibration.",
    FS_K, C_OK_TX, va="top")
x2 = L + ohw + 0.22
box(x2, y - oh, ohw, oh, C_BAD_BG, C_BAD_ED, r=0.08, lw=0.6)
txt(x2 + 0.13, y - 0.20, "Fails any rung", FS_Q, C_BAD_TX, "bold")
txt(x2 + 0.13, y - 0.40,
    "The data CANNOT answer the question \u2014 a statement about the\n"
    "data, not the biology. A non-significant panel is then\n"
    "uninformative rather than negative.",
    FS_K, C_BAD_TX, va="top")
y -= oh

# ------------------------------------------------------------- footer
y -= 0.17
ax.plot([L, R], [y, y], color="#D3DAE1", lw=0.55, zorder=1)
y -= 0.16
txt(L, y, "Generalises beyond the original gene set.", FS_FOOT, C_DARK, "bold", va="top")
y -= 0.26
txt(L, y,
    "The ladder was applied, unchanged, to two further published claims: ferroptosis \u00d7 hepatocellular carcinoma "
    "(reproducible, in-sample AUC 0.9991 vs reported 1.000,\nbut indistinguishable from 500 same-size random gene sets, "
    "empirical p = 0.611) and pyroptosis \u00d7 Parkinson's disease (reported AUC 0.823 not\nrecoverable from public curated "
    "data under the stated fixed weights, 0.6432; the score is dominated by cohort membership, \u03c1 = \u22120.691\npooled "
    "vs \u22120.164 within).",
    FS_FOOT, C_MUT, va="top")

for ext in ("png", "pdf"):
    p = os.path.join(RES, "fig_F1_audit_ladder." + ext)
    fig.savefig(p, dpi=600 if ext == "png" else None, facecolor="white")
    print("wrote", p)
print("size cm: %.3f x %.3f" % (W_CM, H_CM))
