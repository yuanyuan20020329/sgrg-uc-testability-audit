# -*- coding: utf-8 -*-
"""
Part B addendum — OGFOD1 locus sensitivity figure.

fig_partB_ogfod1_locus.png : three-track regional view of chr16:55.6-56.6 Mb
  track 1  UC GWAS (de Lange 2017) -log10 p
  track 2  OGFOD1 Whole_Blood cis-eQTL -log10 p, coloured by r^2 to its lead
  track 3  CES1P1 Whole_Blood cis-eQTL -log10 p, coloured by r^2 to its lead

Key annotations: coloc.abf posteriors, lead-lead LD, the effective number of
independent MR instruments, and the lead-SNP Wald vs naive IVW contrast.
"""
import os
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
RES = os.path.join(BASE, "results")
PANEL = r"D:/ld_ref_1000g/regions/chr16_ogfod1_grch38.tsv"
OUT = os.path.join(RES, "fig_partB_ogfod1_locus.png")

X0, X1 = 55.60e6, 56.60e6
LEAD = {"OGFOD1": ("rs4784670", 56464599), "CES1P1": ("rs12149359", 55832988)}
# Okabe-Ito-derived LD ramp
LD_COLORS = [(0.9, 0.9, 0.9), (0.75, 0.85, 0.95), (0.45, 0.68, 0.90),
             (0.15, 0.45, 0.80), (0.05, 0.25, 0.65), (0.65, 0.20, 0.15)]


def ld_color(r2):
    if r2 is None:
        return (0.78, 0.78, 0.78)
    for cut, c in zip((0.2, 0.4, 0.6, 0.8, 0.95), LD_COLORS[1:]):
        if r2 < cut:
            return c
    return LD_COLORS[-1]


# ---------------- panel dosages ----------------
panel = {}
with open(PANEL, encoding="utf-8") as f:
    f.readline()
    for line in f:
        fld = line.rstrip("\n").split("\t")
        if len(fld) < 4:
            continue
        panel[(int(fld[0]), tuple(sorted((fld[1], fld[2]))))] = \
            np.array([np.nan if x == "" else float(x) for x in fld[4:]])


def r2_with(lead_key, pos, ref, alt):
    a = panel.get(lead_key)
    b = panel.get((pos, tuple(sorted((ref, alt)))))
    if a is None or b is None:
        return None
    m = ~(np.isnan(a) | np.isnan(b))
    if m.sum() < 100:
        return None
    r = float(np.corrcoef(a[m], b[m])[0, 1])
    return r * r


lead_keys = {g: (LEAD[g][1], None) for g in LEAD}   # filled below
# the lead entries themselves need their allele pair — rebuild from GTEx cache
cache = json.load(open(os.path.join(BASE, "data", "raw", "gtex_eqtl",
                                    "OGFOD1__Whole_Blood.json"), encoding="utf-8"))
gtex_var = {}
for r in cache.get("rows", []):
    vid = (r.get("variantId") or "").split("_")
    if len(vid) >= 5:
        gtex_var[r.get("snpId")] = (int(vid[1]), vid[2].upper(), vid[3].upper(),
                                    float(r.get("pValue")), r.get("geneSymbolUpper"))
for g, (rsid, pos) in LEAD.items():
    if rsid in gtex_var:
        v = gtex_var[rsid]
        lead_keys[g] = (v[0], tuple(sorted((v[1], v[2]))))

# ---------------- track data ----------------
# GWAS (all extracted regional SNPs)
gwas_x, gwas_y, gwas_c = [], [], []
with open(os.path.join(RES, "ogfod1_region_gwas.tsv"), encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        try:
            pos = int(r["hm_pos"])
            p = float(r["p_value"])
        except (TypeError, ValueError):
            continue
        if not (X0 <= pos <= X1) or p <= 0:
            continue
        ea, oa = (r.get("hm_effect_allele") or "").upper(), (r.get("hm_other_allele") or "").upper()
        if not ea or not oa:
            continue
        gwas_x.append(pos / 1e6)
        gwas_y.append(-np.log10(p))
        gwas_c.append(ld_color(r2_with(lead_keys["OGFOD1"], pos, ea, oa)))

# eQTL tracks from the raw cache (all SNPs of each gene)
tracks = {}
for gene in ("OGFOD1", "CES1P1"):
    xs, ys, cs = [], [], []
    for rsid, (pos, ref, alt, p, gsym) in gtex_var.items():
        if gsym != gene or not (X0 <= pos <= X1) or p <= 0:
            continue
        xs.append(pos / 1e6)
        ys.append(-np.log10(p))
        cs.append(ld_color(r2_with(lead_keys[gene], pos, ref, alt)))
    tracks[gene] = (np.array(xs), np.array(ys), cs)

# ---------------- figure ----------------
fig, axes = plt.subplots(3, 1, figsize=(7.4, 6.6), dpi=300,
                         sharex=True, gridspec_kw={"height_ratios": [1, 1, 1],
                                                   "hspace": 0.18})

# track 1 — UC GWAS
ax = axes[0]
ax.scatter(gwas_x, gwas_y, s=9, c=gwas_c, edgecolors="none", zorder=2)
ax.set_ylabel("UC GWAS\n$-\\log_{10}p$", fontsize=8)
ax.set_ylim(-0.1, 2.7)
ax.set_title("chr16 OGFOD1 / CES1P1 locus — colocalisation sensitivity analysis",
             fontsize=9.5)
ax.text(0.985, 0.95, "no SNP reaches genome-wide significance\n"
                     "min p = 0.014 (rs3790110)\n"
                     "coloc.abf vs OGFOD1 eQTL: PP.H4 = 0.044 (PP.H1 = 0.951)\n"
                     "coloc.abf vs CES1P1 eQTL: PP.H4 = 0.096 (PP.H1 = 0.901)",
        transform=ax.transAxes, ha="right", va="top", fontsize=6.8,
        bbox=dict(fc="#fdf6e3", ec="#c9b45c", lw=0.6, pad=3))

# track 2 — OGFOD1 eQTL
ax = axes[1]
xs, ys, cs = tracks["OGFOD1"]
ax.scatter(xs, ys, s=9, c=cs, edgecolors="none", zorder=2)
i = int(np.argmax(ys))
ax.scatter([xs[i]], [ys[i]], s=26, marker="D", c="#d95f02", zorder=3)
ax.annotate("rs4784670 (OGFOD1 lead, p=4.1e-9)\nUC GWAS p = 0.118 (Wald ratio)",
            xy=(xs[i], ys[i]), xytext=(xs[i] + 0.13, ys[i] + 1.3),
            fontsize=6.8, arrowprops=dict(arrowstyle="-", lw=0.6, color="#666666"))
ax.set_ylabel("OGFOD1 eQTL\n$-\\log_{10}p$", fontsize=8)
ax.set_ylim(-0.1, 9.6)
ax.text(0.985, 0.06, "10 MR instruments: pairwise |r| median 0.94 — one LD block\n"
                     "effective independent instruments = 1.09\n"
                     "lead Wald p = 0.118  vs  naive IVW p = 1.7e-7 (pseudo-replication)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=6.8,
        bbox=dict(fc="#e8f4f8", ec="#7fb3c8", lw=0.6, pad=3))

# track 3 — CES1P1 eQTL
ax = axes[2]
xs, ys, cs = tracks["CES1P1"]
ax.scatter(xs, ys, s=9, c=cs, edgecolors="none", zorder=2)
i = int(np.argmax(ys))
ax.scatter([xs[i]], [ys[i]], s=26, marker="D", c="#d95f02", zorder=3)
ax.annotate("rs12149359 (CES1P1 lead, p=5.2e-97)\n"
            "r vs OGFOD1 lead = -0.06 (r$^2$ = 0.004) — different LD block",
            xy=(xs[i], ys[i]), xytext=(xs[i] + 0.14, ys[i] - 28),
            fontsize=6.8, arrowprops=dict(arrowstyle="-", lw=0.6, color="#666666"))
ax.set_ylabel("CES1P1 eQTL\n$-\\log_{10}p$", fontsize=8)
ax.set_xlabel("chr16 position (GRCh38, Mb)", fontsize=8.5)
ax.set_ylim(-0.5, 105)
ax.set_xlim(X0 / 1e6, X1 / 1e6)

# gene marks
for ax_ in axes:
    for gene, (rsid, pos) in LEAD.items():
        ax_.axvline(pos / 1e6, color="#999999", lw=0.5, ls=":", zorder=1)
axes[2].text(LEAD["OGFOD1"][1] / 1e6 + 0.008, 92, "OGFOD1", fontsize=7.5,
             style="italic", color="#444444")
axes[2].text(LEAD["CES1P1"][1] / 1e6 + 0.008, 92, "CES1P1", fontsize=7.5,
             style="italic", color="#444444")
axes[0].text(0.012, 0.10, "points coloured by r$^2$ with the OGFOD1 lead (1000G high-coverage EUR, n=633)",
             transform=axes[0].transAxes, fontsize=6.2, color="#555555")

for ax_ in axes:
    ax_.tick_params(labelsize=7)

fig.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close(fig)
print("saved", OUT)
