#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Part C / Stage E — LD-contamination ("neighbour eQTL dominance") diagnostic.

Why this and not colocalisation:
  the GTEx REST endpoint we can reach is
  `/association/singleTissueEqtlByLocation`, i.e. "Get **Significant** Single
  Tissue eQTLs".  It returns only variant-gene pairs that already pass GTEx's
  significance threshold, so the regional SNP set is truncated and a formal
  coloc.abf posterior is not identifiable.  (The full all-pairs file is ~4 GB
  per tissue, has no tabix index in the public bucket, and the v8 all-pairs
  data live in a requester-pays bucket.)

What we CAN do rigorously with the same data:
  a cis-eQTL instrument is only interpretable if the lead SNP of the target gene
  is not simply a marker for a *much stronger* eQTL of a neighbouring gene in
  the same LD block.  If a neighbour within a few hundred kb has an eQTL many
  orders of magnitude stronger, the MR estimate is far more plausibly explained
  by that neighbour (i.e. horizontal pleiotropy / LD contamination) than by the
  target gene.

For every instrument set we therefore report, restricted to the same +/-1 Mb
window:
  * the target gene's lead cis-eQTL SNP and its p-value,
  * the strongest neighbouring gene's lead SNP, its p-value, and the distance
    between the two lead SNPs,
  * delta_log10_p = log10(p_target) - log10(p_neighbour)  (>0 => neighbour is
    stronger), and
  * an explicit contamination flag for |distance| <= 250 kb AND delta >= 5.

Output: results/partC_neighbour_eqtl.tsv
"""
import os
import csv
import json
from collections import defaultdict

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
RES = os.path.join(ROOT, "results")
CACHE = os.path.join(ROOT, "data", "raw", "gtex_eqtl")
OUT = os.path.join(RES, "partC_neighbour_eqtl.tsv")
SETS = json.load(open(os.path.join(RES, "partC_instrument_sets.json"), encoding="utf-8"))

SIG = 5e-8
DIST_KB = 250.0        # proximity cut-off for the flag
DELTA_LOG10 = 5.0      # orders of magnitude by which the neighbour must dominate

# curated, literature-established IBD/UC susceptibility genes (for annotation only)
IBD_GENES = {
    "IL27", "SULT1A1", "SULT1A2", "NOD2", "CARD9", "ATG16L1", "IL23R", "JAK2",
    "STAT3", "IL12B", "IRF5", "TNFSF15", "FCGR2A", "GPR35", "PTGER4", "SMAD3",
    "TNFRSF14", "ICAM1", "CREM", "NR5A2", "SLC22A5", "IL10", "HNF4A", "CCR6",
    "TYK2", "SH2B3", "PIM3", "TNFRSF9", "REL", "TNIP1", "RNF186", "UBE2L3",
    "PRDM1", "IKZF1", "NKX2-3", "NLRC5", "C1orf106", "ZFP90", "MST1", "TCF4",
    "CDH1", "HERC2", "LAMB1", "ZMIZ1", "DAP", "PUS10", "IL2RA", "TNFRSF6B",
    "LSP1", "TNFSF8", "ITLN1", "SLC11A1", "MUC1", "MUC19", "FUT2", "HLA-DRB1",
}


def fnum(x):
    try:
        return float(x)
    except Exception:
        return None


def main():
    rows = []
    for key, info in sorted(SETS.items()):
        gene, tissue = key.split("|")
        cf = os.path.join(CACHE, "%s__%s.json" % (gene, tissue))
        if not os.path.exists(cf):
            continue
        d = json.load(open(cf, encoding="utf-8"))
        win = d.get("window", [0, 0])
        sig = [r for r in d.get("rows", []) if fnum(r.get("pValue")) is not None
               and fnum(r["pValue"]) < SIG]
        if not sig:
            continue
        tgt = [r for r in sig if r.get("geneSymbolUpper") == gene.upper()]
        if not tgt:
            continue
        lead = min(tgt, key=lambda r: r["pValue"])
        byg = defaultdict(list)
        for r in sig:
            if r.get("geneSymbolUpper") != gene.upper():
                byg[r["geneSymbolUpper"]].append(r)
        best = None
        for g, rs in byg.items():
            l = min(rs, key=lambda r: r["pValue"])
            dist = abs(l["pos"] - lead["pos"]) / 1000.0
            delta = (__import__("math").log10(lead["pValue"]) -
                     __import__("math").log10(l["pValue"]))
            cand = {"neighbour": g, "nb_snp": l["snpId"], "nb_pos": l["pos"],
                    "nb_p": l["pValue"], "dist_kb": dist, "delta_log10p": delta,
                    "nb_is_known_ibd": g in IBD_GENES}
            if best is None or cand["delta_log10p"] > best["delta_log10p"]:
                best = cand
        if best is None:
            continue
        best["flag"] = bool(best["dist_kb"] <= DIST_KB and best["delta_log10p"] >= DELTA_LOG10)
        rows.append({
            "gene": gene, "tissue": tissue,
            "lead_snp": lead["snpId"], "lead_pos": lead["pos"], "lead_p": lead["pValue"],
            "n_sig_snps": len(tgt),
            "n_neighbour_genes": len(byg),
            **{k: best[k] for k in ("neighbour", "nb_snp", "nb_pos", "nb_p",
                                    "dist_kb", "delta_log10p", "nb_is_known_ibd", "flag")},
        })

    cols = ["gene", "tissue", "lead_snp", "lead_pos", "lead_p", "n_sig_snps",
            "n_neighbour_genes", "neighbour", "nb_snp", "nb_pos", "nb_p",
            "dist_kb", "delta_log10p", "nb_is_known_ibd", "flag"]
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            out = dict(r)
            # keep p-values at full precision (they span ~10^-153 .. 1);
            # only the geometric columns are rounded for readability.
            for c in ("lead_p", "nb_p"):
                out[c] = "%.6e" % float(r[c])
            for c in ("dist_kb", "delta_log10p"):
                out[c] = round(float(r[c]), 2)
            w.writerow({c: out[c] for c in cols})

    print("%-10s %-18s %-12s %-12s %8s %9s %6s" %
          ("gene", "tissue", "target_p", "neighbour", "dist_kb", "delta_log10", "flag"))
    for r in sorted(rows, key=lambda x: -x["delta_log10p"]):
        print("%-10s %-18s %-12.2e %-12s %8.0f %9.1f %6s" %
              (r["gene"], r["tissue"], r["lead_p"], r["neighbour"],
               r["dist_kb"], r["delta_log10p"], "YES" if r["flag"] else ""))
    nflag = sum(1 for r in rows if r["flag"])
    print("\nsets assessed: %d ; flagged as high LD-contamination risk: %d" % (len(rows), nflag))
    print("WROTE", OUT)


if __name__ == "__main__":
    main()
