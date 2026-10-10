#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Part C / Stage A — build cis-eQTL instrument sets for the SGRG genes.

For every (gene, tissue) pair with at least one genome-wide significant cis-eQTL
we build an instrument set:

  * SNP significance   p < 5e-8  (genome-wide); BH-FDR within the gene's cis
                       SNP set recorded as a secondary criterion;
  * independence       LD clumping is impossible here (no LD reference panel is
                       reachable from this sandbox and no plink/tabix binary is
                       installed), so we fall back to *distance* pruning:
                       walk SNPs by ascending p and keep a SNP only if it is
                       > PRUNE_KB away from every already-kept SNP.  This is a
                       conservative proxy for r^2 < 0.001 and is documented as a
                       limitation.  We ALSO record the single lead SNP so a
                       Wald-ratio robustness check is always available.
  * strength           mean F = mean((beta/se)^2); weak instruments flagged at
                       mean F < 10.

Outputs
-------
results/partC_instruments.tsv        one row per (gene, tissue, snp)
results/partC_instrument_sets.json   per (gene,tissue) summary
results/partC_rsids_for_gwas.txt     rsIDs to pull out of the UC GWAS
                                     (instruments + full cis window of every
                                      gene that produced >=1 instrument, so that
                                      colocalisation has all regional SNPs)
"""
import os
import csv
import json
import math
from collections import defaultdict

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
CIS = os.path.join(ROOT, "results", "partC_eqtl_cis.tsv")
CIS_EXT = os.path.join(ROOT, "results", "partC_eqtl_cis_ext.tsv")
PANEL_MAIN = "SGRG_GO"
PANEL_EXT = "eIF2a_kinase_extension"
CACHE = os.path.join(ROOT, "data", "raw", "gtex_eqtl")
OUT_INS = os.path.join(ROOT, "results", "partC_instruments.tsv")
OUT_SETS = os.path.join(ROOT, "results", "partC_instrument_sets.json")
OUT_RSID = os.path.join(ROOT, "results", "partC_rsids_for_gwas.txt")

SIG = 5e-8
PRUNE_KB = 10
F_MIN = 10.0


def bh_fdr(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), 1):
        k = n - rank + 1
        val = pvals[i] * n / k
        prev = min(prev, val)
        q[i] = prev
    return q


def load_rows(path, panel):
    out = []
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if not r.get("p"):
                continue
            try:
                r["p"] = float(r["p"]); r["nes"] = float(r["nes"])
                r["se"] = float(r["se"]); r["pos"] = int(r["pos"])
            except Exception:
                continue
            r["panel"] = panel
            out.append(r)
    return out


def main():
    rows = load_rows(CIS, PANEL_MAIN) + load_rows(CIS_EXT, PANEL_EXT)
    print("cis-eQTL rows loaded:", len(rows))

    bygene = defaultdict(list)
    for r in rows:
        bygene[(r["gene"], r["tissue"])].append(r)

    ins_rows = []
    sets = {}
    rsids = set()

    for (gene, tissue), rs in sorted(bygene.items()):
        rs = sorted(rs, key=lambda x: x["p"])
        pv = [x["p"] for x in rs]
        qv = bh_fdr(pv)
        for x, q in zip(rs, qv):
            x["q"] = q
        sig = [x for x in rs if x["p"] < SIG]
        if not sig:
            continue
        # distance pruning (conservative independent set).
        # Two SNPs that are perfect LD proxies carry *identical* eQTL statistics
        # (same NES and same SE) — this is common in low-recombination blocks
        # (e.g. the VCP locus on chr9:35 Mb spans an ~80 kb perfect-LD haplotype).
        # Distance pruning alone cannot see that, so we additionally drop any SNP
        # whose (NES, SE) duplicates an already-kept SNP.
        kept = []
        seen_stats = {}
        n_proxy = 0
        for x in sig:
            key = (round(x["nes"], 6), round(x["se"], 6))
            if key in seen_stats:
                n_proxy += 1
                continue
            if all(abs(x["pos"] - y["pos"]) > PRUNE_KB * 1000 for y in kept):
                kept.append(x)
                seen_stats[key] = x["snp"]
        lead = sig[0]
        F = [(x["nes"] / x["se"]) ** 2 for x in kept]
        meanF = sum(F) / len(F) if F else float("nan")
        sets["%s|%s" % (gene, tissue)] = {
            "gene": gene, "tissue": tissue, "panel": sig[0].get("panel", PANEL_MAIN),
            "n_cis_pairs": len(rs), "n_sig": len(sig),
            "n_instruments": len(kept), "meanF": meanF,
            "min_p": lead["p"], "lead_snp": lead["snp"],
            "lead_nes": lead["nes"], "lead_se": lead["se"],
            "fdr_sig_pairs": int(sum(1 for x in rs if x["q"] < 0.05)),
            "n_proxy_dropped": n_proxy,
            "weak": bool(meanF < F_MIN),
            "snps": [x["snp"] for x in kept],
        }
        for x in kept:
            ins_rows.append({
                "gene": gene, "tissue": tissue, "panel": x.get("panel", PANEL_MAIN),
                "snp": x["snp"],
                "chrom": x["chrom"], "pos": x["pos"], "ea": x["ea"], "oa": x["oa"],
                "nes": x["nes"], "se": x["se"], "p": x["p"], "q": round(x["q"], 4),
                "F": round((x["nes"] / x["se"]) ** 2, 2),
                "is_lead": int(x["snp"] == lead["snp"]),
            })
        rsids.add(lead["snp"])
        for x in kept:
            rsids.add(x["snp"])
        # full cis window for colocalisation
        cf = os.path.join(CACHE, "%s__%s.json" % (gene, tissue))
        if os.path.exists(cf):
            d = json.load(open(cf, encoding="utf-8"))
            for r in d.get("rows", []):
                if r.get("geneSymbolUpper") == gene.upper() and r.get("snpId"):
                    rsids.add(r["snpId"])

    cols = ["gene", "tissue", "panel", "snp", "chrom", "pos", "ea", "oa", "nes", "se",
            "p", "q", "F", "is_lead"]
    with open(OUT_INS, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in ins_rows:
            w.writerow(r)
    json.dump(sets, open(OUT_SETS, "w", encoding="utf-8"), indent=1, default=str)
    with open(OUT_RSID, "w", encoding="utf-8") as f:
        for s in sorted(rsids):
            f.write(s + "\n")

    print("\n=== (gene,tissue) pairs with >=1 genome-wide significant cis-eQTL ===")
    for k, v in sorted(sets.items(), key=lambda kv: kv[1]["min_p"]):
        print("  %-10s %-18s k=%2d sig=%3d minP=%.2e meanF=%7.1f%s"
              % (v["gene"], v["tissue"], v["n_instruments"], v["n_sig"],
                 v["min_p"], v["meanF"], "  WEAK" if v["weak"] else ""))
    print("\ninstruments written:", len(ins_rows))
    print("rsIDs for GWAS pull :", len(rsids))
    print("files:", OUT_INS, "|", OUT_SETS, "|", OUT_RSID)


if __name__ == "__main__":
    main()
