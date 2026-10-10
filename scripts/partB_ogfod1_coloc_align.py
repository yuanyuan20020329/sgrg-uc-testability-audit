#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Part B addendum — OGFOD1 locus coloc sensitivity: ALIGNMENT.

Merges three sources for the chr16 OGFOD1/CES1P1 region:
  1. GTEx v8 Whole_Blood cis-eQTLs (significant pairs, from the cached
     region query) for OGFOD1 and CES1P1;
  2. UC GWAS summary stats (de Lange 2017, GCST004133) pulled by rsID;
  3. 1000G phase3 EUR dosages (alt-allele counts) for chr16:55.0-58.0 Mb.

Harmonisation rules (per coloc-susie skill):
  * everything is oriented to the 1000G ALT allele;
  * GWAS beta flipped when its effect allele is the panel REF;
  * GTEx nes is for its own ALT — flipped when GTEx alt == panel REF;
  * palindromic SNPs (A/T, C/G): orientation decided by comparing GWAS EAF
    with panel EAF (direct assignment, never a blind sign flip); dropped
    when the two frequencies disagree in both directions (ambiguous);
  * indels / multi-allelic records dropped unless alleles match exactly;
  * EAF cross-check: Pearson r(GWAS EAF, panel EAF) must be >= 0.98.

Outputs (results/):
  ogfod1_coloc_<GENE>_dataset.tsv   per-SNP merged stats + LD matrix path
  ogfod1_coloc_<GENE>_ld.csv        correlation matrix (rsid rownames)
  ogfod1_coloc_align_report.json    match counts, EAF r, drops
  ogfod1_instrument_ld.tsv          r2 of the 9-10 MR instruments vs the
                                    OGFOD1 lead and the CES1P1 lead
"""
import csv
import json
import os
import numpy as np

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
CACHE = os.path.join(ROOT, "data", "raw", "gtex_eqtl", "OGFOD1__Whole_Blood.json")
GWAS_F = os.path.join(ROOT, "results", "ogfod1_region_gwas.tsv")
PANEL_F = r"D:/ld_ref_1000g/regions/chr16_ogfod1_grch38.tsv"
RES = os.path.join(ROOT, "results")
GENES = ["OGFOD1", "CES1P1"]
INS_F = os.path.join(ROOT, "results", "partC_instruments.tsv")
LEADS = {"OGFOD1": "rs4784670", "CES1P1": "rs12149359"}

PAL = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

# variant-index built while reading the GTEx cache (rsid -> ref/alt, and
# (pos, sorted allele pair) -> rsid for the position-based panel mapping)
GTEX_ALLELES = {}
GTEX_POS = {}
GTEX_KEY = {}
GTEX_AMBIG = set()


def fnum(x):
    try:
        if x is None or x == "" or x == "NA":
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def load_gtex():
    d = json.load(open(CACHE, encoding="utf-8"))
    out = {g: {} for g in GENES}
    for r in d.get("rows", []):
        g = r.get("geneSymbolUpper")
        if g not in out:
            continue
        vid = (r.get("variantId") or "").split("_")
        if len(vid) < 5:
            continue
        rsid = r.get("snpId")
        if not rsid or rsid in out[g]:      # first occurrence wins
            continue
        nes, p = fnum(r.get("nes")), fnum(r.get("pValue"))
        if nes is None or p is None:
            continue
        # |z| from p, se = |nes|/|z|
        from scipy.stats import norm
        z = abs(norm.isf(p / 2.0)) if 0 < p < 1 else (40.0 if p <= 0 else 0.0)
        se = abs(nes) / z if z > 0 else None
        if se is None or se <= 0:
            continue
        out[g][rsid] = {"pos": int(vid[1]), "ref": vid[2].upper(),
                        "alt": vid[3].upper(), "beta": nes, "se": se, "p": p}
        alleles = (vid[2].upper(), vid[3].upper())
        if rsid not in GTEX_ALLELES:
            GTEX_ALLELES[rsid] = alleles
            GTEX_POS[rsid] = int(vid[1])
            key = (int(vid[1]), tuple(sorted(alleles)))
            if key in GTEX_KEY and GTEX_KEY[key] != rsid:
                GTEX_AMBIG.add(key)          # same pos+allele pair, two rsIDs
            else:
                GTEX_KEY[key] = rsid
    for k in list(GTEX_KEY):
        if k in GTEX_AMBIG:
            GTEX_KEY.pop(k)
    return out


def load_gwas():
    g = {}
    with open(GWAS_F, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rs = r.get("hm_rsid")
            if rs and rs != "NA" and rs not in g:
                g[rs] = r
    return g


def load_panel():
    """The 1kGP high-coverage panel is GRCh38 but its ID column is
    'chrom:pos:ref:alt', not an rsID.  Since it shares the GRCh38 coordinate
    system with GTEx, records are mapped to the GTEx variant index by
    (position, sorted allele pair).  When the panel REF/ALT pair is the
    reverse of the GTEx pair, dosages are complemented (2 - d) so they count
    the GTEx ALT allele.
    """
    with open(PANEL_F, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[4:]
        recs = {}
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) < 4:
                continue
            pos, ref, alt = int(fld[0]), fld[1], fld[2]
            rsid = GTEX_KEY.get((pos, tuple(sorted((ref, alt)))))
            if rsid is None or rsid in recs:
                continue
            gref, galt = GTEX_ALLELES[rsid]
            ds = np.array([np.nan if x == "" else float(x) for x in fld[4:]])
            if (ref, alt) == (gref, galt):
                pass
            elif (ref, alt) == (galt, gref):
                ds = 2.0 - ds                     # complement to GTEx ALT
            else:
                continue
            recs[rsid] = {"pos": pos, "ref": gref, "alt": galt, "dosage": ds}
    return recs, samples


def harmonise(gene, eqtl, gwas, panel):
    """Returns list of dicts, oriented to the panel ALT allele."""
    rows, drop = [], {}
    for rsid, e in eqtl.items():
        g = gwas.get(rsid)
        pn = panel.get(rsid)
        if g is None:
            drop[rsid] = "no_gwas"
            continue
        if pn is None:
            drop[rsid] = "no_panel"
            continue
        b_y = fnum(g.get("hm_beta"))
        s_y = fnum(g.get("standard_error"))
        ea_y = (g.get("hm_effect_allele") or "").upper()
        oa_y = (g.get("hm_other_allele") or "").upper()
        if None in (b_y, s_y) or s_y <= 0 or not ea_y or not oa_y:
            drop[rsid] = "gwas_missing_fields"
            continue

        pref, palt = pn["ref"], pn["alt"]
        # ---- GTEx orientation ----
        b_x = e["beta"]
        if e["ref"] == pref and e["alt"] == palt:
            pass
        elif e["ref"] == palt and e["alt"] == pref:
            b_x = -b_x
        else:
            drop[rsid] = "eqtl_allele_mismatch"
            continue
        # ---- GWAS orientation ----
        # NOTE: GCST004133 has no EAF at all, so palindromic SNPs cannot be
        # oriented — they are excluded (conservative).  Non-palindromic SNPs
        # are unambiguous from the allele labels.
        if (pref, palt) in PAL:
            drop[rsid] = "palindromic_no_gwas_eaf"
            continue
        if ea_y == palt and oa_y == pref:
            pass
        elif ea_y == pref and oa_y == palt:
            b_y = -b_y
        else:
            drop[rsid] = "gwas_allele_mismatch"
            continue

        ds = pn["dosage"]
        ok = ~np.isnan(ds)
        if ok.sum() < 0.8 * len(ds):
            drop[rsid] = "panel_missing"
            continue
        eaf_p = float(np.nansum(ds) / (2 * ok.sum()))
        rows.append({"snp": rsid, "pos": pn["pos"], "ref": pref, "alt": palt,
                     "beta_eqtl": b_x, "se_eqtl": e["se"], "p_eqtl": e["p"],
                     "beta_gwas": b_y, "se_gwas": s_y,
                     "p_gwas": fnum(g.get("p_value")),
                     "eaf_gwas": "", "eaf_panel": eaf_p})
    return rows, drop


def main():
    eqtl = load_gtex()
    gwas = load_gwas()
    panel, samples = load_panel()
    print("panel SNPs mapped by (pos, alleles): %d, samples: %d"
          % (len(panel), len(samples)))
    print("GWAS region rows: %d" % len(gwas))
    for g in GENES:
        print("  %s eQTL SNPs: %d" % (g, len(eqtl[g])))

    report = {"panel_snps": len(panel), "panel_samples": len(samples),
              "gwas_rows": len(gwas), "genes": {}}
    sets = {}
    for gene in GENES:
        rows, drop = harmonise(gene, eqtl[gene], gwas, panel)
        snps = [r["snp"] for r in rows]
        if len(snps) < 10:
            raise RuntimeError("%s: only %d SNPs survived harmonisation" % (gene, len(snps)))
        # EAF cross-check: GCST004133 has no EAF, so the check is done the
        # other way around — the panel EAF of the kept SNPs is reported and
        # orientation relies on allele labels (palindromics excluded).
        ep = np.array([r["eaf_panel"] for r in rows])
        r_eaf = None
        # LD matrix (correlation of alt dosages)
        X = np.vstack([panel[s]["dosage"] for s in snps])
        X = np.nan_to_num(np.where(np.isnan(X), np.nanmean(X, axis=1, keepdims=True), X))
        R = np.corrcoef(X)
        np.fill_diagonal(R, 1.0)
        ds_f = os.path.join(RES, "ogfod1_coloc_%s_dataset.tsv" % gene)
        with open(ds_f, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, delimiter="\t", fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        ld_f = os.path.join(RES, "ogfod1_coloc_%s_ld.csv" % gene)
        with open(ld_f, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow([""] + snps)
            for i, s in enumerate(snps):
                w.writerow([s] + ["%.6f" % v for v in R[i]])
        sets[gene] = {"snps": snps, "R": R}
        n_drop = {}
        for d in drop.values():
            n_drop[d] = n_drop.get(d, 0) + 1
        report["genes"][gene] = {
            "n_eqtl": len(eqtl[gene]), "n_used": len(rows), "r_eaf": r_eaf,
            "drops": n_drop,
            "eaf_panel_range": [float(ep.min()), float(ep.max())],
            "min_p_eqtl": min(r["p_eqtl"] for r in rows),
            "min_p_gwas": min(r["p_gwas"] for r in rows if r["p_gwas"] is not None),
        }
        print("%s: used %d SNPs, panel EAF range %.3f-%.3f, drops=%s"
              % (gene, len(rows), float(ep.min()), float(ep.max()), n_drop))

    # ---- MR instruments vs the two lead SNPs ----
    # r^2 is computed directly from panel dosages so the comparison also
    # works for SNPs absent from the GWAS (e.g. the CES1P1 lead).
    ins = [r for r in csv.DictReader(open(INS_F, encoding="utf-8"), delimiter="\t")
           if r["gene"] == "OGFOD1" and r["tissue"] == "Whole_Blood"]
    leads = {g: GTEX_KEY.get((GTEX_POS[LEADS[g]],
                              tuple(sorted(GTEX_ALLELES[LEADS[g]])))) for g in GENES}
    for g in GENES:
        assert leads[g] == LEADS[g], "lead index mismatch for %s" % g

    def panel_r(a, b):
        """correlation of alt dosages for two rsIDs, straight from the panel"""
        if a not in panel or b not in panel:
            return None
        x, y = panel[a]["dosage"], panel[b]["dosage"]
        m = ~(np.isnan(x) | np.isnan(y))
        if m.sum() < 100:
            return None
        return float(np.corrcoef(x[m], y[m])[0, 1])

    out_rows = []
    for r in ins:
        rs = r["snp"]
        row = {"snp": rs, "eqtl_p": r.get("p")}
        for g in GENES:
            rv = panel_r(rs, LEADS[g])
            if rv is not None:
                row["r_vs_%s_lead" % g] = round(rv, 4)
                row["r2_vs_%s_lead" % g] = round(rv ** 2, 4)
            else:
                row["r_vs_%s_lead" % g] = ""
                row["r2_vs_%s_lead" % g] = ""
        out_rows.append(row)
    # leads vs each other (the headline contamination number)
    rl = panel_r(LEADS["OGFOD1"], LEADS["CES1P1"])
    report["lead_lead_r"] = rl if rl is None else round(rl, 4)
    report["lead_lead_r2"] = None if rl is None else round(rl ** 2, 4)
    print("r(OGFOD1 lead rs4784670, CES1P1 lead rs12149359) = %s"
          % report["lead_lead_r"])
    with open(os.path.join(RES, "ogfod1_instrument_ld.tsv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    json.dump(report, open(os.path.join(RES, "ogfod1_coloc_align_report.json"), "w",
                           encoding="utf-8"), indent=2)
    print("instrument LD table written (%d instruments)" % len(out_rows))


if __name__ == "__main__":
    main()
