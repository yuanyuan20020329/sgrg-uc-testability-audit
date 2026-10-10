# -*- coding: utf-8 -*-
"""
Part C LD audit — re-clump every instrument set on measured LD, then re-run IVW.

Part C selected instruments by physical distance (>10 kb apart).  That is not LD
clumping: the 10 OGFOD1 instruments, spread over ~160 kb, turned out to sit in
one LD block (effective independent instruments 1.09), so a null lead-SNP Wald
(p = 0.118) became a spurious IVW p = 1.7e-7.

This script repeats the check for every gene x tissue row in Part C:
  * effective number of independent instruments (eigenvalue method)
  * greedy LD clumping at r^2 < 0.001 / 0.01 / 0.1 (GTEx eQTL p ascending)
  * IVW re-estimated on the clumped set (single SNP -> Wald ratio)
  * allele harmonisation against the 1kGP GRCh38 EUR panel (oriented to panel ALT);
    palindromic SNPs are dropped because GCST004133 harmonised carries no EAF

Inputs : D:/ld_ref_1000g/audit/{gene}.tsv      (from partC_ld_audit_fetch.py)
         results/partC_instruments.tsv
         results/partC_audit_gwas.tsv
         results/partC_mr_geneset.tsv          (nominal IVW, for comparison)
Output : results/partC_ld_audit.tsv
"""
import collections
import csv
import math
import os

import numpy as np
from scipy import stats

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
AUDIT = r"D:/ld_ref_1000g/audit"
INS_F = os.path.join(ROOT, "results", "partC_instruments.tsv")
GWAS_F = os.path.join(ROOT, "results", "partC_audit_gwas.tsv")
MR_F = os.path.join(ROOT, "results", "partC_mr_geneset.tsv")
OUT_F = os.path.join(ROOT, "results", "partC_ld_audit.tsv")
PAL = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
THRESH = (("0.001", 0.001), ("0.01", 0.01), ("0.1", 0.1))


def fnum(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def load_panel(gene):
    path = os.path.join(AUDIT, "%s.tsv" % gene)
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        n_sam = len(header) - 4
        recs = {}
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) < 4:
                continue
            try:
                pos = int(fld[0])
            except ValueError:
                continue
            recs[pos] = {"ref": fld[1].upper(), "alt": fld[2].upper(),
                         "dosage": np.array([np.nan if x == "" else float(x)
                                             for x in fld[4:4 + n_sam]])}
    return recs


def load_gwas():
    out = {}
    with open(GWAS_F, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rs = r.get("hm_rsid")
            if not rs or rs == "NA":
                continue
            b = fnum(r.get("hm_beta"))
            if b is None:
                b = fnum(r.get("beta"))
            se = fnum(r.get("standard_error"))
            if b is None or se is None or se <= 0:
                continue
            out[rs] = {"pos": int(r["hm_pos"]),
                       "ea": (r["hm_effect_allele"] or "").upper(),
                       "oa": (r["hm_other_allele"] or "").upper(),
                       "b": b, "se": se, "p": fnum(r.get("p_value"))}
    return out


def clump(snps, r2map, pvals, thr):
    """Greedy clumping: SNP order by ascending eQTL p; keep if r2 <= thr vs all kept."""
    kept = []
    for s in sorted(snps, key=lambda x: pvals[x]):
        if all(r2map.get((s, k), 0.0) <= thr for k in kept):
            kept.append(s)
    return kept


def ivw(bx, by, sy):
    w = 1.0 / (sy ** 2)
    den = float(np.sum(w * bx ** 2))
    b = float(np.sum(w * bx * by)) / den
    se = math.sqrt(1.0 / den)
    return b, se, 2 * stats.norm.sf(abs(b / se))


def main():
    gwas = load_gwas()
    print("GWAS rows usable: %d" % len(gwas), flush=True)

    nominal = {}
    with open(MR_F, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            nominal[(r["gene"], r["tissue"])] = fnum(r.get("IVW_p"))

    by_gt = collections.defaultdict(list)
    with open(INS_F, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            by_gt[(r["gene"], r["tissue"])].append(r)

    panels, rows = {}, []
    notes = collections.Counter()

    for (gene, tissue) in sorted(by_gt):
        ins = by_gt[(gene, tissue)]
        if gene not in panels:
            panels[gene] = load_panel(gene)
        panel = panels[gene]
        if panel is None:
            notes["panel_missing:%s" % gene] += 1
            continue

        # ---------- model-matrix build (all instruments found in the panel) ----------
        snps, X, pval = [], [], {}
        dropped = collections.Counter()
        for r in ins:
            rs, pos = r["snp"], int(r["pos"])
            p = fnum(r.get("p"))
            if pos not in panel:
                dropped["not_in_panel"] += 1
                continue
            if p is None:
                dropped["no_eqtl_p"] += 1
                continue
            snps.append(rs)
            X.append(panel[pos]["dosage"])
            pval[rs] = p
        if not snps:
            notes["no_snp_%s/%s" % (gene, tissue)] += 1
            continue

        X = np.vstack(X)
        n = len(snps)
        R = np.eye(n)
        for i in range(n):
            for j in range(i + 1, n):
                xi, xj = X[i], X[j]
                m = ~np.isnan(xi) & ~np.isnan(xj)
                if m.sum() < 30 or np.nanstd(xi[m]) == 0 or np.nanstd(xj[m]) == 0:
                    c = 0.0
                else:
                    c = float(np.corrcoef(xi[m], xj[m])[0, 1])
                R[i, j] = R[j, i] = c
        R2 = R ** 2
        r2map = {(snps[i], snps[j]): float(R2[i, j])
                 for i in range(n) for j in range(n) if i != j}
        off = R2[np.triu_indices(n, 1)] if n > 1 else np.array([0.0])
        ev = np.linalg.eigvalsh(R)
        ev = ev[ev > 1e-9]
        eff = float((ev.sum() ** 2) / (ev ** 2).sum()) if len(ev) else float("nan")

        # ---------- allele harmonisation, oriented to panel ALT ----------
        ali = {}
        for r in ins:
            rs, pos = r["snp"], int(r["pos"])
            if rs not in pval or pos not in panel:
                continue
            pn, g = panel[pos], gwas.get(rs)
            ref_, alt_ = pn["ref"], pn["alt"]
            if (ref_, alt_) in PAL:
                dropped["palindromic"] += 1
                continue
            if g is None:
                dropped["gwas_missing"] += 1
                continue
            ea, oa = (r["ea"] or "").upper(), (r["oa"] or "").upper()
            nes, se_x = fnum(r.get("nes")), fnum(r.get("se"))
            if nes is None or se_x is None or se_x <= 0:
                dropped["bad_eqtl"] += 1
                continue
            if {ea, oa} != {ref_, alt_}:
                dropped["eqtl_allele_mismatch"] += 1
                continue
            if {g["ea"], g["oa"]} != {ref_, alt_}:
                dropped["gwas_allele_mismatch"] += 1
                continue
            ali[rs] = {"bx": nes if ea == alt_ else -nes, "sx": se_x,
                       "by": g["b"] if g["ea"] == alt_ else -g["b"],
                       "sy": g["se"], "p": pval[rs]}

        res = {"gene": gene, "tissue": tissue, "n_nominal": n, "n_aligned": len(ali),
               "eff_indep": round(eff, 2),
               "median_r2": round(float(np.median(off)), 4) if len(off) else "",
               "max_r2": round(float(off.max()), 4) if len(off) else "",
               "drop": "|".join("%s:%d" % kv for kv in sorted(dropped.items())),
               "IVW_p_nominal": nominal.get((gene, tissue), "")}

        ok = list(ali)
        for tag, thr in THRESH:
            keep = clump(ok, r2map, pval, thr)
            res["n_clump_r2_%s" % tag] = len(keep)
            if not keep:
                continue
            if len(keep) == 1:
                d = ali[keep[0]]
                b = d["by"] / d["bx"]
                se = d["sy"] / abs(d["bx"])
                p = 2 * stats.norm.sf(abs(b / se))
            else:
                b, se, p = ivw(np.array([ali[s]["bx"] for s in keep]),
                               np.array([ali[s]["by"] for s in keep]),
                               np.array([ali[s]["sy"] for s in keep]))
            res["IVW_b_r2_%s" % tag] = round(float(b), 4)
            res["IVW_se_r2_%s" % tag] = round(float(se), 4)
            res["IVW_p_r2_%s" % tag] = float("%.4g" % p)
        rows.append(res)

    # ---------- BH-FDR on the conventional clumping (r2 < 0.001) ----------
    valid = sorted([(i, r["IVW_p_r2_0.001"]) for i, r in enumerate(rows)
                    if isinstance(r.get("IVW_p_r2_0.001"), float)], key=lambda x: x[1])
    m = len(valid)
    for rank, (i, p) in enumerate(valid, 1):
        rows[i]["FDR_r2_0.001"] = round(min(1.0, p * m / rank), 6)

    cols = ["gene", "tissue", "n_nominal", "n_aligned", "eff_indep", "median_r2", "max_r2",
            "n_clump_r2_0.001", "n_clump_r2_0.01", "n_clump_r2_0.1",
            "IVW_p_nominal", "IVW_p_r2_0.001", "IVW_p_r2_0.01", "IVW_p_r2_0.1",
            "IVW_b_r2_0.001", "IVW_se_r2_0.001", "FDR_r2_0.001", "drop"]
    with open(OUT_F, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x.get("IVW_p_r2_0.001")
                                             if isinstance(x.get("IVW_p_r2_0.001"), float) else 1)):
            w.writerow(r)

    key = lambda x: (x.get("IVW_p_r2_0.001") if isinstance(x.get("IVW_p_r2_0.001"), float) else 1)
    print("\nrows: %d | notes: %s\n" % (len(rows), dict(notes)))
    print("%-9s %-16s %4s %4s %5s %6s %7s %11s %11s" %
          ("gene", "tissue", "nom", "clp", "effIV", "medR2", "maxR2", "p_nominal", "p_clumped"))
    for r in sorted(rows, key=key):
        print("%-9s %-16s %4d %4s %5s %6s %7s %11s %11s" % (
            r["gene"], r["tissue"], r["n_nominal"], r.get("n_clump_r2_0.001", ""),
            r["eff_indep"], r["median_r2"], r["max_r2"],
            ("%.3g" % r["IVW_p_nominal"]) if isinstance(r["IVW_p_nominal"], float) else "-",
            ("%.3g" % r["IVW_p_r2_0.001"]) if isinstance(r.get("IVW_p_r2_0.001"), float) else "-"))
    print("\n-> %s" % OUT_F)


if __name__ == "__main__":
    main()
