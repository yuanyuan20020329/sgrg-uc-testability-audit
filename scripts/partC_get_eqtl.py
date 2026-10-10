#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Part C / Step 1 — Extract GTEx v8 cis-eQTLs for the 104 SGRG genes.

Why a region query: the GTEx API v2 gene-based endpoint
(`/association/singleTissueEqtl?gencodeId=...`) silently returns an empty set
(the documented pitfall). The *region* endpoint
(`/association/singleTissueEqtlByLocation`, TSS +/- 1 Mb) does return data and
carries everything we need for a cis-eQTL instrument:

    snpId       -> rsID
    variantId   -> "chr5_151481475_G_A_b38"  (GRCh38; REF and ALT alleles)
    pos         -> GRCh38 position
    nes         -> normalised effect size  == the beta we use as the instrument effect
    pValue      -> association p-value     == used to back out the SE
                   (GTEx tests the slope with a t-test; |z| = |nes/se| and
                    p = 2*Phi(-|z|)  =>  se = |nes| / |z|)
    geneSymbolUpper == the gene

The API is idempotent: raw JSON for every (gene, tissue) pair is cached under
data/raw/gtex_eqtl/, so re-running is nearly free.

Output: results/partC_eqtl_cis.tsv
"""
import os
import sys
import json
import time
import csv
import math
import urllib.parse
import urllib.request

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
CACHE = os.path.join(ROOT, "data", "raw", "gtex_eqtl")
OUT = os.path.join(ROOT, "results", "partC_eqtl_cis.tsv")
GENESET = os.path.join(ROOT, "results", "sgrg_gene_set_clean.tsv")
BASE = "https://gtexportal.org/api/v2"
TISSUES = ["Colon_Sigmoid", "Colon_Transverse", "Brain_Cortex", "Whole_Blood"]
WINDOW = 1_000_000          # +/- 1 Mb around TSS (standard cis window)
UA = "Mozilla/5.0 (compatible; research-script/1.0)"


def http_json(url, tries=4, timeout=90):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode("utf-8", "replace"))
        except Exception as e:
            if i == tries - 1:
                return {"__error__": "%s" % e}
            time.sleep(1.5 * (i + 1))
    return {}


def load_genes():
    syms = []
    with open(GENESET, encoding="utf-8-sig") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            s = (row.get("GeneSymbol") or "").strip()
            if s:
                syms.append(s)
    return sorted(set(syms))


def gene_coord(sym):
    d = http_json("%s/reference/gene?geneId=%s" % (BASE, urllib.parse.quote(sym)))
    if d.get("__error__") or not d.get("data"):
        return None
    g = d["data"][0]
    return {"symbol": g["geneSymbol"], "chrom": g["chromosome"], "tss": int(g["tss"]),
            "start": int(g["start"]), "end": int(g["end"]), "gid": g["gencodeId"],
            "strand": g.get("strand")}


def region_eqtl(chrom, start, end, tissue):
    url = ("%s/association/singleTissueEqtlByLocation?chromosome=%s&start=%d&end=%d"
           "&tissueSiteDetailId=%s" % (BASE, chrom, start, end, tissue))
    d = http_json(url)
    if d.get("__error__"):
        return None, d["__error__"]
    return d.get("singleTissueEqtl", []), None


def z_from_p(p):
    """two-sided p -> |z|"""
    if p is None:
        return None
    if p >= 1:
        return 0.0
    if p <= 0:
        return 40.0
    # inverse normal CDF via Peter Acklam's rational approximation
    return abs(_invnorm(p / 2.0))


def _invnorm(q):
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if q < plow:
        qq = math.sqrt(-2 * math.log(q))
        return (((((c[0]*qq+c[1])*qq+c[2])*qq+c[3])*qq+c[4])*qq+c[5]) / \
               ((((d[0]*qq+d[1])*qq+d[2])*qq+d[3])*qq+1)
    if q > phigh:
        qq = math.sqrt(-2 * math.log(1 - q))
        return -(((((c[0]*qq+c[1])*qq+c[2])*qq+c[3])*qq+c[4])*qq+c[5]) / \
               ((((d[0]*qq+d[1])*qq+d[2])*qq+d[3])*qq+1)
    qq = q - 0.5
    rr = qq * qq
    return (((((a[0]*rr+a[1])*rr+a[2])*rr+a[3])*rr+a[4])*rr+a[5])*qq / \
           (((((b[0]*rr+b[1])*rr+b[2])*rr+b[3])*rr+b[4])*rr+1)


def parse_variant(vid):
    """chr5_151481475_G_A_b38 -> (chrom, pos, ref, alt)"""
    parts = vid.split("_")
    if len(parts) < 5:
        return None
    return parts[0], int(parts[1]), parts[2].upper(), parts[3].upper()


def fetch_task(task):
    """task = (sym, gc, tissue). Returns list of parsed cis-eQTL rows."""
    sym, gc, tissue = task
    cache_f = os.path.join(CACHE, "%s__%s.json" % (sym, tissue))
    data = None
    if os.path.exists(cache_f) and os.path.getsize(cache_f) > 2:
        try:
            data = json.load(open(cache_f, encoding="utf-8"))
        except Exception:
            data = None
    if data is None:
        lo = max(0, gc["tss"] - WINDOW)
        hi = gc["tss"] + WINDOW
        arr, err = region_eqtl(gc["chrom"], lo, hi, tissue)
        if arr is None:
            return sym, tissue, None, err
        data = {"gene": sym, "tissue": tissue, "chrom": gc["chrom"],
                "tss": gc["tss"], "gid": gc["gid"], "window": [lo, hi],
                "rows": arr}
        tmp = cache_f + ".tmp"
        json.dump(data, open(tmp, "w", encoding="utf-8"))
        os.replace(tmp, cache_f)
    out = []
    for r in data["rows"]:
        if r.get("geneSymbolUpper") != sym.upper():
            continue
        pv = parse_variant(r.get("variantId", ""))
        if not pv:
            continue
        chrom, pos, ref, alt = pv
        p = r.get("pValue")
        nes = r.get("nes")
        z = z_from_p(p)
        se = (abs(nes) / z) if (z and z > 0 and nes is not None) else None
        out.append({
            "gene": sym, "gid": gc["gid"], "tissue": tissue,
            "snp": r.get("snpId"), "chrom": chrom, "pos": pos,
            "ref": ref, "alt": alt, "ea": alt, "oa": ref,
            "nes": nes, "se": se, "p": p, "snp_pos_api": r.get("pos"),
        })
    return sym, tissue, out, None


def main():
    from concurrent.futures import ThreadPoolExecutor, as_completed
    workers = int(os.environ.get("EQTL_WORKERS", "6"))
    os.makedirs(CACHE, exist_ok=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    genes = load_genes()
    print("SGRG genes:", len(genes), flush=True)

    t0 = time.time()
    coords = {}
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(gene_coord, s): s for s in genes}
        for fu in as_completed(futs):
            s = futs[fu]
            try:
                gc = fu.result()
            except Exception:
                gc = None
            if gc:
                coords[s] = gc
    print("coords resolved: %d/%d  [%.0fs]" % (len(coords), len(genes), time.time() - t0),
          flush=True)

    tasks = [(s, coords[s], t) for s in genes if s in coords for t in TISSUES]
    rows = []
    nerr = 0
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fetch_task, tk): tk for tk in tasks}
        for fu in as_completed(futs):
            tk = futs[fu]
            try:
                sym, tissue, out, err = fu.result()
            except Exception as e:
                nerr += 1
                print("   ERR %s" % e, flush=True)
                continue
            done += 1
            if err:
                nerr += 1
                print("   %s/%s ERR %s" % (sym, tissue, err), flush=True)
                continue
            rows.extend(out)
            if done % 25 == 0:
                print("  %d/%d (gene,tissue) done | %d cis rows | %.0fs"
                      % (done, len(tasks), len(rows), time.time() - t0), flush=True)
    print("fetched %d/%d tasks, %d errors, %d cis rows [%.0fs]"
          % (done, len(tasks), nerr, len(rows), time.time() - t0), flush=True)

    cols = ["gene", "gid", "tissue", "snp", "chrom", "pos", "ref", "alt",
            "ea", "oa", "nes", "se", "p", "snp_pos_api"]
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    print("WROTE %s (%d rows)" % (OUT, len(rows)), flush=True)

    # quick summary: genome-wide significant genes per tissue
    print("\n=== genome-wide significant (p<5e-8) cis-eQTL, per gene x tissue ===")
    from collections import defaultdict
    sig = defaultdict(list)
    for r in rows:
        try:
            if float(r["p"]) < 5e-8:
                sig[(r["gene"], r["tissue"])].append(float(r["p"]))
        except Exception:
            pass
    for (g, t), ps in sorted(sig.items()):
        print("  %-10s %-18s n=%3d  minP=%.2e" % (g, t, len(ps), min(ps)))


if __name__ == "__main__":
    main()
