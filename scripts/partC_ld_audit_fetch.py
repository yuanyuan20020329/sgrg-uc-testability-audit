# -*- coding: utf-8 -*-
"""
Part C instrument LD audit — fetch 1kGP GRCh38 EUR genotypes for every gene region.

Motivation
----------
Part C selected MR instruments by *physical distance* (>10 kb apart) instead of
LD clumping.  For OGFOD1 this produced 10 "independent" instruments that are in
fact a single LD block (effective independent instruments = 1.09): the lead-SNP
Wald ratio is null (p = 0.118) while naively pooling the 10 correlated SNPs gave
IVW p = 1.7e-7.  Part C also reports MR-positive results for 9 genes whose
instrument counts span 1-47 SNPs, and *none* of them passed colocalisation
(PP.H4 max 0.108 across 46 gene x tissue rows) - i.e. the two methods inside
Part C already disagreed.

This script pulls the same external LD reference for all 35 Part C gene regions
so every instrument set can be re-clumped on measured LD.

Reference panel: 1kGP 20220422 high-coverage phased SNV/INDEL callset
(3202 samples, GRCh38, tabix-indexed on EBI FTP), restricted to EUR (n=633).
GRCh38 positions match GTEx v8 / GWAS Catalog harmonised coordinates directly,
so no liftOver is needed (none is available in this sandbox).

Output
------
D:/ld_ref_1000g/audit/{gene}.tsv
    pos, ref, alt, vcf_id, <EUR sample columns>   (alt-allele dosages, 0/1/2)
Only variants requested by results/partC_instruments.tsv are written.
D:/ld_ref_1000g/audit/_fetch_log.json
    per-gene byte ranges, downloaded sizes, kept/missing variants, failures.
"""
import collections
import csv
import json
import os
import sys
import time
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ld_bgzf_tabix import read_tbi, region_file_range, decompress_range

BASE_DIR = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
            "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/")
TBI_DIR = r"D:/ld_ref_1000g"
POPS = os.path.join(TBI_DIR, "g1k_3202_pops.txt")
OUTDIR = os.path.join(TBI_DIR, "audit")
SAMPLES_CACHE = os.path.join(TBI_DIR, "panel_samples.json")
INS_F = r"C:/Users/H-P/WorkBuddy/应激颗粒/results/partC_instruments.tsv"
PAD = 250_000
HEAD_PROBE_BYTES = 4_000_000          # 400 kB missed #CHROM on some contigs
TIMEOUT = 1800


def panel_name(chrom):
    """chrX is published as ...phased_panel.v2.vcf.gz; all others have no .v2."""
    base = "1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel" % chrom
    if chrom == "X":
        base += ".v2"
    return base


def panel_url(chrom):
    return BASE_DIR + panel_name(chrom) + ".vcf.gz"


def fetch_range(url, start, stop, timeout=TIMEOUT, tries=6):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=%d-%d" % (start, stop)})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:                      # noqa: BLE001
            last = e
            time.sleep(5 + 10 * i)
    raise RuntimeError("range fetch failed: %s" % last)


def ensure_tbi(chrom):
    p = os.path.join(TBI_DIR, "chr%s_1kgp38.tbi" % chrom)
    if os.path.exists(p) and os.path.getsize(p) > 10_000:
        return p
    url = panel_url(chrom) + ".tbi"
    print("  [tbi] downloading chr%s index ..." % chrom, flush=True)
    data = fetch_range(url, 0, 50_000_000)
    with open(p, "wb") as fh:
        fh.write(data)
    print("  [tbi] chr%s -> %.1f MB" % (chrom, len(data) / 1e6), flush=True)
    return p


def panel_samples(url):
    """Column-10+ sample IDs from the #CHROM line (cached once for all contigs)."""
    if os.path.exists(SAMPLES_CACHE):
        with open(SAMPLES_CACHE, encoding="utf-8") as fh:
            return json.load(fh)
    head = decompress_range(fetch_range(url, 0, HEAD_PROBE_BYTES), 0).decode("utf-8", "replace")
    samples = None
    for line in head.split("\n"):
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            break
    if not samples:
        raise RuntimeError("no #CHROM in %d-byte head probe" % HEAD_PROBE_BYTES)
    with open(SAMPLES_CACHE, "w", encoding="utf-8") as fh:
        json.dump(samples, fh)
    return samples


def load_eur():
    eur = []
    with open(POPS, encoding="utf-8") as f:
        f.readline()
        for line in f:
            fld = line.strip().split()          # whitespace-separated (CRLF)
            if len(fld) >= 7 and fld[6] == "EUR":
                eur.append(fld[1])
    if not eur:
        raise RuntimeError("EUR sample list empty (pops file parsing)")
    return eur


def dosage(gt):
    if not gt or gt[0] == ".":
        return None
    a = gt.split(":")[0].replace("|", "/").split("/")
    if len(a) != 2 or "." in a:
        return None
    try:
        return int(a[0]) + int(a[1])
    except ValueError:
        return None


def norm_chrom(c):
    return str(c).replace("chr", "").replace("CHR", "")


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    eur = load_eur()
    print("EUR samples in pops file: %d" % len(eur), flush=True)

    # ---- group instruments by gene ----
    by_gene = collections.defaultdict(list)
    with open(INS_F, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            by_gene[r["gene"]].append(r)

    regions = {}
    for g, rs in sorted(by_gene.items()):
        chrom = collections.Counter(norm_chrom(r["chrom"]) for r in rs).most_common(1)[0][0]
        ps = [int(r["pos"]) for r in rs]
        regions[g] = {"chrom": chrom, "start": max(1, min(ps) - PAD),
                      "end": max(ps) + PAD,
                      "snps": sorted({r["snp"] for r in rs})}
    print("gene regions to fetch: %d" % len(regions), flush=True)

    log = {"panel": "1kGP 20220422 high-coverage phased SNV_INDEL_SV (GRCh38)",
           "eur_n": len(eur), "pad_bp": PAD, "genes": {}, "failed": []}

    for g, reg in sorted(regions.items()):
        chrom, start, end = reg["chrom"], reg["start"], reg["end"]
        outp = os.path.join(OUTDIR, "%s.tsv" % g)
        if os.path.exists(outp):
            print("[skip] %-9s already fetched" % g, flush=True)
            continue
        try:
            t0 = time.time()
            tbi = read_tbi(ensure_tbi(chrom))
            names = tbi.get("names", [])
            ref_name = chrom if chrom in names else ("chr%s" % chrom if "chr%s" % chrom in names
                                                     else names[0])
            rng = region_file_range(tbi, ref_name, start - 1, end)
            url = panel_url(chrom)

            samples = panel_samples(url)
            idx = [samples.index(x) for x in eur if x in samples]
            if not idx:
                raise RuntimeError("no EUR sample matched")

            data = fetch_range(url, rng[0], rng[1] - 1)
            txt = decompress_range(data, 0).decode("utf-8", "replace")

            want_pos = {int(r["pos"]) for r in by_gene[g]}
            kept, seen = 0, set()
            with open(outp, "w", newline="", encoding="utf-8") as fh:
                fh.write("\t".join(["pos", "ref", "alt", "vcf_id"] +
                                   [samples[i] for i in idx]) + "\n")
                for line in txt.split("\n"):
                    if not line or line.startswith("##") or line.startswith("#CHROM"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) != len(samples) + 9:
                        continue
                    try:
                        pos = int(f[1])
                    except ValueError:
                        continue
                    if pos not in want_pos or "," in f[4]:
                        continue
                    ds = [dosage(f[9 + i]) for i in idx]
                    ok = sum(1 for d in ds if d is not None)
                    if ok < 0.8 * len(ds):
                        continue
                    seen.add(pos)
                    fh.write("\t".join([str(pos), f[3].upper(), f[4].upper(), f[2]] +
                                       ["" if d is None else str(d) for d in ds]) + "\n")
                    kept += 1

            missing = sorted(want_pos - seen)
            log["genes"][g] = {"chrom": chrom, "start": start, "end": end,
                               "byte_range": rng, "download_mb": round(len(data) / 1e6, 2),
                               "eur_matched": len(idx), "kept": kept,
                               "requested": len(want_pos), "missing_pos": missing,
                               "seconds": round(time.time() - t0, 1)}
            print("[ok]   %-9s chr%-2s %d-%d | %.1f MB | kept %d/%d | miss %d | %.0fs"
                  % (g, chrom, start, end, len(data) / 1e6, kept, len(want_pos),
                     len(missing), time.time() - t0), flush=True)
        except Exception as e:                      # noqa: BLE001
            log["failed"].append({"gene": g, "chrom": chrom, "error": str(e)})
            print("[FAIL] %-9s %s" % (g, e), flush=True)

    with open(os.path.join(OUTDIR, "_fetch_log.json"), "w", encoding="utf-8") as fh:
        json.dump(log, fh, ensure_ascii=False, indent=2)
    print("\nfailed genes: %s" % [x["gene"] for x in log["failed"]], flush=True)
    print("log -> %s" % os.path.join(OUTDIR, "_fetch_log.json"), flush=True)


if __name__ == "__main__":
    main()
