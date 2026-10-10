#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Part C / Step 2 — Pull the SNPs we need out of the UC GWAS harmonised file.

The harmonised file (de Lange 2017 UC; GCST004133) is ~10M rows / 299 MB gzip,
so we stream it exactly once and keep only the rsIDs we care about:

  * every cis-eQTL SNP of a candidate SGRG gene (for colocalisation), and
  * the LD-pruned instrument SNPs (subset of the above).

Matching is done on rsID (`hm_rsid`), which sidesteps the build mismatch between
GTEx v8 (GRCh38) and this GWAS (GRCh37).

Usage:  python partC_extract_gwas.py <rsid_list.txt> <out.tsv>
"""
import sys
import os
import csv
import gzip
import time

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
GWAS = os.path.join(ROOT, "data", "raw", "gwas", "GCST004133_uc_harmonised.h.tsv.gz")

KEEP = ["hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele",
        "hm_beta", "hm_effect_allele_frequency", "hm_code",
        "other_allele", "effect_allele", "beta", "standard_error", "p_value",
        "chromosome", "base_pair_location", "effect_allele_frequency"]


def main():
    rsid_file, out = sys.argv[1], sys.argv[2]
    want = set()
    with open(rsid_file, encoding="utf-8") as f:
        for line in f:
            s = line.strip().split("\t")[0]
            if s and s != "snpid":
                want.add(s)
    print("rsIDs wanted:", len(want), flush=True)

    n = 0
    hits = 0
    t0 = time.time()
    with gzip.open(GWAS, "rt", encoding="utf-8", errors="replace") as f, \
            open(out, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(KEEP)
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {c: hdr.index(c) for c in KEEP if c in hdr}
        missing = [c for c in KEEP if c not in hdr]
        if missing:
            print("NOTE header missing:", missing, flush=True)
        order = [c for c in KEEP if c in idx]
        for line in f:
            n += 1
            if n % 2000000 == 0:
                print("  scanned %d M rows, %d hits, %.0fs"
                      % (n // 1000000, hits, time.time() - t0), flush=True)
            i = line.find("\t")
            if i < 0:
                continue
            j = line.find("\t", i + 1)
            if j < 0:
                continue
            rs = line[i + 1:j]
            if rs in want:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < len(hdr):
                    continue
                w.writerow([parts[idx[c]] if c in idx else "" for c in order])
                hits += 1
    print("scanned %d rows, wrote %d hits -> %s  (%.0fs)"
          % (n, hits, out, time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
