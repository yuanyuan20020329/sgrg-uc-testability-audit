#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Part C / Stage 0-ext — cis-eQTL for the eIF2alpha-kinase pathway extension.

EIF2AK4 (GCN2) / EIF2AK2 (PKR) / EIF2AK3 (PERK) phosphorylate eIF2alpha, which is
the upstream trigger of stress-granule assembly.  They are NOT annotated to
GO:0010494 ("cytoplasmic stress granule", a cellular-component term) and are
therefore absent from the strict 104-gene SGRG set — but EIF2AK4 was the top
pre-specified candidate in the project's own pre-check report, so it is analysed
here as a clearly separated *pathway extension*, never mixed into the primary
SGRG panel.

Reuses the exact same fetch/parse machinery as partC_get_eqtl.py.
"""
import os
import sys
import csv
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import partC_get_eqtl as G

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT = os.path.join(ROOT, "results", "partC_eqtl_cis_ext.tsv")
EXT_GENES = ["EIF2AK4", "EIF2AK2", "EIF2AK3", "G3BP1"]   # G3BP1 = core marker, brain-only signal


def main():
    workers = 6
    genes = EXT_GENES
    coords = {}
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(G.gene_coord, s): s for s in genes}
        for fu in as_completed(futs):
            s = futs[fu]
            try:
                gc = fu.result()
            except Exception:
                gc = None
            if gc:
                coords[s] = gc
                print("coord ok:", s, gc["chrom"], gc["tss"], flush=True)
            else:
                print("coord FAIL:", s, flush=True)

    tasks = [(s, coords[s], t) for s in genes if s in coords for t in G.TISSUES]
    rows = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(G.fetch_task, tk): tk for tk in tasks}
        for fu in as_completed(futs):
            sym, tissue, out, err = fu.result()
            if err:
                print("ERR", sym, tissue, err, flush=True)
                continue
            rows.extend(out)
            sig = [r for r in out if float(r["p"]) < 5e-8]
            print("%-9s %-18s pairs=%4d sig=%3d%s"
                  % (sym, tissue, len(out), len(sig),
                     ("  minP=%.2e" % min(float(r["p"]) for r in sig)) if sig else ""),
                  flush=True)

    cols = ["gene", "gid", "tissue", "snp", "chrom", "pos", "ref", "alt",
            "ea", "oa", "nes", "se", "p", "snp_pos_api"]
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    print("WROTE", OUT, len(rows))


if __name__ == "__main__":
    main()
