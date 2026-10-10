# -*- coding: utf-8 -*-
"""
Fetch a GRCh38 1kGP high-coverage EUR panel for the OGFOD1/CES1P1 region.

Why: the phase3 20130502 release is GRCh37 AND has '.' in its ID column, while
GTEx v8 variant IDs are GRCh38 — no reachable liftOver service exists in this
sandbox.  The 20220422 3202-sample NYGC callset is already GRCh38 and carries
rsIDs, so positions/alleles match GTEx directly.

Region: chr16:55,300,000-57,600,000 (GRCh38) — the OGFOD1 cis window
(55,451,490-57,451,490) plus padding.

Output: D:/ld_ref_1000g/regions/chr16_ogfod1_grch38.tsv
        pos, ref, alt, rsid + EUR dosages (alt-allele counts)
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ld_bgzf_tabix import read_tbi, region_file_range, decompress_range
import urllib.request

BASE = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
        "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
        "1kGP_high_coverage_Illumina.chr16.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
TBI = r"D:/ld_ref_1000g/chr16_1kgp38.tbi"
POPS = r"D:/ld_ref_1000g/g1k_3202_pops.txt"
OUT = r"D:/ld_ref_1000g/regions/chr16_ogfod1_grch38.tsv"
START, END = 55_300_000, 57_600_000
HEAD_PROBE_BYTES = 400_000


def fetch_range(url, start, stop, timeout=3600):
    req = urllib.request.Request(url, headers={"Range": "bytes=%d-%d" % (start, stop)})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def main():
    tbi = read_tbi(TBI)
    names = tbi.get("names", [])
    ref_name = "chr16" if "chr16" in names else ("16" if "16" in names else names[0])
    rng = region_file_range(tbi, ref_name, START - 1, END)
    mb = (rng[1] - rng[0]) / 1e6
    print("contig name: %s | byte range: %s (%.1f MB)" % (ref_name, rng, mb))

    # EUR sample list
    eur = []
    with open(POPS, encoding="utf-8") as f:
        f.readline()
        for line in f:
            fld = line.strip().split()          # file is whitespace-separated
            if len(fld) >= 7 and fld[6] == "EUR":
                eur.append(fld[1])
    print("EUR samples:", len(eur))
    if not eur:
        raise RuntimeError("EUR sample list is empty (pops file parsing?)")

    # header (#CHROM) — the region chunk lacks it
    head = decompress_range(fetch_range(BASE, 0, HEAD_PROBE_BYTES), 0).decode("utf-8", "replace")
    samples = None
    for line in head.split("\n"):
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            break
    if not samples:
        raise RuntimeError("no #CHROM line in head probe")
    idx = [samples.index(x) for x in eur if x in samples]
    print("matched EUR in VCF:", len(idx))
    if not idx:
        raise RuntimeError("no EUR sample matched the VCF header")

    t0 = time.time()
    data = fetch_range(BASE, rng[0], rng[1] - 1)
    print("downloaded %.1f MB in %.0fs" % (len(data) / 1e6, time.time() - t0))
    txt = decompress_range(data, 0).decode("utf-8", "replace")

    eur_used = [samples[i] for i in idx]
    n_kept = 0
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        fh.write("\t".join(["pos", "ref", "alt", "rsid"] + eur_used) + "\n")
        for line in txt.split("\n"):
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != len(samples) + 9:
                continue
            try:
                pos = int(f[1])
            except ValueError:
                continue
            if pos < START or pos > END:
                continue
            alt = f[4]
            if "," in alt:
                continue
            # GT is the first subfield of FORMAT in this callset
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
            ds = [dosage(f[9 + i]) for i in idx]
            ok = sum(1 for d in ds if d is not None)
            if ok < 0.8 * len(ds):
                continue
            af = sum(d for d in ds if d is not None) / (2 * ok)
            if af < 0.01 or af > 0.99:
                continue
            fh.write("\t".join([str(pos), f[3].upper(), alt.upper(), f[2]] +
                               ["" if d is None else str(d) for d in ds]) + "\n")
            n_kept += 1
    print("kept %d variants -> %s" % (n_kept, OUT))


if __name__ == "__main__":
    main()
