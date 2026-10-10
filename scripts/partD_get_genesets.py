# -*- coding: utf-8 -*-
"""取基因集：MSigDB Hallmark + QuickGO 把 SGRG 面板拆成胞质/核应激颗粒两亚组"""
import json
import os
import ssl
import socket
import time
import urllib.request

socket.setdefaulttimeout(60)
CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE
os.makedirs("data/raw/genesets", exist_ok=True)


def fetch(url, binary=False):
    r = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(r, context=CTX) as f:
        b = f.read()
    return b if binary else b.decode("utf-8")


# ---- 1. Hallmark GMT ----
p = "data/raw/genesets/h.all.v2023.2.Hs.symbols.gmt"
if not os.path.exists(p):
    txt = fetch("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/"
                "h.all.v2023.2.Hs.symbols.gmt")
    open(p, "w", encoding="utf-8").write(txt)
print("[ok] hallmark:", p, os.path.getsize(p), "bytes")

# ---- 2. QuickGO：两个 GO term 的 human 基因 ----
out = {}
for go, name in [("GO:0010494", "cytoplasmic_stress_granule"),
                 ("GO:0035061", "nuclear_stress_granule")]:
    genes, page = set(), 1
    while page <= 60:
        url = ("https://www.ebi.ac.uk/QuickGO/services/annotation/search"
               "?goId=%s&taxonId=9606&geneProductType=protein&limit=100&page=%d" % (go, page))
        j = json.loads(fetch(url))
        res = j.get("results", [])
        for r in res:
            sym = r.get("symbol")
            if sym:
                genes.add(sym)
        total_pages = (j.get("pageInfo") or {}).get("total", 1)
        if not res or page >= total_pages:
            break
        page += 1
        time.sleep(0.15)
        if page % 10 == 0:
            print("   ... %s page %d / %d, 累计 %d 基因" % (go, page, total_pages, len(genes)))
    out[name] = sorted(genes)
    print("[ok] %s %s -> %d 基因" % (go, name, len(genes)))

json.dump(out, open("results/partD_go_sgrg_subsets.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)

# 与 HGNC 面板对照
panel = [l.split("\t")[0] for l in
         open("results/sgrg_gene_set_clean.tsv", encoding="utf-8").read().splitlines()[1:] if l.strip()]
cyto = set(out["cytoplasmic_stress_granule"]) & set(panel)
nuc = set(out["nuclear_stress_granule"]) & set(panel)
print("\n面板 %d 基因中：胞质 SG %d、核 SG %d、交集 %d、两者皆无 %d"
      % (len(panel), len(cyto), len(nuc), len(cyto & nuc), len(set(panel) - cyto - nuc)))
print("胞质:", sorted(cyto))
print("核:", sorted(nuc))
print("皆无:", sorted(set(panel) - cyto - nuc))
