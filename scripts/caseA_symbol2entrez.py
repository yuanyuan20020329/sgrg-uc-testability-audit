# -*- coding: utf-8 -*-
"""
symbol -> NCBI Entrez ID 映射（用 mygene.info，本沙箱可达）。
Gemma 的 /expressions/genes/ 端点只吃 Entrez ID 或 Ensembl ID，不吃 symbol。
输出：results/caseA_ferroptosis_entrez.tsv
"""
import urllib.request, urllib.parse, json, time, os

GENES = r"C:\Users\H-P\WorkBuddy\应激颗粒\results\caseA_ferroptosis_genes.txt"
OUT = r"C:\Users\H-P\WorkBuddy\应激颗粒\results\caseA_ferroptosis_entrez.tsv"


def get(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            return json.loads(urllib.request.urlopen(req, timeout=60).read())
        except Exception as e:
            last = e
            time.sleep(1.5 * (i + 1))
    raise last


def main():
    genes = [l.strip() for l in open(GENES, encoding="utf-8") if l.strip()]
    rows, missing = [], []
    for i, g in enumerate(genes):
        q = urllib.parse.quote(g)
        url = ("https://mygene.info/v3/query?q=symbol:%s&species=human"
               "&fields=symbol,name,entrezgene,ensembl.gene&size=5" % q)
        try:
            d = get(url)
            hit = None
            for h in d.get("hits", []):
                if str(h.get("symbol", "")).upper() == g.upper() and h.get("entrezgene"):
                    hit = h
                    break
            if hit is None and d.get("hits"):
                h = d["hits"][0]
                if h.get("entrezgene"):
                    hit = h
            if hit:
                eg = str(hit["entrezgene"])
                rows.append((g, eg, hit.get("symbol", ""), hit.get("name", "")))
            else:
                missing.append(g)
        except Exception as e:
            missing.append("%s(ERR %s)" % (g, repr(e)[:50]))
        if (i + 1) % 15 == 0:
            print("  %d/%d mapped" % (i + 1, len(genes)), flush=True)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("query_symbol\tentrez_id\tmatched_symbol\tname\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    print("mapped %d/%d -> %s" % (len(rows), len(genes), OUT))
    if missing:
        print("unmapped:", missing)


if __name__ == "__main__":
    main()
