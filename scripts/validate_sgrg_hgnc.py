# -*- coding: utf-8 -*-
"""
用 HGNC 官方 REST API 校验 SGRG 基因集符号（GO:0010494 U GO:0035061, Homo sapiens）。
目的：EBI QuickGO 返回的 `symbol` 字段混杂了 UniProt accession（如 O00361、A8K5C4）
      与计算推断条目，直接用作基因面板会虚高分母。此处仅保留 HGNC 批准的官方符号。
输出：results/sgrg_gene_set_clean.tsv  (symbol, evidence_codes, hgnc_id, ensembl_id, source_go)
      results/sgrg_rejected.tsv        (被剔除的符号及原因)
"""
import json, os, re, glob, urllib.request, urllib.parse
import concurrent.futures as cf

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
RAW = os.path.join(BASE, "data", "raw")
OUT = os.path.join(BASE, "results")
os.makedirs(OUT, exist_ok=True)

# ---- 1. 汇总 GO 注释中的 symbol -> evidence ----
sym_ev = {}
for f in glob.glob(os.path.join(RAW, "qgo", "*.json")):
    try:
        d = json.load(open(f, encoding="utf-8"))
    except Exception:
        continue
    for r in d.get("results", []) or []:
        s = (r.get("symbol") or "").strip()
        if not s:
            continue
        ev = r.get("goEvidence") or r.get("evidenceCode") or "?"
        sym_ev.setdefault(s, set()).add(ev)

print(f"[input] unique symbols from QuickGO = {len(sym_ev)}")

# ---- 2. 逐个查 HGNC ----
def hgnc_lookup(sym):
    url = "https://rest.genenames.org/fetch/symbol/" + urllib.parse.quote(sym)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=25) as r:
            d = json.load(r)
        resp = d.get("response", {})
        if resp.get("numFound", 0) > 0:
            doc = resp["docs"][0]
            return sym, True, doc.get("hgnc_id", ""), doc.get("ensembl_gene_id", "")
        return sym, False, "", ""
    except Exception as e:
        return sym, None, f"ERR:{type(e).__name__}", ""

syms = sorted(sym_ev)
results = []
with cf.ThreadPoolExecutor(max_workers=6) as ex:
    for i, res in enumerate(ex.map(hgnc_lookup, syms), 1):
        results.append(res)
        if i % 40 == 0:
            print(f"  ... {i}/{len(syms)}", flush=True)

ok, bad, err = [], [], []
for sym, found, hid, ens in results:
    if found is True:
        ok.append((sym, hid, ens))
    elif found is False:
        bad.append(sym)
    else:
        err.append(sym)

print(f"[HGNC] approved = {len(ok)} | rejected = {len(bad)} | lookup-error = {len(err)}")
if err:
    print("  errors:", err)

# ---- 3. 输出 ----
p1 = os.path.join(OUT, "sgrg_gene_set_clean.tsv")
with open(p1, "w", encoding="utf-8") as f:
    f.write("GeneSymbol\tHR_evidence\tHGNC_id\tEnsembl_gene_id\n")
    for sym, hid, ens in ok:
        ev = ",".join(sorted(sym_ev[sym]))
        f.write(f"{sym}\t{ev}\t{hid}\t{ens}\n")
print("[write]", p1)

p2 = os.path.join(OUT, "sgrg_rejected.tsv")
with open(p2, "w", encoding="utf-8") as f:
    f.write("RejectedSymbol\tReason\n")
    for s in sorted(bad):
        reason = "not_an_HGNC_approved_symbol"
        if re.match(r"^[A-Z][0-9]", s) and len(s) >= 5:
            reason = "likely_UniProt_accession"
        elif s.startswith("hCG_") or s.startswith("DKFZp"):
            reason = "non-standard_clone/ORF_name"
        elif "_" in s:
            reason = "non-standard_name"
        f.write(f"{s}\t{reason}\n")
print("[write]", p2)

summary = {
    "go_terms": ["GO:0010494 (cytoplasmic stress granule)", "GO:0035061 (nuclear stress granule)"],
    "taxon": "Homo sapiens (NCBITaxon:9606)",
    "symbols_from_quickgo": len(sym_ev),
    "hgnc_approved": len(ok),
    "rejected": len(bad),
    "approved_symbols": [s for s, _, _ in ok],
    "rejected_symbols": sorted(bad),
}
json.dump(summary, open(os.path.join(OUT, "sgrg_gene_set_clean.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=2)
print("[done]")
