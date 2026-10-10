# -*- coding: utf-8 -*-
"""把 Gemma 铁死亡表达矩阵与样本标签拼成 analysis-ready 长表。

输出：
  results/caseA_expr_7673.tsv         行=样本, 列=sample,label,gene symbols...（宽表）
  results/caseA_sample_labels_7673.tsv 行=样本, sample,label,gsm,gemma_name,platform
"""
import sys, os, json, csv, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from caseA_gemma import all_samples

ROOT = r"C:\Users\H-P\WorkBuddy\应激颗粒"
EXPR = os.path.join(ROOT, "data", "raw", "caseA", "gemma", "ferroptosis_expr_7673.json")
OUTL = os.path.join(ROOT, "results", "caseA_sample_labels_7673.tsv")
OUTE = os.path.join(ROOT, "results", "caseA_expr_7673.tsv")
DID = 7673


def main():
    d = json.load(open(EXPR, encoding="utf-8"))
    gel = d["data"][0]["geneExpressionLevels"]
    # sample order from first gene vector
    sample_order = list(gel[0]["vectors"][0]["bioAssayExpressionLevels"].keys())
    # gene x sample matrix
    mat = {}
    for g in gel:
        sym = g.get("geneOfficialSymbol")
        vecs = g.get("vectors") or []
        if not vecs:
            continue
        # consolidate=average -> one vector
        vals = vecs[0]["bioAssayExpressionLevels"]
        mat[sym] = [vals.get(s) for s in sample_order]

    rows = all_samples(DID)
    meta = {}
    for r in rows:
        smp = r.get("sample") or {}
        ch = {c["category"]: c["value"] for c in (smp.get("characteristics") or [])}
        meta[r.get("name")] = {
            "label": ch.get("tissue") or ch.get("BioSource") or "?",
            "gsm": (r.get("accession") or {}).get("accession"),
            "gemma_name": smp.get("name"),
            "platform": (r.get("arrayDesign") or {}).get("shortName"),
            "ajcc": ch.get("ajcc stage"),
        }

    matched = sum(1 for s in sample_order if s in meta)
    print("samples in expr: %d ; matched to metadata: %d" % (len(sample_order), matched))

    with open(OUTL, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "label", "gsm", "gemma_name", "platform", "ajcc"])
        for s in sample_order:
            m = meta.get(s, {})
            w.writerow([s, m.get("label", "?"), m.get("gsm", ""), m.get("gemma_name", ""),
                        m.get("platform", ""), m.get("ajcc", "")])

    genes = [g for g in mat if g]
    with open(OUTE, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "label"] + genes)
        for i, s in enumerate(sample_order):
            m = meta.get(s, {})
            w.writerow([s, m.get("label", "?")] + [mat[g][i] for g in genes])

    print("wrote", OUTL)
    print("wrote", OUTE, "genes:", len(genes), "samples:", len(sample_order))
    cnt = collections.Counter(meta.get(s, {}).get("label") for s in sample_order)
    print("label counts:", dict(cnt))


if __name__ == "__main__":
    main()
