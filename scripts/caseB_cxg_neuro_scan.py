# -*- coding: utf-8 -*-
"""案例 B：在 CELLxGENE 全库缓存（397 集合）中筛查神经退行相关集合及其来源结构。

L1（对照来源）审计的零成本入口：读元数据即可判定某集合的「对照/正常」细胞
是否为本研究原产（is_primary_data=True）还是复用外部数据（False）。

输出：
  results/caseB_cxg_neuro.tsv      —— 神经相关集合 × 数据集明细
  results/caseB_cxg_neuro.json     —— 聚合摘要
"""
import json
import glob
import os
import collections

NEURO_DISEASE = ["alzheimer", "parkinson", "dementia", "lewy", "amyotrophic",
                 "multiple sclerosis", "huntington", "frontotemporal", "neurodegener"]
BRAIN_TISSUE = ["brain", "cortex", "substantia nigra", "hippocamp", "midbrain",
                "striatum", "cerebell", "spinal cord", "dorsolateral prefrontal",
                "temporal lobe", "frontal", "neuron", "glia", "microglia"]
# CELLxGENE 的「对照」标注词
CTL = {"normal", "control", "healthy", "na", "unknown"}


def labels(field):
    """field 是 [{'label':..,'ontology_term_id':..}] 形式。"""
    out = []
    for f in field or []:
        if isinstance(f, dict):
            out.append((f.get("label") or "", f.get("ontology_term_id") or ""))
    return out


def main():
    rows = []
    for p in glob.glob("data/raw/cxg_coll_detail/*.json"):
        try:
            c = json.load(open(p, encoding="utf-8"))
        except Exception:
            continue
        cname = c.get("name") or ""
        cid = c.get("collection_id") or os.path.basename(p)
        for ds in c.get("datasets") or []:
            dis = labels(ds.get("disease"))
            tis = labels(ds.get("tissue"))
            disl = [d[0].lower() for d in dis]
            tisl = [t[0].lower() for t in tis]
            neuro = any(any(k in d for k in NEURO_DISEASE) for d in disl) or \
                    any(any(k in t for k in BRAIN_TISSUE) for t in tisl)
            if not neuro:
                continue
            cc = ds.get("cell_count") or 0
            pc = ds.get("primary_cell_count") or 0
            rows.append({
                "collection_id": cid,
                "collection_name": cname,
                "dataset_id": ds.get("dataset_id"),
                "title": (ds.get("title") or "")[:160],
                "assay": ";".join(a[0] for a in labels(ds.get("assay"))),
                "organism": ";".join(a[0] for a in labels(ds.get("organism"))),
                "disease": ";".join(d[0] for d in dis),
                "tissue": ";".join(t[0] for t in tis),
                "cell_count": cc,
                "primary_cell_count": pc,
                "frac_primary": round(pc / cc, 4) if cc else None,
                "is_primary_data": ds.get("is_primary_data"),
                "n_donors": len(ds.get("donor_id") or []),
                "doi": c.get("doi") or "",
                "collection_url": c.get("collection_url") or "",
                "raw_data_location": ds.get("raw_data_location") or "",
            })
        # end datasets
    # end files

    # write TSV
    cols = ["collection_id", "collection_name", "dataset_id", "title", "assay",
            "organism", "disease", "tissue", "cell_count", "primary_cell_count",
            "frac_primary", "is_primary_data", "n_donors", "doi", "raw_data_location"]
    with open("results/caseB_cxg_neuro.tsv", "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols) + "\n")

    # aggregate
    human = [r for r in rows if "homo sapiens" in (r["organism"] or "").lower()]
    agg = {
        "n_neuro_datasets": len(rows),
        "n_human": len(human),
        "n_collections": len(set(r["collection_id"] for r in rows)),
        "total_cells": sum(r["cell_count"] for r in human),
        "non_primary_cells": sum(r["cell_count"] - r["primary_cell_count"] for r in human),
        "disease_counts": dict(collections.Counter(
            r["disease"] for r in human).most_common(30)),
        "is_primary_data_counts": dict(collections.Counter(
            str(r["is_primary_data"]) for r in human).most_common()),
        "by_disease_nonprimary": {},
    }
    agg["share_non_primary_cells"] = round(
        agg["non_primary_cells"] / agg["total_cells"], 4) if agg["total_cells"] else None

    # per-disease: total cells, primary cells (L1 核心指标)
    dd = collections.defaultdict(lambda: [0, 0])
    for r in human:
        dd[r["disease"]][0] += r["cell_count"]
        dd[r["disease"]][1] += r["primary_cell_count"]
    agg["by_disease_nonprimary"] = {
        k: {"cells": v[0], "primary": v[1],
            "frac_primary": round(v[1] / v[0], 4) if v[0] else None}
        for k, v in sorted(dd.items(), key=lambda x: -x[1][0])[:40]
    }

    json.dump(agg, open("results/caseB_cxg_neuro.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)

    print("神经相关数据集: %d（人类 %d），涉及集合 %d 个" % (
        len(rows), len(human), agg["n_collections"]))
    print("人类细胞总数 %s，非原产 %s（%.2f%%）" % (
        agg["total_cells"], agg["non_primary_cells"],
        100 * (agg["share_non_primary_cells"] or 0)))
    print()
    print("is_primary_data 分布:", agg["is_primary_data_counts"])
    print()
    print("== 按疾病：细胞数 / frac_primary（L1 核心指标）==")
    for k, v in list(agg["by_disease_nonprimary"].items())[:25]:
        print("  %-45s cells=%-9d frac_primary=%.3f" % (k[:45], v["cells"], v["frac_primary"] or 0))
    print()
    print("== 神经相关集合清单（按细胞数） ==")
    seen = set()
    for r in sorted(human, key=lambda x: -x["cell_count"]):
        key = (r["collection_id"], r["disease"])
        if key in seen:
            continue
        seen.add(key)
        print("  %-42s | %-38s | cells=%-8d fp=%.2f | %s" % (
            r["collection_name"][:42], r["disease"][:38], r["cell_count"],
            r["frac_primary"] or 0, r["tissue"][:30]))
    print("\nwrote results/caseB_cxg_neuro.tsv / .json")


if __name__ == "__main__":
    main()
