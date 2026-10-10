"""Europe PMC bibliometric calibration for a meta-research / resource-audit paper.

Purpose: find real published precedents for 'audit of public omics resources'
and report which journals actually host them, so the user can calibrate a
realistic submission target instead of guessing.

Usage: python epmc_meta_research_calib.py
"""
import json
import time
import urllib.parse
import urllib.request
from collections import Counter

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

QUERIES = {
    "A_sc_dataset_quality_audit": (
        '(TITLE:"single-cell" OR TITLE:"single cell") '
        'AND (TITLE:"dataset" OR TITLE:"datasets" OR TITLE:"data") '
        'AND (TITLE:"confounding" OR TITLE:"confounded" OR TITLE:"batch" '
        'OR TITLE:"pitfall" OR TITLE:"pitfalls" OR TITLE:"quality" '
        'OR TITLE:"reproducibility")'
    ),
    "B_public_data_reuse_audit": (
        '(TITLE:"public" AND (TITLE:"datasets" OR TITLE:"data")) '
        'AND (TITLE:"reuse" OR TITLE:"reusability" OR TITLE:"reanalysis" '
        'OR TITLE:"reproducibility" OR TITLE:"integrity" OR TITLE:"duplicate")'
    ),
    "C_meta_research_omics": (
        '("meta-research" OR "meta research" OR "metaresearch") '
        'AND (ABSTRACT:"single-cell" OR ABSTRACT:"genomics" '
        'OR ABSTRACT:"bioinformatics" OR ABSTRACT:"transcriptomics" '
        'OR ABSTRACT:"public datasets")'
    ),
    "D_uc_scrna_public": (
        '(ABSTRACT:"ulcerative colitis" OR ABSTRACT:"inflammatory bowel disease") '
        'AND (ABSTRACT:"single-cell RNA" OR ABSTRACT:"scRNA-seq") '
        'AND (ABSTRACT:"publicly available" OR ABSTRACT:"public datasets" '
        'OR ABSTRACT:"reanalysis" OR ABSTRACT:"integrated")'
    ),
    "E_depth_confounding": (
        '(TITLE:"sequencing depth" OR TITLE:"library size" OR TITLE:"cell number") '
        'AND (TITLE:"confounding" OR TITLE:"confounder" OR TITLE:"artefact" '
        'OR TITLE:"artifact" OR TITLE:"spurious")'
    ),
}


def fetch(query, page_size=100):
    params = {
        "query": query,
        "format": "json",
        "pageSize": str(page_size),
        "resultType": "lite",
    }
    url = BASE + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "meta-research-calib/1.0"})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return json.loads(fh.read().decode("utf-8"))


def main():
    out = {}
    for name, q in QUERIES.items():
        try:
            data = fetch(q)
        except Exception as exc:  # noqa: BLE001
            print(f"[{name}] FAILED: {exc}")
            continue
        hits = data.get("hitCount")
        results = data.get("resultList", {}).get("result", [])
        journals = Counter(
            (r.get("journalTitle") or "NO_JOURNAL") for r in results
        )
        rows = []
        for r in results:
            rows.append(
                {
                    "title": (r.get("title") or "")[:150],
                    "journal": r.get("journalTitle") or "",
                    "year": r.get("pubYear") or "",
                    "cites": r.get("citedByCount") or 0,
                    "doi": r.get("doi") or "",
                    "pmid": r.get("pmid") or "",
                }
            )
        rows.sort(key=lambda x: -x["cites"])
        out[name] = {"query": q, "hitCount": hits, "journals": journals, "rows": rows}
        print("=" * 78)
        print(f"[{name}]  hitCount = {hits}  (fetched {len(results)})")
        print(f"  query: {q}")
        print("  top journals:")
        for j, c in journals.most_common(12):
            print(f"    {c:>3}  {j}")
        print("  most cited papers in window:")
        for r in rows[:12]:
            print(f"    {r['cites']:>5}  {r['year']}  {r['journal'][:38]:<38} {r['title'][:70]}")
        time.sleep(0.4)
    with open("epmc_meta_research_calib.json", "w", encoding="utf-8") as fh:
        json.dump(
            {k: {"query": v["query"], "hitCount": v["hitCount"], "rows": v["rows"]}
             for k, v in out.items()},
            fh,
            ensure_ascii=False,
            indent=1,
        )
    print("\nwrote epmc_meta_research_calib.json")


if __name__ == "__main__":
    main()
