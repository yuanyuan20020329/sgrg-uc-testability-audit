# -*- coding: utf-8 -*-
"""
案例 A · Step 0b：在 Gemma 中系统筛选「肿瘤 vs 癌旁正常」肝癌队列。

输出：results/caseA_gemma_cohorts.tsv
评价维度：
  n_total, n_tumor, n_adjacent, 患者配对数, batchEffect, 平台, 表达数据可下（HEAD 探测）
"""
import sys, os, json, collections, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from caseA_gemma import list_datasets, get_dataset, sample_labels, CACHE

OUT = r"C:\Users\H-P\WorkBuddy\应激颗粒\results\caseA_gemma_cohorts.tsv"

TUMOR_KEYS = ("tumor", "tumour", "carcinoma", "cancer", "hcc", "neoplasm")
NORM_KEYS = ("adjacent", "non-tumor", "nontumor", "non-tumour", "normal", "peritumor",
             "paratumor", "cirrhotic", "healthy")

USEFUL_BATCH = ("NO_BATCH_EFFECT_SUCCESS", "BATCH_CORRECTED_SUCCESS", "SINGLE_BATCH_SUCCESS",
                "PROBLEMATIC_BATCH_INFO_FAILURE", "BATCH_EFFECT_FAILURE", "SINGLETON_BATCHES_FAILURE",
                "NO_BATCH_INFO")


def classify(label):
    if not label:
        return None
    l = label.lower()
    # 先判正常（"adjacent non-tumor liver" 同时含 tumor 字样）
    if any(k in l for k in NORM_KEYS):
        return "normal"
    if any(k in l for k in TUMOR_KEYS):
        return "tumor"
    return None


def main():
    print("listing HCC datasets ...", flush=True)
    ds = list_datasets("hepatocellular carcinoma", limit=100)
    print("  total listed:", len(ds), flush=True)

    human = [d for d in ds if (d.get("taxon") or {}).get("commonName") == "human"
             or (d.get("taxon") or {}).get("scientificName") == "Homo sapiens"]
    print("  human:", len(human), flush=True)

    # 预筛：n 在 20..1000，且 batchEffect 是可用枚举
    cands = []
    for d in human:
        n = d.get("bioAssayCount") or 0
        be = d.get("batchEffect")
        if 20 <= n <= 1000 and be in USEFUL_BATCH:
            cands.append(d)
    print("  after size/batch prefilter:", len(cands), flush=True)

    rows = []
    for i, d in enumerate(cands):
        did = d.get("id")
        name = d.get("shortName")
        try:
            srows = sample_labels(did)
        except Exception as e:
            print("  %s sample ERR %s" % (name, repr(e)[:80]), flush=True)
            continue
        grp = collections.Counter(classify(r["label"]) for r in srows)
        n_t = grp.get("tumor", 0)
        n_n = grp.get("normal", 0)
        if n_t < 10 or n_n < 10:
            continue
        # 患者配对：同名 patient_key 同时出现在 tumor 与 normal
        byk = collections.defaultdict(set)
        for r in srows:
            g = classify(r["label"])
            if g:
                byk[r["patient_key"]].add(g)
        paired = sum(1 for k, v in byk.items() if len(v) == 2)
        plat = collections.Counter(r["platform"] for r in srows)
        rows.append({
            "dataset": name, "id": did, "n_total": len(srows),
            "n_tumor": n_t, "n_normal": n_n, "n_paired_patients": paired,
            "batchEffect": d.get("batchEffect"),
            "platforms": ";".join("%s:%d" % (k, v) for k, v in plat.most_common()),
            "name": (d.get("name") or "")[:70],
        })
        print("  [%d/%d] %-13s n=%4d (T%4d/N%4d) paired=%3d %s" % (
            i + 1, len(cands), name, len(srows), n_t, n_n, paired, d.get("batchEffect")), flush=True)

    rows.sort(key=lambda r: (-r["n_paired_patients"], -r["n_total"]))
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT, "rows:", len(rows))


if __name__ == "__main__":
    main()
