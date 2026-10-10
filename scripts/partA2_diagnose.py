# -*- coding: utf-8 -*-
"""诊断：为什么多数队列取不到「疾病 vs 对照」的存算 DE / 算不出主轴可分性。"""
import collections
import glob
import json
import os

CACHE = r"C:/Users/H-P/WorkBuddy/应激颗粒/data/raw/gemma_a2"

cat_by_factor = collections.Counter()
datasets_no_dis = []
datasets_ok = []
for p in glob.glob(os.path.join(CACHE, "an_*.json")):
    did = os.path.basename(p)[3:-5]
    d = json.load(open(p, encoding="utf-8"))
    if not d.get("data"):
        datasets_no_dis.append((did, "NO_ANALYSIS"))
        continue
    found = False
    for a in d["data"]:
        for fid, fvs in (a.get("factorValuesUsed") or {}).items():
            cats = set()
            for fv in fvs:
                for c in (fv.get("characteristics") or []):
                    cats.add((c.get("category") or "?").strip().lower())
                    cat_by_factor[(c.get("category") or "?").strip().lower()] += 1
            if cats & {"disease", "disease state", "disease status", "diagnosis", "clinical diagnosis",
                       "condition", "disease outcome"}:
                found = True
                datasets_ok.append(did)
                break
        if found:
            break
    if not found:
        datasets_no_dis.append((did, "NO_DISEASE_FACTOR"))

print("有疾病因子分析的队列：%d" % len(set(datasets_ok)))
print("没有的：%d" % len(datasets_no_dis))
print("\n各特征类别出现次数（top 30）：")
for k, v in cat_by_factor.most_common(30):
    print("  %-32s %d" % (k, v))

print("\n抽样 12 个「无疾病因子」的数据集，看它们的因子值长什么样：")
shown = 0
for did, why in datasets_no_dis:
    if why != "NO_DISEASE_FACTOR" or shown >= 12:
        continue
    d = json.load(open(os.path.join(CACHE, "an_%s.json" % did), encoding="utf-8"))
    a = d["data"][0]
    desc = []
    for fid, fvs in (a.get("factorValuesUsed") or {}).items():
        cats = sorted({(c.get("category") or "?") for fv in fvs for c in (fv.get("characteristics") or [])})
        vals = sorted({str(fv.get("description") or fv.get("factorValue") or "")[:28] for fv in fvs})[:4]
        desc.append("fid=%s cats=%s n=%d vals=%s" % (fid, cats, len(fvs), vals))
    print("  gid=%-7s name=%-22s %s" % (did, a.get("name"), " || ".join(desc)[:260]))
    shown += 1

print("\n没有存算分析的数据集数：%d" % sum(1 for _, w in datasets_no_dis if w == "NO_ANALYSIS"))
