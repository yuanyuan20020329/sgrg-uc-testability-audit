# -*- coding: utf-8 -*-
"""诊断 2：物种构成、样本标签可得性、可用的「后果」指标覆盖率。"""
import collections
import glob
import json
import os
import re

CACHE = r"C:/Users/H-P/WorkBuddy/应激颗粒/data/raw/gemma_a2"
idx = json.load(open(os.path.join(CACHE, "cohort_index.json"), encoding="utf-8"))


def taxon_of(e):
    t = e.get("taxon") or {}
    return t.get("commonName") if isinstance(t, dict) else str(t)


print("[1] 列表条目字段：", list(next(iter(idx.values())).keys()))
tax = collections.Counter()
for k, e in idx.items():
    tax[taxon_of(e)] += 1
print("[2] taxon 分布：", dict(tax))

print("[3] 人类数据集里，有多少能算出「UC vs 对照」主轴可分性 / 有疾病 DE：")
DISEASE = {"disease", "disease state", "disease status", "diagnosis", "clinical diagnosis", "condition", "disease outcome"}
UC = re.compile(r"ulcerative colitis|\buc\b", re.I)
CTRL = re.compile(r"control|normal|healthy|unaffected|non-inflam|reference", re.I)

n_human = n_auc = n_de = 0
rows = []
for k, e in idx.items():
    if taxon_of(e) != "human":
        continue
    n_human += 1
    did = k
    # 样本标签
    lab = collections.defaultdict(set)
    for p in glob.glob(os.path.join(CACHE, "samp_%s_*.json" % did)):
        d = json.load(open(p, encoding="utf-8"))
        for ent in (d.get("data") or []):
            s = ent.get("sample") or ent
            for c in (s.get("characteristics") or []):
                cat = (c.get("category") or "").strip().lower()
                lab[cat].add(str(c.get("value")))
    dis_vals = sorted(v for c in lab if c in DISEASE for v in lab[c])
    uc = [v for v in dis_vals if UC.search(v) and "non-involved" not in v.lower()]
    ct = [v for v in dis_vals if CTRL.search(v) and not UC.search(v)]
    # DE
    de = json.load(open(os.path.join(CACHE, "an_%s.json" % did), encoding="utf-8"))
    has_dis_de = False
    for a in (de.get("data") or []):
        for fid, fvs in (a.get("factorValuesUsed") or {}).items():
            cats = {(c.get("category") or "").lower() for fv in fvs for c in (fv.get("characteristics") or [])}
            if cats & DISEASE and a.get("resultSets"):
                has_dis_de = True
    rows.append((did, e.get("shortName"), len(dis_vals), bool(uc), bool(ct), has_dis_de,
                 ";".join(lab.get("disease", set()))[:60]))
    n_auc += bool(uc and ct)
    n_de += has_dis_de

print("  人类数据集 %d 个；有 UC+对照标签 %d；有疾病 DE 分析 %d" % (n_human, n_auc, n_de))
print("\n  有 UC/对照标签的人类数据集：")
for r in rows:
    if r[3] and r[4]:
        print("   gid=%-7s %-10s 疾病标签数=%d DE=%s  %s" % (r[0], r[1], r[2], r[5], r[6]))

print("\n  人类数据集里既无 UC 也无对照标签的前 15 个：")
for r in [x for x in rows if not (x[3] and x[4])][:15]:
    print("   gid=%-7s %-10s 标签=%d DE=%s  %s" % (r[0], r[1], r[2], r[5], r[6]))
