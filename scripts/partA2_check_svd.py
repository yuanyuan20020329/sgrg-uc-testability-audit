# -*- coding: utf-8 -*-
"""校验：能否用 Gemma 的 SVD 载荷（不入库表达矩阵）复现「PC1 疾病可分性」。

阳性对照 = GSE38713（5873）：我们此前用原始矩阵算得 PC1 单独预测疾病的 AUC = 0.92、
仅批次 AUC = 0.77（批次与疾病 Cramér's V = 0.71）。若 SVD 载荷能复现，方法成立。
"""
import json
import os

import numpy as np

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
svd = json.load(open(os.path.join(ROOT, "data/raw/gemma_a2/svd_5873.json"), encoding="utf-8"))["data"]
samp = json.load(open(os.path.join(ROOT, "data/raw/gemma_5873_samples.json"), encoding="utf-8"))
if isinstance(samp, dict):
    samp = samp.get("data") or []

# 逐样本特征
lab = {}
for e in samp:
    s = e.get("sample") or e
    bid = s.get("id") or (s.get("bioAssay") or {}).get("id")
    m = {}
    for c in (s.get("characteristics") or []):
        m.setdefault((c.get("category") or "").lower(), set()).add(c.get("value"))
    lab[bid] = m

ids = svd["bioAssayIds"]
V = np.array(svd["vmatrix"], dtype=float)
var = np.array(svd["variances"], dtype=float)
print("SVD: %d 样本 × %d 主成分；PC1 方差占比 %.3f，PC1-3 累计 %.3f" % (V.shape[0], V.shape[1], var[0], var[:3].sum()))

def col(keys):
    out = []
    for bid in ids:
        m = lab.get(bid, {})
        v = None
        for k in keys:
            if k in m:
                v = sorted(x for x in m[k] if x)
                break
        out.append(";".join(v) if v else None)
    return np.array(out, dtype=object)

dis = col(["disease", "disease state", "diagnosis"])
bat = col(["block", "batch", "scan date", "label", "scan_date"])
print("\n疾病标签:", {v: int((dis == v).sum()) for v in set(dis) if v})
print("批次标签:", {v: int((bat == v).sum()) for v in set(bat) if v})

def auc(score, y):
    y = np.asarray(y, dtype=float)
    if y.sum() == 0 or y.sum() == len(y):
        return np.nan
    order = np.argsort(score)
    ranks = np.empty(len(score), dtype=float)
    ranks[order] = np.arange(1, len(score) + 1)
    n1, n0 = y.sum(), len(y) - y.sum()
    return (ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

uc = np.array([("ulcerative colitis" in str(x)) or (str(x).lower().startswith("uc")) for x in dis])
ctrl = np.array([("control" in str(x).lower()) or ("reference" in str(x).lower()) or ("normal" in str(x).lower()) for x in dis])
keep = uc | ctrl
y = uc[keep].astype(int)
print("\n可比样本：%d（UC %d / 对照 %d）" % (keep.sum(), uc.sum(), ctrl.sum()))
if keep.sum() > 10 and 0 < y.sum() < len(y):
    for pc in range(min(3, V.shape[1])):
        a = auc(V[keep, pc], y)
        print("  PC%d 单轴预测疾病 AUC = %.3f（方向无关，取 max(a,1-a) = %.3f）" % (pc + 1, a, max(a, 1 - a)))
    # 批次可分性（对照我们算过的 0.77）
    bl = bat[keep]
    for b in set(bl):
        if b and (bl == b).sum() >= 5:
            yb = (bl == b).astype(int)
            print("  批次「%s」在 PC1 上的 AUC = %.3f" % (b, max(auc(V[keep, 0], yb), 1 - auc(V[keep, 0], yb))))
