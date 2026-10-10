# -*- coding: utf-8 -*-
"""案例 B：目标论文（Sci Rep 2024, PD pyroptosis, PMC11574289）细节抽取。

抽取决定审计可行性的要素：
  · 各 GSE 的角色（发现 / 验证 / 单细胞）与样本量
  · 对照组定义（健康对照来自哪里）
  · 基因集来源与面板大小
  · 报告的核心数字（AUC / n 基因）
  · 数据可得性
"""
import re
import sys

path = sys.argv[1] if len(sys.argv) > 1 else "data/raw/caseB/paper/PMC11574289_fullText.xml"
x = open(path, encoding="utf-8").read()


def clean(s):
    s = re.sub(r"<[^>]+>", " ", s)
    return re.sub(r"\s+", " ", s).strip()


txt = clean(x)


def ctx(kw, n=3, w=420, flags=re.I):
    out = []
    for m in list(re.finditer(re.escape(kw), txt, flags))[:n]:
        i = m.start()
        out.append(txt[max(0, i - w): i + w])
    return out


print("### 数据集与角色 ###")
for kw in ["GSE140231", "GSE157783", "GSE20163", "GSE204796", "GSE49036",
           "GSE7621", "GSE8397", "GSE20164"]:
    c = ctx(kw, n=1, w=300)
    if c:
        print("-- %s --" % kw)
        print("   ", c[0][:340])
print()

print("### 基因集 / 面板 ###")
for kw in ["pyroptosis-related genes", "PRGs", "GSEA", "consensus", "LASSO",
           "SVM", "random forest", "hub gene", "marker gene"]:
    for s in ctx(kw, n=2, w=320):
        print("[%s] %s" % (kw, s[:340]))
    print()

print("### 核心数字 ###")
for kw in ["AUC", "sensitivity", "specificity", "GSEA", "immune"]:
    for s in ctx(kw, n=3, w=260):
        print("[%s] %s" % (kw, s[:300]))
    print()

print("### 数据可得性 ###")
for s in ctx("Data availability", n=1, w=600):
    print(s)
