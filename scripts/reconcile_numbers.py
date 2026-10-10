# -*- coding: utf-8 -*-
"""数字—来源逐项核对（§C5 规则的可复用实现，v2 严格版）。

规则：英文稿正文（Title–§6）出现的任何数字 token，必须能在中文附注 §B 中找到
（允许 0.5% / 0.02 以内的舍入），或在白名单（章节号、年份、平台 ID、基因组版本号）内。

v2 相比 v1 的改动：
  * 剥离引用标记 [n] / [n,m] / [n–m]，避免把文献编号当数据；
  * 千分位归一化（12,366 -> 12366），两侧都做；
  * 剥离基因组版本 GRCh37/38/39、hg19/38；
  * 舍入容差改为相对 0.5%（下限 0.02），可容纳 57.35 -> 57.4、78.3 -> 78。
"""
import re

EN = "C:/Users/H-P/WorkBuddy/应激颗粒/results/主论文初稿_v1_EN.md"
CN = "C:/Users/H-P/WorkBuddy/应激颗粒/results/主论文初稿_v1_CN附注.md"

en = open(EN, encoding="utf-8").read()
cn = open(CN, encoding="utf-8").read()


def normalise(text, strip_citations=False):
    if strip_citations:
        # citation markers: [12] [12,13] [12–17] [3,4]
        text = re.sub(r"\[[\d,\s\u2010-\u2015]+\]", " ", text)
    text = re.sub(r"10\.\d{4,5}/[^\s)\]}]+", " ", text)  # DOIs (identifiers, not results)
    text = re.sub(r"GRCh3[789]|hg1[89]", " ", text)
    text = re.sub(r"\d{1,2}:\d{2}", " ", text)          # clock times (23:12, 20:38)
    text = re.sub(r"(?<=\d),(?=\d{3}\b)", "", text)   # 12,366 -> 12366
    return text


body = normalise(en.split("## 7. References")[0], strip_citations=True)

cn_b = cn[cn.find("## §B"):cn.find("## §C")]
cn_b = normalise(cn_b)

NUM = re.compile(r"\d+(?:\.\d+)?")
body_nums = set(NUM.findall(body))
cn_vals = [float(x) for x in NUM.findall(cn_b)]

WHITELIST = {
    # 章节号
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
    "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31",
    "2.0", "2.1", "2.2", "2.3", "2.4", "2.5", "2.6", "2.7", "2.8", "6.4",
    # 年份 & 工具版本 & 阈值
    "2002", "2005", "2007", "2013", "2016", "2017", "2018", "2019", "2020", "2021", "2022",
    "2023", "2024", "2025", "2026", "3.13", "4.6", "5.2", "0.001", "0.01", "0.05", "0.1",
    # 平台 ID
    "004133", "0010494", "0035061", "5873", "8385", "8122", "22619", "19738", "53987", "80655",
    "133060", "38713", "36376", "7673", "157783", "140231", "204796", "7621", "20164", "20163",
    "8397", "49036", "65372", "25794", "1746", "22952",
    # 1.4e-5 是 §B3「p = 1.44e-5」的两位有效数字写法
    "1.4",
}


def ok(n):
    if n in WHITELIST:
        return True
    v = float(n)
    tol = max(0.02, 0.005 * abs(v))
    return any(abs(v - c) <= tol for c in cn_vals)


missing = {n: None for n in body_nums if not ok(n)}
for n in list(missing):
    pat = re.compile(r"(?<![\d.])" + re.escape(n) + r"(?![\d.])")
    m = pat.search(body)
    if m:
        missing[n] = body[max(0, m.start() - 95):m.start() + 45].replace("\n", " ")

print("body numeric tokens : %d" % len(body_nums))
print("CN §B numeric values: %d" % len(cn_vals))
print()
if missing:
    print("NOT FOUND in CN appendix §B (%d):" % len(missing))
    for n in sorted(missing, key=float):
        print("  %-9s ...%s..." % (n, missing[n]))
else:
    print("OK — every number in the English body is in the CN appendix §B or whitelisted.")
