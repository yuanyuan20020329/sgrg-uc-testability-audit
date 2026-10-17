# -*- coding: utf-8 -*-
"""对照词表覆盖与敏感性核验（对应审稿人 R3-1、R3-5）。

核验内容
--------
1) 唯一疾病标签总数、被词表命中的标签与数据集数、逐分支命中数；
2) 误分类清单（`^normal` 命中 `normal pressure hydrocephalus`）；
3) 词表微改的敏感性：① 行首锚定 vs 子串搜索；② 移除某一分支；③ 整表移除（= 全部记为疾病）；
   三种口径下 `arm_class` 发生变化的行数；
4) 「近失命中」——标签文本含对照关键词、但未被行首锚定命中的标签。

★ 字段结构：`diseases` 是用 `;` 分隔的标签列表；单个标签内部可能再用 `||` 连接**共病组合**
  （如 `dementia || Alzheimer disease`）。本脚本按 `;` 切分，与发布表、审计脚本口径一致。

输出：out/vocab_audit/vocab_audit.json
"""
from __future__ import annotations

import collections
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "out" / "vocab_audit"
OUT.mkdir(parents=True, exist_ok=True)

# 与正文完全相同的词表正则
BASE = r"^(normal|healthy|control|ctrl|non[- ]?disease|unaffected|wild[- ]?type|wt|na[iï]ve|saline|sham|baseline)\b"
PAT = re.compile(BASE, re.I)
KEY = re.compile(r"normal|healthy|control|ctrl|unaffected|non[\s\-]?disease|wild[\s\-]?type|"
                 r"na[iï]ve|saline|sham|baseline", re.I)
BRANCHES = {
    "normal": r"^normal\b", "healthy": r"^healthy\b", "control": r"^control\b",
    "ctrl": r"^ctrl\b", "non-disease": r"^non[- ]?disease\b", "unaffected": r"^unaffected\b",
    "wild-type": r"^wild[- ]?type\b", "wt": r"^wt\b", "naive": r"^na[iï]ve\b",
    "saline": r"^saline\b", "sham": r"^sham\b", "baseline": r"^baseline\b",
}


def arm_class_of(labels: list[str], pat: re.Pattern | None, mode: str = "match") -> str:
    """给定标签列表与对照判定正则，返回 ctl_only / dis_only / mixed。

    mode="match"  行首锚定（与发布表一致）
    mode="search" 子串搜索（放松锚定，等于把词表改成子串规则）
    pat=None      整表移除，全部记为疾病标签
    """
    ctl = dis = 0
    for x in labels:
        if pat is not None and (pat.search(x) if mode == "search" else pat.match(x)):
            ctl += 1
        else:
            dis += 1
    if ctl and dis:
        return "mixed"
    if ctl:
        return "ctl_only"
    return "dis_only"


def main() -> None:
    D = pd.read_csv(ROOT / "results" / "partD_resfeas_libwide_datasets.tsv", sep="\t")
    D["lab"] = D.diseases.fillna("").apply(lambda s: [x.strip() for x in s.split(";") if x.strip()])
    D["pub_arm"] = D.arm_class

    counter: collections.Counter = collections.Counter()
    for ls in D.lab:
        for x in ls:
            counter[x] += 1
    n_labels = len(counter)

    matched = {k: v for k, v in counter.items() if PAT.match(k)}
    near = {k: v for k, v in counter.items() if KEY.search(k) and not PAT.match(k)}

    # 逐分支
    per_branch = {}
    for b, rx in BRANCHES.items():
        p = re.compile(rx, re.I)
        ls = [k for k in counter if p.match(k)]
        per_branch[b] = {"n_labels": len(ls), "n_datasets": sum(counter[k] for k in ls),
                         "examples": ls[:6]}

    # 敏感性
    sub = re.compile(KEY.pattern, re.I)          # 子串搜索（放松行首锚定）
    D["arm_sub"] = D.lab.apply(lambda ls: arm_class_of(ls, sub, mode="search"))
    D["arm_none"] = D.lab.apply(lambda ls: arm_class_of(ls, None))
    no_normal = re.compile(BASE.replace("normal|", ""), re.I)
    D["arm_nonormal"] = D.lab.apply(lambda ls: arm_class_of(ls, no_normal))

    changed_sub = int((D.arm_sub != D.pub_arm).sum())
    changed_none = int((D.arm_none != D.pub_arm).sum())
    changed_nonormal = int((D.arm_nonormal != D.pub_arm).sum())

    # 放松锚定后新命中的标签（行首锚定下未被命中、子串搜索命中）
    relaxed_hits = sorted(set(k for k in counter if sub.search(k) and not PAT.match(k)))
    flip_rows = D.loc[D.arm_sub != D.pub_arm, ["dataset_id", "diseases", "pub_arm", "arm_sub"]] \
        if "dataset_id" in D.columns else D.loc[D.arm_sub != D.pub_arm, ["diseases", "pub_arm", "arm_sub"]]
    flip_rows.to_csv(OUT / "sensitivity_arm_class_flips.tsv", sep="\t", index=False)

    rep = {
        "n_datasets": int(len(D)),
        "n_unique_labels": n_labels,
        "n_labels_matched_by_vocabulary": len(matched),
        "n_datasets_matched_by_vocabulary": int(sum(
            1 for ls in D.lab if any(PAT.match(x) for x in ls))),
        "matched_labels": matched,
        "per_branch": per_branch,
        "misclassified_labels": {k: v for k, v in matched.items() if k != "normal"},
        "near_miss_labels": near,
        "n_near_miss_labels": len(near),
        "sensitivity": {
            "prefix_anchor_vs_substring": {
                "n_rows_arm_class_changed": changed_sub,
                "newly_matched_labels": relaxed_hits,
            },
            "drop_normal_branch": {"n_rows_arm_class_changed": changed_nonormal},
            "drop_vocabulary_entirely": {"n_rows_arm_class_changed": changed_none},
        },
        "arm_class_published": D.pub_arm.value_counts().to_dict(),
        "arm_class_under_substring_vocabulary": D.arm_sub.value_counts().to_dict(),
    }
    (OUT / "vocab_audit.json").write_text(
        json.dumps(rep, ensure_ascii=False, indent=2), encoding="utf-8")

    print("数据集 %d 行；唯一疾病标签 %d 个" % (len(D), n_labels))
    print("被词表命中的标签 %d 个；覆盖数据集 %d 个（%.1f%%）"
          % (len(matched), rep["n_datasets_matched_by_vocabulary"],
             100 * rep["n_datasets_matched_by_vocabulary"] / len(D)))
    print("\n逐分支命中：")
    for b, d in per_branch.items():
        print("   %-12s 标签 %2d 个 / 数据集 %5d 个" % (b, d["n_labels"], d["n_datasets"]))
    print("\n误分类清单:")
    for k, v in rep["misclassified_labels"].items():
        print("   %-45s %d 个数据集" % (k, v))
    print("\n近失命中（含关键词但未被行首锚定）: %d 个" % len(near))
    for k, v in sorted(near.items(), key=lambda x: -x[1]):
        print("   %-60s %d" % (k[:60], v))
    print("\n敏感性（arm_class 发生变化的行数）:")
    print("   行首锚定 → 子串搜索            : %d  （新命中标签：%s）" % (changed_sub, relaxed_hits))
    print("   移除 normal 分支               : %d" % changed_nonormal)
    print("   整表移除（全部记为疾病标签）    : %d" % changed_none)
    print("\n发布表 arm_class 分布:", rep["arm_class_published"])
    print("子串规则下 arm_class 分布:", rep["arm_class_under_substring_vocabulary"])


if __name__ == "__main__":
    main()
