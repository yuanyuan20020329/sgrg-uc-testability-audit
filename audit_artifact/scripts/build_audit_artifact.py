# -*- coding: utf-8 -*-
"""把第十一轮审稿所依赖的 GEO 审计中间结果整理成一个自洽、可哈希的工件包。

为什么要做
----------
审稿人指出：本稿的核心主张是 provenance 与可复现，而支撑「截断遗漏了
GSE107593」这一结论的中间文件只存在于作者工作目录中。本脚本把这些文件
从 `out/` 汇总到 `resource_paper/audit_artifact/`，加上 README、逐文件
SHA-256 清单与一个 zip，使其成为一个可以直接上传到 Zenodo 的独立记录。

★ 设计原则
- 不复制 `deposit/` 的任何字节：工件包是**追加**证据，不改动已发布数据包。
- 每个文件都在 README 中写明来源、生成方与「能核验什么／不能核验什么」。
- 凡本版未保留的中间量（2026-10-15 的 UID 清单），在 README 中显式列为缺口，
  不假装完整。
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import partD_resfeas_geo as G  # noqa: E402

AUD = ROOT / "out" / "geo_retrieval_audit"
VOC = ROOT / "out" / "vocab_audit"
QSE = ROOT / "out" / "query_sensitivity"
DST = ROOT / "resource_paper" / "audit_artifact"
VER = "1.0.0"
RUN = "2026-10-17"

# 样本级 SOFT 只读 scratch，绝不写回归档缓存
G.CACHE = str(AUD / "soft_cache_scratch")


def rmrf(p: Path) -> None:
    if p.is_dir():
        shutil.rmtree(p)


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    rmrf(DST)
    for sub in ("geo", "rules", "queries", "scripts"):
        (DST / sub).mkdir(parents=True, exist_ok=True)

    # ---------- GEO 检索 ----------
    arch = pd.read_csv(ROOT / "results" / "partD_resfeas_geo_series_raw.tsv", sep="\t")
    arch = arch[["uid", "gse", "title", "n_samples", "pdat"]] if "summary" in arch.columns else arch
    assert len(arch) == 521, len(arch)
    arch.to_csv(DST / "geo" / "uids_archived_20261009.tsv", sep="\t", index=False)

    c600 = pd.read_csv(AUD / "gse_cap600_human.tsv", sep="\t")
    full = pd.read_csv(AUD / "gse_full_human.tsv", sep="\t")
    c600.to_csv(DST / "geo" / "uids_retmax600_20261017.tsv", sep="\t", index=False)
    full.to_csv(DST / "geo" / "uids_retmax5000_20261017.tsv", sep="\t", index=False)

    u600 = [x.strip() for x in (AUD / "uids_cap600_union.txt").read_text().split() if x.strip()]
    ufull = [x.strip() for x in (AUD / "uids_full_union.txt").read_text().split() if x.strip()]
    (DST / "geo" / "uids_union_retmax600_20261017.txt").write_text(
        "\n".join(u600) + "\n", encoding="utf-8")
    (DST / "geo" / "uids_union_retmax5000_20261017.txt").write_text(
        "\n".join(ufull) + "\n", encoding="utf-8")

    # 双向差集（成员层面，逐条）
    set_arch, set600, setfull = set(arch.gse), set(c600.gse), set(full.gse)
    d1 = sorted(set_arch - set600)
    d2 = sorted(set600 - set_arch)
    d3 = sorted(setfull - set600)
    pd.DataFrame({"gse": d1}).to_csv(
        DST / "geo" / "diff_archived_minus_600_20261017.tsv", sep="\t", index=False)
    c600[c600.gse.isin(d2)].to_csv(
        DST / "geo" / "diff_600_minus_archived_20261017.tsv", sep="\t", index=False)

    trunc = pd.read_csv(AUD / "gse_truncated_by_retmax600.tsv", sep="\t")
    assert len(trunc) == 60 == len(d3), (len(trunc), len(d3))
    trunc.to_csv(DST / "geo" / "truncated_by_retmax600_20261017.tsv", sep="\t", index=False)

    # 逐式命中与取回
    summ = json.loads((AUD / "summary.json").read_text(encoding="utf-8"))
    qc = summ["query_counts"]["full"]
    rows = []
    for q in ("hts", "spatial", "scnuc", "array"):
        n = qc[q]
        rows.append({"query": q, "n_hits": n,
                     "n_retrieved_retmax600": min(n, 600),
                     "n_retrieved_retmax5000": n,
                     "hit_the_cap": n > 600})
    pd.DataFrame(rows).to_csv(DST / "geo" / "per_query_hits_20261017.tsv", sep="\t", index=False)

    # 截断组 vs 已纳入组的构成对照
    inc, tru = summ["included_cap600"], summ["truncated_full_minus_cap600"]
    comp = [
        ("系列数", inc["n_series"], tru["n_series"]),
        ("样本数中位数", inc["n_samples"]["median"], tru["n_samples"]["median"]),
        ("样本数 Q1", inc["n_samples"]["q1"], tru["n_samples"]["q1"]),
        ("样本数 Q3", inc["n_samples"]["q3"], tru["n_samples"]["q3"]),
        ("样本数范围下限", inc["n_samples"]["min"], tru["n_samples"]["min"]),
        ("样本数范围上限", inc["n_samples"]["max"], tru["n_samples"]["max"]),
        ("样本数 > 400 的系列", inc["n_samples_gt400"], tru["n_samples_gt400"]),
        ("年份中位数", inc["year"]["median"], tru["year"]["median"]),
        ("年份范围下限", inc["year"]["min"], tru["year"]["min"]),
        ("年份范围上限", inc["year"]["max"], tru["year"]["max"]),
        ("标题/摘要含 UC", inc["label_UC"], tru["label_UC"]),
        ("标题/摘要含 CD", inc["label_CD"], tru["label_CD"]),
        ("标题/摘要含 IBD", inc["label_IBD"], tru["label_IBD"]),
        ("含任一 IBD 标签", inc["label_UC_or_CD_or_IBD"], tru["label_UC_or_CD_or_IBD"]),
        ("组织=结肠", inc["site_colon"], tru["site_colon"]),
        ("组织=回肠", inc["site_ileum"], tru["site_ileum"]),
        ("组织=直肠", inc["site_rectum"], tru["site_rectum"]),
        ("组织未识别", inc["site_none"], tru["site_none"]),
        ("模态代理=mRNA", inc["molecule"]["mRNA"], tru["molecule"]["mRNA"]),
        ("模态代理=chromatin", inc["molecule"]["chromatin"], tru["molecule"]["chromatin"]),
        ("模态代理=ncRNA", inc["molecule"]["ncRNA"], tru["molecule"]["ncRNA"]),
        ("模态代理=other", inc["molecule"]["other"], tru["molecule"]["other"]),
        ("像单细胞或空间", inc["is_sc_or_spatial"], tru["is_sc_or_spatial"]),
        ("含离体模型线索", inc["in_vitro_hint"], tru["in_vitro_hint"]),
    ]
    pd.DataFrame(comp, columns=["indicator", "included_retmax600", "truncated_full_minus_600"]).to_csv(
        DST / "geo" / "truncation_comparison.tsv", sep="\t", index=False)

    # SOFT 同源规则复核
    soft = pd.read_csv(AUD / "truncated_soft_audit.tsv", sep="\t")
    soft.to_csv(DST / "geo" / "truncated_soft_recheck.tsv", sep="\t", index=False)

    # ---------- GSE107593 逐样本证据 ----------
    gse = "GSE107593"
    samples = G.parse_soft(gse, 400)
    assert isinstance(samples, list) and len(samples) == 48, type(samples)

    recs = []
    donor_states: dict[str, list[str]] = {}
    pair_site: dict[tuple, list[str]] = {}
    for s in samples:
        txt = " ".join([s["title"], s["source"], s["chars"]])
        if G.RE_NONINFL.search(txt):
            st = "noninflamed"
        elif G.RE_INFLAM.search(txt):
            st = "inflamed"
        elif G.RE_HEALTHY.search(txt):
            st = "healthy"
        else:
            st = "unknown"
        dn = G.donor_token(s["title"], s["source"])
        site = G.site_of(txt)
        donor_states.setdefault(dn or "<%s#%s>" % (st, s["gsm"]), []).append(st)
        pair_site.setdefault((dn or s["gsm"], site or "na"), []).append(st)
        recs.append({"gsm": s["gsm"], "sample_title": s["title"],
                     "source_name_ch1": s["source"], "donor_token": dn or "",
                     "site": site or "", "state": st,
                     "disease_token_UC": bool(G.RE_UC.search(txt))})
    per = pd.DataFrame(recs)

    donor_df = pd.DataFrame(
        [{"donor_token": k,
          "n_samples": len(v),
          "states": ";".join(sorted(set(v))),
          "is_paired_donor": {"inflamed", "noninflamed"} <= set(v)}
         for k, v in sorted(donor_states.items())])
    site_df = pd.DataFrame(
        [{"donor_token": k[0], "site": k[1], "n_samples": len(v),
          "states": ";".join(sorted(set(v))),
          "is_paired_same_site": k[1] != "na" and {"inflamed", "noninflamed"} <= set(v)}
         for k, v in sorted(pair_site.items())])

    per.to_csv(DST / "geo" / "gse107593_sample_level.tsv", sep="\t", index=False)
    donor_df.to_csv(DST / "geo" / "gse107593_donor_states.tsv", sep="\t", index=False)
    site_df.to_csv(DST / "geo" / "gse107593_donor_site_pairs.tsv", sep="\t", index=False)

    n_paired = int(donor_df.is_paired_donor.sum())
    n_paired_site = int(site_df.is_paired_same_site.sum())

    # 四道闸门逐项判定（复用发布流水线的判定函数，不手写）
    row = trunc[trunc.gse == gse].iloc[0]
    ser = G.parse_soft_series(gse)
    blob = str(row.title) + " " + (str(row.get("summary")) if pd.notna(row.get("summary")) else "")
    mol, res = G.assay_class(ser.get("series_type", ""), row.gdsType, blob)
    invitro = bool(G.RE_INVITRO.search(blob))
    gates = pd.DataFrame([
        {"gate": "G1_mRNA", "value": mol == "mRNA",
         "evidence": f"!Series_type -> molecule={mol}"},
        {"gate": "G2_UC", "value": bool(row.lbl_UC),
         "evidence": "系列疾病标签 token 含 UC"},
        {"gate": "G3_powered_control", "value": n_paired >= 6,
         "evidence": f"闸门输入为 n_paired_donors（同供体，不要求同部位）={n_paired}，阈值 >= 6"},
        {"gate": "G4_in_vivo", "value": not invitro,
         "evidence": "未命中离体模型正则" if not invitro else "命中离体模型正则"},
        {"gate": "（敏感性）若改用 n_paired_same_site", "value": n_paired_site >= 6,
         "evidence": f"同供体且同部位 = {n_paired_site}；闸门若改用该口径则 G3 为假"},
    ])
    gates["gse"] = gse
    gates.to_csv(DST / "geo" / "gse107593_gate_evaluation.tsv", sep="\t", index=False)

    # ---------- 词表 ----------
    (DST / "rules").mkdir(exist_ok=True)
    shutil.copyfile(VOC / "vocab_audit.json", DST / "rules" / "vocab_audit.json")
    vj = json.loads((VOC / "vocab_audit.json").read_text(encoding="utf-8"))
    pd.DataFrame([{"branch": k, **{kk: vv for kk, vv in v.items() if kk != "examples"},
                   "examples": ";".join(v.get("examples", []))}
                  for k, v in vj["per_branch"].items()]).to_csv(
        DST / "rules" / "vocab_branch_hits.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"sensitivity_variant": "行首锚定 -> 子串搜索",
         "n_rows_arm_class_changed": vj["sensitivity"]["prefix_anchor_vs_substring"]["n_rows_arm_class_changed"]},
        {"sensitivity_variant": "移除 normal 分支",
         "n_rows_arm_class_changed": vj["sensitivity"]["drop_normal_branch"]["n_rows_arm_class_changed"]},
        {"sensitivity_variant": "整表移除（全部记为疾病标签）",
         "n_rows_arm_class_changed": vj["sensitivity"]["drop_vocabulary_entirely"]["n_rows_arm_class_changed"]},
    ]).to_csv(DST / "rules" / "arm_class_sensitivity.tsv", sep="\t", index=False)

    # ---------- 检索式敏感性 ----------
    for f in ("query_sensitivity_summary.json", "query_sensitivity_triaged_uids.tsv"):
        shutil.copyfile(QSE / f, DST / "queries" / f)

    # ---------- 脚本 ----------
    for f in ("geo_retrieval_audit.py", "geo_truncation_soft.py",
              "vocab_audit.py", "geo_query_sensitivity.py"):
        shutil.copyfile(ROOT / "out" / f, DST / "scripts" / f)
    shutil.copyfile(ROOT / "out" / "build_audit_artifact.py", DST / "scripts" / "build_audit_artifact.py")

    # ---------- README ----------
    readme = f"""# sgrg GEO audit artifact — v{VER}

本工件包汇总《面向溃疡性结肠炎及相关炎性肠病公共转录组元数据的来源可溯性与预设筛选资格标注资源》
（下称"主稿"）在第十轮与第十一轮审稿应答中生成的全部 **GEO 检索审计中间结果**。
它是对已发布数据包 **1.0.2 版（10.5281/zenodo.23414533）** 的**追加证据**，不修改、不替换其中任何字节。

- 主稿版本：v3 → v4（2026-10-17 修订）
- 生成日期：{RUN}
- 许可：CC BY 4.0
- 逐文件校验：`MANIFEST.sha256`（SHA-256，按相对路径）

## 1 这些文件回答什么问题

主稿报告的以下几类结论，**不能**从 1.0.2 数据包离线重建，只能由重跑检索得到：

1. 四条检索式的逐式命中数与并集 UID 数；
2. 521 / 522 / 523 / 582 / 583 五个计数各自的来源与成员级双向差集；
3. `retmax = 600` 截断丢弃的 60 条系列清单及其与已纳入组的构成对照；
4. 这 60 条按发布表同源 `!Series_type` 规则复核的分子类别与分辨率；
5. `GSE107593` 满足全部四道闸门的逐样本依据；
6. 对照词表的分支命中与敏感性；
7. 检索式词表变体的增量与分诊结果。

本包即为上述七项的落盘产物。

## 2 文件清单

### 2.1 `geo/` — 检索与截断

| 文件 | 内容 | 行数 |
|---|---|---:|
| `uids_archived_20261009.tsv` | 归档流程 2026-10-09 运行**自身产出**的 521 条人类 GSE（uid, gse, title, n_samples, pdat） | 521 |
| `uids_retmax600_20261017.tsv` | 2026-10-17 复算，`retmax = 600` 口径的人类 GSE 全字段表 | 523 |
| `uids_retmax5000_20261017.tsv` | 2026-10-17 复算，`retmax = 5000` 全口径 | 583 |
| `uids_union_retmax600_20261017.txt` | 上述 600 口径下四式并集去重后的 UID（每行一个） | 1,269 |
| `uids_union_retmax5000_20261017.txt` | 全口径四式并集去重后的 UID | 2,734 |
| `diff_archived_minus_600_20261017.tsv` | 归档 521 − 复算 523（**仅表头，0 行**） | 0 |
| `diff_600_minus_archived_20261017.tsv` | 复算 523 − 归档 521（`GSE339323`、`GSE350397`） | 2 |
| `per_query_hits_20261017.tsv` | 四条检索式的命中数与取回数 | 4 |
| `truncated_by_retmax600_20261017.tsv` | 全口径多出的 60 条人类 GSE | 60 |
| `truncation_comparison.tsv` | 已纳入组（523）与被截断组（60）的 24 项构成指标 | 24 |
| `truncated_soft_recheck.tsv` | 60 条的 `!Series_type` 复核：分子类别、分辨率、候选配对计数 | 60 |
| `gse107593_sample_level.tsv` | `GSE107593` 48 个样本逐条解析结果（供体 token、部位、状态） | 48 |
| `gse107593_donor_states.tsv` | 该系列按供体汇总的状态集合与是否构成配对 | {len(donor_df)} |
| `gse107593_donor_site_pairs.tsv` | 该系列按（供体, 部位）汇总的配对判定 | {len(site_df)} |
| `gse107593_gate_evaluation.tsv` | 该系列四道闸门的逐项取值与依据 | 4 |

**成员级双向差集只在"归档 521"与"复算 523"之间完整成立**：两份清单都在本包内，`archived − 600` 为 0 行、
`600 − archived` 为 2 行，可逐条核对。

**`GSE107593` 的关键细节：两个配对口径给出不同结果。** 该系列按供体汇总的候选配对为
**{n_paired} 对**，而按"同供体且同部位"计只有 **{n_paired_site} 对**——差别来自供体 157、1057、1192
的炎症与非炎症样本取自不同部位（见 `gse107593_sample_level.tsv` 的 `site` 列，其值取自
`!Sample_characteristics_ch1` 的 `location:` 字段）。**台账 `G3` 的输入是供体级计数 `n_paired_donors`**
（`scripts/partD_resfeas_ledger.py`：`g3 = n_pr >= PAIR_MIN`，`PAIR_MIN = 6`），故该系列以 6 对恰好达到阈值；
若闸门改用同部位口径（3 对），它**不会**通过。因此"取全口径后会成为台账第三个通过行"这一反事实陈述，
成立的前提是**沿用已发布的闸门输入定义**，而不是"该系列有 6 对同部位配对"。

### 2.2 `rules/` — 对照词表

| 文件 | 内容 |
|---|---|
| `vocab_audit.json` | 332 个唯一标签的词表命中、误分类、近失命中与三项敏感性统计 |
| `vocab_branch_hits.tsv` | 词表 12 个分支各自的命中标签数与数据集数 |
| `arm_class_sensitivity.tsv` | 三种词表改动下 `arm_class` 发生变化的行数 |

### 2.3 `queries/` — 检索式敏感性

| 文件 | 内容 |
|---|---|
| `query_sensitivity_summary.json` | 各词表变体的命中数、与基线的重叠与增量；仅对增量 ≤ 阈值者做了分诊 |
| `query_sensitivity_triaged_uids.tsv` | 已分诊变体的增量 UID 清单 |

### 2.4 `scripts/` — 生成上述文件的脚本

`geo_retrieval_audit.py`（检索复算与截断对照）、`geo_truncation_soft.py`（60 条 SOFT 复核）、
`vocab_audit.py`（词表）、`geo_query_sensitivity.py`（检索式敏感性）、`build_audit_artifact.py`（本包打包）。
脚本依赖 `scripts/partD_resfeas_geo.py`（与发布表同一份检索/解析代码），该文件属**代码归档**
（10.5281/zenodo.23412772）而非本包。

## 3 本包**不能**核验什么（重要）

1. **2026-10-15 的 UID 清单未保留。** 该时点的 522（`retmax = 600`）与 582（全口径）**只有计数**，
   其 UID 清单在当次运行后未落盘，本包不含。因此 521 − 522 = 0 与 522 − 521 = 1（`GSE350397`）
   这两个成员级断言**无法由本包复核**，只能由计数与首次公开日期归因得到。主稿已就此改写表述。
2. **源库当时的状态不可回溯。** 上述所有清单都是特定时点的快照；重跑会得到不同结果
   （2026-10-15 → 2026-10-17 两天内即漂移 1 条）。
3. **供体配对仍是字符串规则输出。** `gse107593_*` 三张表给的是**规则如何得出这 6 对**，
   不是人工金标准；若规则本身在别处产生假阳性（主稿表 5 的两个实例），本表不能排除同类问题。
4. **模态判定是规则性的。** `truncated_soft_recheck.tsv` 的 `molecule_soft` / `resolution_soft`
   由 `!Series_type` 与平台推出，非人工逐条判读。
5. **不覆盖 CELLxGENE 侧。** 该侧可由 1.0.2 数据包内的冻结快照离线重建，无需本包。

## 4 复现方式

```bash
# 1) 复算检索与截断对照（需可访问 NCBI E-utilities）
python scripts/geo_retrieval_audit.py
# 2) 对 60 条被截断系列做同源 SOFT 复核
python scripts/geo_truncation_soft.py
# 3) 词表覆盖与敏感性（只读 1.0.2 数据包的四个 TSV）
python scripts/vocab_audit.py
# 4) 检索式词表敏感性（需 NCBI）
python scripts/geo_query_sensitivity.py
```

脚本以主稿仓库的 `scripts/partD_resfeas_geo.py` 为解析核心，故须在仓库根目录下运行。
`geo_truncation_soft.py` 把样本级 SOFT 缓存重定向到 `out/geo_retrieval_audit/soft_cache_scratch`，
不写回归档缓存。

## 5 与已发布数据包的关系

- 数据包 1.0.2（10.5281/zenodo.23414533）的 14 个文件与 1.0.1 **逐字节相同**；本包的加入不改变这一点。
- 本包若被合并进数据包，按主稿的语义化版本规则属 **minor** 升级（新增文件、不改既有字段语义）。
"""

    (DST / "README.md").write_text(readme, encoding="utf-8")

    # ---------- MANIFEST ----------
    files = sorted(p for p in DST.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256")
    lines = [f"{sha256(p)}  {p.relative_to(DST).as_posix()}" for p in files]
    (DST / "MANIFEST.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # ---------- zip ----------
    zpath = ROOT / "resource_paper" / f"sgrg-geo-audit-artifact_v{VER}.zip"
    if zpath.exists():
        zpath.unlink()
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for p in sorted(DST.rglob("*")):
            if p.is_file():
                z.write(p, f"sgrg-geo-audit-artifact_v{VER}/{p.relative_to(DST).as_posix()}")

    print("== 打包完成 ==")
    print("目录:", DST)
    print("文件数:", len(files), "| 总字节:", sum(p.stat().st_size for p in files))
    print("zip:", zpath.name, zpath.stat().st_size, "B")
    print(f"GSE107593: 配对供体 {n_paired} / 同部位配对 {n_paired_site} / molecule={mol}")
    for g in gates.itertuples():
        print(f"   {g.gate:20s} = {g.value}   ({g.evidence})")


if __name__ == "__main__":
    main()
