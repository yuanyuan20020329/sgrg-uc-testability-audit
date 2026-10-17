# sgrg GEO audit artifact — v1.0.0

本工件包汇总《面向溃疡性结肠炎及相关炎性肠病公共转录组元数据的来源可溯性与预设筛选资格标注资源》
（下称"主稿"）在后续各轮审稿应答中生成的全部 **GEO 检索审计中间结果**。
它是对已发布数据包 **1.0.2 版（10.5281/zenodo.23414533）** 的**追加证据**，不修改、不替换其中任何字节。

- 主稿版本：v3（2026-10-17 修订，对应 docx `中文稿草稿_v12`）
- 生成日期：2026-10-17
- 公开位置：代码仓库 `yuanyuan20020329/sgrg-uc-testability-audit` 的标签 **`audit-artifact-v1.0.0`**
  下的 `audit_artifact/` 目录（`https://github.com/yuanyuan20020329/sgrg-uc-testability-audit/tree/audit-artifact-v1.0.0/audit_artifact`）
- **独立历史根**：该标签指向 `audit-artifact` 分支上的一个**无父提交**，其树中只有 `audit_artifact/`。
  因此克隆该标签**只得到本包**，不含代码归档的任何文件；本包在取得方式上等同于一个独立仓库。
- 归档下载：`https://github.com/yuanyuan20020329/sgrg-uc-testability-audit/archive/refs/tags/audit-artifact-v1.0.0.zip`
- 许可：CC BY 4.0
- 逐文件校验：`MANIFEST.sha256`（SHA-256，按相对路径；共 26 行，不含清单自身）

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
| `gse107593_donor_states.tsv` | 该系列按供体汇总的状态集合与是否构成配对 | 30 |
| `gse107593_donor_site_pairs.tsv` | 该系列按（供体, 部位）汇总的配对判定 | 34 |
| `gse107593_gate_evaluation.tsv` | 该系列四道闸门的逐项取值与依据 | 4 |

**成员级双向差集只在"归档 521"与"复算 523"之间完整成立**：两份清单都在本包内，`archived − 600` 为 0 行、
`600 − archived` 为 2 行，可逐条核对。

**`GSE107593` 的关键细节：两个配对口径给出不同结果。** 该系列按供体汇总的候选配对为
**6 对**，而按"同供体且同部位"计只有 **3 对**——差别来自供体 157、1057、1192
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
- **本包不是 1.0.2 数据包的一部分**，也不在 1.0.2 的数据 DOI 覆盖范围内；它单独公开于代码仓库的固定标签下，
  因此引用时应同时给出标签名与文件相对路径，并核对 `MANIFEST.sha256`。
- 1.0.2 数据包 `datapackage.json` 中 `G3_powered_control` 的字段描述仍写作 *design-eligibility gate*，
  与主稿采用的“未验证的设计线索闸门”不一致；该措辞将在下一数据版本（v1.0.3）统一修正，
  属**语义修正、不改任何取值**（见主稿“局限与已知问题”第 2 条）。

## 6 如何取得与核对本包

```bash
# 取得（只含本包，不含代码归档）
git clone --branch audit-artifact-v1.0.0 --depth 1 \
  https://github.com/yuanyuan20020329/sgrg-uc-testability-audit.git sgrg-audit
cd sgrg-audit/audit_artifact

# 或直接下载 zip
# https://github.com/yuanyuan20020329/sgrg-uc-testability-audit/archive/refs/tags/audit-artifact-v1.0.0.zip

# 核对（在 audit_artifact/ 目录下）
sha256sum -c MANIFEST.sha256     # Linux
shasum -a 256 -c MANIFEST.sha256 # macOS
```

26 行全部 `OK` 即与主稿引用的字节一致。若需核对主稿结论，建议顺序：
先看 `geo/diff_*.tsv` 与 `geo/uids_*.tsv`（成员关系），再看 `geo/truncation_comparison.tsv`
（截断组构成），最后看 `geo/gse107593_*.tsv`（单条反事实实例）。
