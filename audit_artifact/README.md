# sgrg GEO audit artifact — v1.0.2

本工件包汇总《面向溃疡性结肠炎及相关炎性肠病公共分子数据集元数据的来源追溯与规则化筛选标注资源》
（下称"主稿"）在后续各轮审稿应答中生成的全部 **GEO 检索审计中间结果**。
它是对已发布数据包 **1.0.3 版（10.5281/zenodo.23436090）** 的**追加证据**，不修改、不替换其中任何字节。

- 主稿版本：v3（2026-10-17 修订，对应 docx `中文稿草稿_v17`）
- 生成日期：2026-10-17
- **本版（v1.0.2）相对 v1.0.1 的变更：只改 2 个文件，不增删文件**（仍 31 个文件）。
  `scripts/gdstype_recompute.py` 修正了一处输入缺陷——v1.0.1 的脚本写成
  `r.get("series_type", "")`，而归档 raw 表**没有** `series_type` 列（其第 6 列名为 `gdsType`），
  于是 `series_type` 也被传成空串，**同时改变了两个输入**，不再是"只改变 `gdstype` 一个输入"
  的对照实验。v1.0.2 改为从系列级 SOFT 缓存取 `!Series_type` 的真值（与归档流程同一函数），
  并把取到的值写进 `geo/gdstype_recompute_20261017.tsv` 新增的 `series_type` 列，供读者核对
  第三个输入确实非空。**重算结论不变**：423 行的 `molecule` 与 `resolution` 仍无一改变，
  521 行中 `changed` 为真的行数仍为 0。
- **本版新增的一条机理性解释**：521 行中 `!Series_type` 与 `gdstype` 的取值**逐行相同**
  （仅分隔符 `| ` 与 `; ` 不同，各 22 种取值；规范化后 521/521 相等）。这正是"`gdstype` 为空
  而结果不变"的原因——不是因为该输入无关紧要，而是 `!Series_type` 已携带完全相同的信息。
  脚本重跑时会打印该计数。
- **v1.0.1 相对 v1.0.0 的变更：新增 4 个文件**——`geo/skipped98_modality_recheck.tsv`
  （98 个未进入发布表的系列按发布表同一 `assay_class()` 规则重判的逐条结果）、
  `geo/gdstype_recompute_20261017.tsv`（按正确键名 `gdstype` 重取并重算 521 条的对照表）、
  以及两者的生成脚本 `scripts/skipped98_modality_recheck.py` 与 `scripts/gdstype_recompute.py`。
  **其余 27 个文件逐字节未动**，`audit-artifact-v1.0.0` 标签仍然有效、内容不变。
- 公开位置：代码仓库 `yuanyuan20020329/sgrg-uc-testability-audit` 的标签 **`audit-artifact-v1.0.2`**
  下的 `audit_artifact/` 目录（`https://github.com/yuanyuan20020329/sgrg-uc-testability-audit/tree/audit-artifact-v1.0.2/audit_artifact`）
- **独立历史根**：该标签指向 `audit-artifact` 分支上的一个**无父提交**，其树中只有 `audit_artifact/`。
  因此克隆该标签**只得到本包**，不含代码归档的任何文件；本包在取得方式上等同于一个独立仓库。
  v1.0.0、v1.0.1、v1.0.2 指向同一分支上的三个不同提交，前两个标签的树仍在、内容未被覆盖。
- 归档下载：`https://github.com/yuanyuan20020329/sgrg-uc-testability-audit/archive/refs/tags/audit-artifact-v1.0.2.zip`
- 许可：CC BY 4.0
- 逐文件校验：`MANIFEST.sha256`（SHA-256，按相对路径；共 30 行，不含清单自身）

## 1 这些文件回答什么问题

主稿报告的以下几类结论，**不能**从 1.0.3 数据包离线重建，只能由重跑检索得到：

1. 四条检索式的逐式命中数与并集 UID 数；
2. 521 / 522 / 523 / 582 / 583 五个计数各自的来源与成员级双向差集；
3. `retmax = 600` 截断丢弃的 60 条系列清单及其与已纳入组的构成对照；
4. 这 60 条按发布表同源规则复核的分子类别与分辨率；
5. `GSE107593` 满足全部四道闸门的逐样本依据；
6. 对照词表的分支命中与敏感性；
7. 检索式词表变体的增量与分诊结果；
8. 未进入发布表的 98 个系列的同源规则复核，以及 `gdstype` 输入未生效的影响评估（v1.0.1 新增，
   生成脚本的输入缺陷在 v1.0.2 修正）。

本包即为上述八项的落盘产物。

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
| `truncated_soft_recheck.tsv` | 60 条的复核：按发布表同一 `assay_class()` 给出的分子类别、分辨率与候选配对计数 | 60 |
| `skipped98_modality_recheck.tsv` | **v1.0.1 新增**：未进入发布表的 98 个系列的逐条模态复核（发布表规则与早期弱规则并置，附差异标记） | 98 |
| `gdstype_recompute_20261017.tsv` | **v1.0.1 新增／v1.0.2 修正**：521 条的对照——`series_type` 与 `gdstype` 两个输入的原值、发布值、重算值与是否改变 | 521 |
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

**98 个未进入发布表的系列：两个规则并置，不合并成一个结论（`skipped98_modality_recheck.tsv`）。**
该表逐条给出四列判定：`molecule_true_rule` / `resolution_true_rule` 为**发布表同一函数** `assay_class()`
在 `!Series_type`、`gdsType` 与"标题＋摘要前 400 字符"上的输出；`molecule_old_rule` / `resolution_old_rule`
为早期一次复核所用的**弱规则**（只拼 `!Series_type` 与标题、分辨率不区分阵列与测序），`disagree` 标记两者是否不同。
读法须注意三点：

1. **分子类别上两条规则完全一致**（mRNA 86、chromatin 7、ncRNA 5）；86 行标记为 `disagree`，**全部**来自
   分辨率粒度——弱规则只给 `bulk`，发布表规则再把它拆成 `bulk-array` 49 与 `bulk-seq` 37。因此 `disagree`
   一列**不是**"规则互相矛盾"的计数，不能读作 86 条判读不一致。
2. **无一条为单细胞或空间系列**（`resolution_true_rule` 的取值集合只有 `bulk-array`、`bulk-seq`、
   `chromatin`、`ncRNA`）；这是**该规则内**的结论，不是对源库的绝对陈述。
3. **归档 raw 表名为 `gdsType` 的一列在表内全为空**，原因见 §3 第 6 条：归档流程读的是 esummary
   的 `gdsType` 键，而接口实际返回 `gdstype`，故该项输入在归档中从未携带信息。本表复现的是
   **发布当时的实际输入**；用正确键名补取该值、并同时取回 `!Series_type` 的真值后重算，
   98 条的构成与 423 行的已发布取值均未改变（§3 第 6 条）。

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
`skipped98_modality_recheck.py`（**v1.0.1 新增**，98 条未解析系列的同源复核）、
`gdstype_recompute.py`（**v1.0.1 新增、v1.0.2 修正输入**，521 条的 `gdstype` 重算对照）、
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
4. **模态判定是规则性的。** `truncated_soft_recheck.tsv` 与 `skipped98_modality_recheck.tsv` 的
   分子类别与分辨率由 `!Series_type`、标题与摘要推出，非人工逐条判读。
5. **不覆盖 CELLxGENE 侧。** 该侧可由 1.0.3 数据包内的冻结快照离线重建，无需本包。
6. **`gdsType` 输入在归档中未生效，但已确认不影响已发布取值。** 生成脚本读 esummary 的
   `gdsType` 键，而该接口返回的键名是 `gdstype`，因此归档 raw 表的这一列 521 行**全部为空**，
   `assay_class()` 的第二个输入从未携带信息——这是实现缺陷，已写入主稿“局限与已知问题”。
   为量化其影响，v1.0.1 按正确键名重新取回全部 521 条的 `gdstype`，并用**同一函数**重算：
   423 行的 `molecule` 与 `resolution` **无一改变**，98 条未解析系列的构成也不变
   （mRNA 86／chromatin 7／ncRNA 5；bulk-array 49／bulk-seq 37）。

   **v1.0.2 修正了该对照实验本身的一处缺陷。** v1.0.1 的脚本写成 `r.get("series_type", "")`，
   但归档 raw 表并没有 `series_type` 列（它的第 6 列名为 `gdsType`，且恒为空），所以
   `series_type` 也被传成空串——**同时改变了两个输入**，不能作为"只改变 `gdstype` 一个输入"
   的对照。v1.0.2 改为调 `partD_resfeas_geo.parse_soft_series()` 从系列级 SOFT 缓存读
   `!Series_type` 的真值（与归档流程同一函数、同一缓存），并把该值写入输出表新增的
   `series_type` 列。修正后重跑：`empty series_type = 0`、`empty gdstype = 0`、423 行
   **`changed` 仍为 0**，521 行的 `molecule`／`resolution` 与 v1.0.1 表**逐行相同**。

   **为什么"一项输入为空"结果却不变**：输出表可逐行核对——521 行中 `!Series_type` 与
   `gdstype` 的取值**逐行相同**，仅分隔符不同（`!Series_type` 用 ` | `、`gdstype` 用 `; `），
   两者各只有 22 种取值；规范化分隔符后 **521/521 相等**。因此结论不是"`gdstype` 无关紧要"，
   而是**归档实际使用的 `!Series_type` 已经携带了与 `gdstype` 完全相同的信息**——缺陷之所以
   没有后果，是因为存在冗余；这属于幸运，不能推广为"该输入可以省略"。逐条对照见
   `geo/gdstype_recompute_20261017.tsv`，生成脚本为 `scripts/gdstype_recompute.py`；
   该脚本需访问 NCBI 补取 `gdstype`，重跑所得会随源库刷新而变。
   **本包的两张模态表复现的仍是发布当时的实际输入**（`gdsType` 为空），不是补算后的取值。

## 4 复现方式

**运行前提（v1.0.2 起写清楚，因为 v1.0.1 的说明有歧义）。** 本包 `scripts/` 下的脚本是在
**代码归档仓库**的 `scripts/` 目录里开发的：除 `gdstype_recompute.py` 外，其余脚本都用
`ROOT = Path(__file__).resolve().parents[1]` 定位仓库根，再从 `ROOT / "scripts"` 引入
`partD_resfeas_geo.py`。若直接在本包的 `audit_artifact/scripts/` 下运行，`parents[1]` 会指向
`audit_artifact/`，那里没有 `partD_resfeas_geo.py`，脚本会导入失败。因此：

- **取得代码归档（10.5281/zenodo.23412772）并解包**，把本包 `scripts/*.py` 复制到其
  `scripts/` 目录，然后在该仓库根目录下按下面 1→6 的顺序运行；
- `gdstype_recompute.py` 例外：v1.0.2 起它接受 `--repo` 指定仓库根，**可从任意目录运行**，
  并可保留 `audit_artifact/scripts/` 这一原始位置。

```bash
# 1) 复算检索与截断对照（需可访问 NCBI E-utilities）
python scripts/geo_retrieval_audit.py
# 2) 对 60 条被截断系列做同源 SOFT 复核
python scripts/geo_truncation_soft.py
# 3) 对 98 个未进入发布表的系列做同源复核（v1.0.1 新增；只读 SOFT 缓存，不联网）
python scripts/skipped98_modality_recheck.py
# 4) 重算 521 条并与发布值比对：gdstype 按正确键名重取，series_type 取原 SOFT 的 !Series_type
#    （v1.0.1 新增、v1.0.2 修正输入；需 NCBI 补取 gdstype，需系列级 SOFT 缓存）
python audit_artifact/scripts/gdstype_recompute.py --repo . --out gdstype_recompute_20261017.tsv
# 5) 词表覆盖与敏感性（只读 1.0.3 数据包的四个 TSV）
python scripts/vocab_audit.py
# 6) 检索式词表敏感性（需 NCBI）
python scripts/geo_query_sensitivity.py
```

第 4 步的脚本用 `--repo` 指定代码归档的工作树根，因此它要求该根目录下同时存在
`results/partD_resfeas_geo_series_raw.tsv`、`resource_paper/deposit/partD_resfeas_geo_series.tsv`
与 `data/raw/geo_soft_series/`。第 2、3、4 步还需要系列级或样本级 SOFT 缓存
（`data/raw/geo_soft_series`、`data/raw/geo_soft`），该缓存**不在本包内**，也不在数据包或
代码归档内；缺缓存时脚本会逐条报告缺失的 GSE。`geo_truncation_soft.py` 把样本级 SOFT 缓存
重定向到 `out/geo_retrieval_audit/soft_cache_scratch`，不写回归档缓存。

## 5 与已发布数据包的关系

- 数据包 1.0.3（10.5281/zenodo.23436090）为**纯文档同步版**：其四个表格记录与 1.0.2、1.0.1 **逐字节相同**，
  差异只在 `CHANGELOG.md`、`README.md`、`datapackage.json` 与 `manifest.json` 的描述文字；本包的加入不改变这一点。
- 本包若被合并进数据包，按主稿的语义化版本规则属 **minor** 升级（新增文件、不改既有字段语义）。
- **本包不是 1.0.3 数据包的一部分**，也不在 1.0.3 的数据 DOI 覆盖范围内；它单独公开于代码仓库的固定标签下，
  因此引用时应同时给出标签名与文件相对路径，并核对 `MANIFEST.sha256`。
- 1.0.2 数据包 `datapackage.json` 中 `G3_powered_control` 的字段描述曾写作 *design-eligibility gate*，
  与主稿采用的“未验证的设计线索闸门”不一致；该措辞已在 **1.0.3 版**统一改为
  *unvalidated rule-based design signal*，属**语义修正、不改任何取值**（见主稿“局限与已知问题”第 2 条）。
- 1.0.3 版字典中 `depth_ratio` 的一句描述仍写作 *control arm is shallower*，把基因检出复杂度误写为测序深度；
  该措辞**尚未修正**，主稿“局限与已知问题”第 8 条已如实披露，修正需另一版文档补丁。

## 6 如何取得与核对本包

```bash
# 取得（只含本包，不含代码归档）
git clone --branch audit-artifact-v1.0.2 --depth 1 \
  https://github.com/yuanyuan20020329/sgrg-uc-testability-audit.git sgrg-audit
cd sgrg-audit/audit_artifact

# 或直接下载 zip
# https://github.com/yuanyuan20020329/sgrg-uc-testability-audit/archive/refs/tags/audit-artifact-v1.0.2.zip

# 核对（在 audit_artifact/ 目录下）
sha256sum -c MANIFEST.sha256     # Linux
shasum -a 256 -c MANIFEST.sha256 # macOS
```

30 行全部 `OK` 即与主稿引用的字节一致。若需核对主稿结论，建议顺序：
先看 `geo/diff_*.tsv` 与 `geo/uids_*.tsv`（成员关系），再看 `geo/truncation_comparison.tsv`
（截断组构成），最后看 `geo/gse107593_*.tsv`（单条反事实实例）。
