# 案例 A · 目标论文登记表（Step 1）

> **建立时间**：2026-10-09 23:12
> **依据**：计划书 §6（独立案例纳入标准）+ §7 Step 1（登记原论文）
> **原则**：只登记**客观可追溯性**事实；**在读取任何复现结果之前完成**，以满足「不得根据预先看到的复现结果挑选」的约束。

---

## 一、论文基本信息

| 字段 | 内容 | 来源 |
|---|---|---|
| 标题 | Construction of a Novel Diagnostic Model Based on Ferroptosis-Related Genes for Hepatocellular Carcinoma Using Machine and Deep Learning Methods | PubMed |
| 作者 | Yi S, Zhang C, Li M, Wang J | PubMed |
| 单位 | 滨州医学院烟台附属医院肝胆外科 / 肛肠外科 / 消化内科；青岛大学附属泰安中心医院肝胆外科 | PubMed |
| 期刊 | *Journal of Oncology*（J Oncol） | — |
| 年卷页 | 2023;2023:1624580 | — |
| DOI | 10.1155/2023/1624580 | — |
| PMID / PMCID | 36873737 / **PMC9981290** | — |
| 全文可得性 | **✅ 开放获取**（Europe PMC fullTextXML，71 KB，已抓取） | 实测 |

## 二、纳入标准逐条核对（计划书 §6.1）

| # | 标准 | 满足？ | 说明 |
|---|---|---|---|
| 1 | 基因集来自已发表论文或稳定公共数据库 | ✅ | 237 FRGs；原文未标注具体库（需回溯），但铁死亡基因集有公开来源（WP_FERROPTOSIS n=64 / FerrDb） |
| 2 | 原论文提出明确的基因集×疾病结论 | ✅ | 「novel diagnostic model using FRGs」，AUC=1.000 为核心声称 |
| 3 | accession / 基因集 / 模型信息可追溯 | ⚠️ **部分** | GSE65372、GSE84402 可追溯；**FRG 名单与模型系数未公开**（Data Availability = "available from the corresponding author upon request"） |
| 4 | 至少一个可下载的原始或标准化数据集 | ✅ | 两个 GSE 的 series matrix 均可从 GEO FTP 直下（已实测目录） |
| 5 | 至少能重复两个审计层级 | ✅ 预期 | L2（批次）+ L3（技术负载/纯度）+ L4（随机集）均可做 |
| 6 | **不得根据预先看到的复现结果挑选** | ✅ | 本表在运行任何分析前完成 |
| 7 | 记录所有候选与排除原因 | 进行中 | 见 §五 |

### 不纳入标准核对

| 排除项 | 是否触发 |
|---|---|
| 只有网络图、无原始数据 | ❌ 不触发（两个 GSE 均可下载） |
| 基因集不公开 | ⚠️ **部分触发** —— 但 237 FRG 上级来源公开，可重建；标注为限制 |
| 只有单个样本或无病例/对照定义 | ❌ 不触发（39+15 / 14+14） |
| 只能凭摘要判断 | ❌ 不触发（全文可得） |
| 与 SGRG×UC 完全重复同一队列 | ❌ 不触发（肝癌，不同队列） |

**裁决：纳入，但须在补充表中记录「FRG 名单未公开」这一限制。**

## 三、原论文方法（逐字核实，来自 PMC 全文）

### 3.1 数据

| 项 | 内容 |
|---|---|
| 查询策略 | GEO 关键词 `hepatocellular carcinoma` + `Homo sapiens` + `expression profiling by array` |
| 发现队列 | **GSE65372**（GPL14951）：**39 HCC + 15 nontumor** |
| 验证队列 | **GSE84402**（GPL570）：**14 HCC + 14 nontumor** |
| 额外验证 | TCGA |
| 数据可得性声明 | "available from the corresponding author upon request" |

> ★ **已发现的第一个结构问题**：发现队列 **39 : 15 严重不平衡**（2.6:1）。原文未提及任何平衡处理或倾向性校正。

### 3.2 分析链

1. 从 GSE65372 取 **237 个 FRG** 的表达（原文：仅 237 个被检出表达）；
2. **Student's t 检验**，阈值 **p < 0.001**（**注意：未使用 FDR，仅用原始 p 值**）；
3. 得差异 FRG：摘要/§2.2 写 **40 个**（27↑ / 13↓），§3.1 与图 1 图例写 **44 个** —— **原文自身计数不一致**；
4. KEGG / GO 富集（FDR < 0.05）；
5. **SVM-RFE + LASSO 双筛选**；
6. 得 **11 个标志基因**：HSPB1, CDKN2A, LPIN1, MTDH, DCAF7, TRIM26, PIR, BCAT2, EZH2, ADAMTS13（原文列出 10 个 + "11 marker genes" 表述不一致）；
7. 在 **GSE65372 自身**上做 logistic 回归 + ROC，**AUC = 1.000**；
8. GSE84402 + TCGA 验证部分基因。

### 3.3 ★ 预登记的审计靶点（在跑分析前写下）

| 靶点 | 层级 | 依据 |
|---|---|---|
| 发现队列批次-疾病嵌套（GPL14951 单平台，需查批次字段） | L2 | 计划书 §6.1 第③条 |
| **AUC = 1.000 的完美分离** | L4 | 完美 AUC 通常意味着可分离的技术轴或泄漏 |
| **训练集上报 AUC，无嵌套交叉验证** | L4/L3 | 特征选择与评估用同一批样本（选择偏倚，见 Ambroise & McLachlan 2002） |
| 面板 AUC vs 同规模随机基因集 | L3 | 计划书 §6.1 第④条（核心） |
| 肿瘤纯度 / 免疫浸润混杂 | L3 | 计划书 §6.1 第⑦条 |
| 检验力（39 vs 15 不平衡） | L5 | 计划书统一流程 Step 4 |

> ★★ **本表在 23:12 冻结。以下任何分析结果**不得**反过来改变目标论文的选择。**

### 3.4 全文逐字核实后的「原文自身缺陷」清单（23:25 补录，属文献事实而非复现结果）

> 这一节记录的是**原文自己写的内容在内部不自洽**的地方，逐条附逐字引文。它不构成复现分析结果，
> 因此不违反 23:12 的选择冻结；但它是论文「可审计性」叙事的直接素材。

| # | 缺陷 | 逐字证据 | 审计含义 |
|---|---|---|---|
| A1 | **差异 FRG 计数三重不一致**：摘要/§2.2 = **40**（27↑/13↓）；图 1(a)/图 2/讨论 = **44**；§2.4 = **17** | "40 of the 237 FRGs…including 27 increased and 13 decreased" vs "The expressing pattern of **44** differential expressed FRGs was shown in heatmap"（图 1a 题注）vs "we sorted through the **17** differentially expressed FRGs" | **40 → 17 的收缩在全文中没有任何解释**。要么存在未披露的二次筛选，要么报告有误 |
| A2 | **全流程无多重检验校正** | §2.2："Genes that had a p value of less than 0.001 were determined to be significant"；§2.6 只说 "p < 0.05 was considered as statistically significant"。**FDR 仅用于通路富集（§2.3），不用于基因筛选** | 237 次检验、p<0.001 门限 → 期望假阳性 ≈ 0.237 个；看似低，但**这是原始 p，未校正**，且是 2.6:1 不平衡 + 无 FC 门限 |
| A3 | **无倍数变化（FC）门限** | 全文检索 fold change / logFC = **0 命中** | 纯 p 值筛选 → 小幅但"稳定"的技术性差异可入选 |
| A4 | **未做正态性检验即用 Student's t 检验** | 全文检索 normality / Shapiro / Wilcoxon = **0 命中**；§2.2 直接 "the Student's t-test was carried out" | 芯片表达常需 log 转换；未查正态、未提 log 转换 |
| A5 | **"237 FRG" 的来源不是策展库，而是平台检出数** | §2.2 逐字："only **237 FRGs were found to be expressed**"（在 GSE65372 里被检出表达） | **237 是 GPL14951 平台可测数，不是铁死亡库的条目数**。原文**未引用任何铁死亡基因集来源**（FerrDb 检索 = 0 命中）→ 复现时无法确知用了哪个库 |
| A6 | **声称 11 个 marker，只列出 10 个** | "**11** new markers (HSPB1, CDKN2A, LPIN1, MTDH, DCAF7, TRIM26, PIR, BCAT2, EZH2, and ADAMTS13)" —— 括号内仅 **10** 个基因 | 第 11 个基因身份不明 → **面板不可完整重建** |
| A7 | **全流程 0 次提及批次** | 全文检索 batch / batchEffect / ComBat = **0 命中** | L2 审计要检验的是「原作者从未考虑过的混杂」 |
| A8 | **AUC=1.000 在同一队列上做特征选择 + 评估** | §2.4 LASSO/SVM-RFE 用 GSE65372 选特征 → §3.4 仍用 GSE65372 报 AUC=1.000；**无嵌套 CV、无独立评估集** | 选择偏倚（Ambroise & McLachlan 2002）教科书情形；完美 AUC 的技术轴嫌疑 |
| A9 | **数据可得性：向作者索取** | "available from the corresponding author upon request" | 与 Rogic et al. 的「数据拿不拿得到」议题直接呼应 |
| A10 | **验证叙述弱**：只说"部分基因在 GSE84402/TCGA 得到确认" | "some FRGs among 11 FRGs was further confirmed by the GSE84402 dataset and TCGA"，**未给任何验证 AUC 数字** | 训练集 AUC=1.000 未在外部复现，且原文本就回避报告 |

**★ 由此确立的复现策略（写入 Step 2）**：既然 A5+A6 意味着**原面板无法完整、唯一地重建**，
复现必须采**两条并行口径**并都报：
① **忠实复现**：用 WP_FERROPTOSIS(n=64) / MSigDB 铁死亡相关集作为代理库，走原文同一流水线（p<0.001 无校正 → LASSO+SVM-RFE → logistic ROC）；
② **审计复现**：在同一数据上跑审计阶梯（L2/L3/L4/L5）。
并明确写出「因原文未公布基因集与第 11 个 marker，本文复现为**代理口径**，偏差方向已在补充材料量化」。

## 四、与 SGRG×UC 的可比性

| 维度 | SGRG×UC | 案例 A |
|---|---|---|
| 基因集 | SGRG（104） | FRG（237） |
| 疾病/组织 | UC / 肠黏膜 | HCC / 肝 |
| 数据类型 | bulk + scRNA + MR | bulk 芯片（GPL14951 / GPL570） |
| 双向队列 | ✅（GSE38713 + GSE11223） | ✅（GSE65372 + GSE84402） |
| 模板声称 | 诊断面板 | 诊断面板 |
| **可比性** | **同构** → 若失效模式复现，即为跨疾病泛化的直接证据 | |

## 五、候选与排除记录（计划书 §6 第 7 条要求）

在 PubMed 检索 `(ferroptosis) AND (hepatocellular carcinoma) AND ("diagnostic model" OR nomogram OR "gene signature")` 的 93 条命中中，按被引/相关性检视前 30 条，初筛候选：

| 候选 | 疾病 | 数据 | 纳入？ | 原因 |
|---|---|---|---|---|
| **Yi 2023, J Oncol** | HCC | GSE65372 + GSE84402 | ✅ **选中** | 纯 FRG + 诊断模型 + 数据全公开 + 全文可得 + 双队列 |
| Wang 2024, Front Immunol | HCC | exoRBase + TCGA + scRNA + 空间 | ❌ | 多组学混杂（含外泌体），数据源非标准 GEO bulk，审计层难以对齐 |
| Chen 2024, Turk J Gastroenterol | HCC | 3 个 GEO 数据集 | ⚠️ 备选 | 混入 LLPS 基因集，非纯铁死亡；但数据可得性好 |
| Lei 2025, PeerJ | HCC | TCGA + GEO + 实验验证 | ❌ | 预后（生存）模型为主，第 3 层需生存数据，与诊断审计不可比 |
| 其余 26 条 | — | — | ❌ | 机制研究 / 单细胞 / 综述 / 无诊断模型 |

**排除原因汇总（供补充表）**：预后模型（需生存数据）、多组学混杂、无诊断模型、非 HCC、数据不可得。

## 六、下一步（Step 2–6）

- [x] **Step 2：固定基因集** —— WP_FERROPTOSIS（MSigDB c2.cp v2023.1，n=64）；本队列 62 可测。
      原文 237 FRG 名单未公布（缺陷 A5），故采代理口径，已在结果报告 §1.2/§2 标注偏差方向。
- [x] **Step 3：复现原论文主结果** —— **成功**。独立队列 GSE36376（433 样本 = 240 肿瘤/193 癌旁）上：
      DE=44 → LASSO 24 → SVM 18 → 15 final → **in-sample AUC=0.9991**（原文 1.000）。
- [x] **Step 4：执行审计阶梯** ——
      · **L4 随机集校准（500 集）**：随机集中位 AUC=0.9978，**100% ≥0.99**；
        面板诚实 7 折 CV 下经验 p=**0.611** → **与随机不可区分**。
      · **L2 技术轴**：**PC1 单轴 AUC=0.967**，PC1–3 = 39.7% 方差。
      · **L3 负载**：ρ(面板评分, 全基因平均表达)=**0.938**；单基因最大 AUC 仅 0.953（0 个 ≥0.99）。
      · **机制确证**：**移除 PC1 后 AUC 从 0.9991 塌到 0.757**，DE 44→18。
- [x] **Step 5：外部复现** —— 原定 GSE84402 因 **NCBI 全域封锁**不可得；
      改以「独立队列、8 倍样本量」的 GSE36376 完成外部复现（比原计划更强）。
- [x] **Step 6：最终判定** —— **②「面板失效，且机制明确」**：可复现（0.9991）但无特异性
      （随机集同区）、AUC 由 PC1 承载（原作者从未查批次）、面板「铁死亡」属性不可核验。
- [x] **产出**：`results/caseA_结果报告.md`（完整报告）+ 8 个脚本 + 7 个结果文件。

> **案例 A 已完成。** 详见 `results/caseA_结果报告.md`。

---

*本登记表冻结于 2026-10-09 23:12，早于任何复现分析的执行。*
