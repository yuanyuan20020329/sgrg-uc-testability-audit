# 案例 B 目标论文登记表（冻结）

> ★★ **本表于 2026-10-10 20:38 冻结**，冻结时点**尚未运行任何复现或审计统计**。
> 冻结后**不得**因任何后续结果更改目标论文的选择。
> 本表中「已核验事实」一栏均为**对论文文本与公开元数据的记录**（等同案例 A 的做法），不属于分析结果。

---

## 一、目标论文

| 项 | 内容 |
|---|---|
| 标题 | Integrated analysis of single-cell RNA sequencing and bulk transcriptome data identifies a pyroptosis-associated diagnostic model for Parkinson's disease |
| 作者 | Wang L, Qin Y, Song J, Xu J, Quan W, Su H, Zeng H, Zhang J, Li J, Chen J |
| 期刊 / 年 / 卷页 | **Scientific Reports** 2024; **14**(1):28548 |
| DOI | `10.1038/s41598-024-80185-9` |
| PMID / PMCID | **39558055** / **PMC11574289** |
| 开放获取 | 是（全文 XML 已存 `data/raw/caseB/paper/PMC11574289_fullText.xml`） |
| 疾病 / 基因集族 | 帕金森病（PD） / **焦亡（pyroptosis）** |
| 模态 | scRNA-seq + bulk 微阵列（双模态） |

**核心声称**：由 2 个焦亡相关基因 **POLR2K + TIMM8B** 构成诊断评分
`diagnostic score = −0.029 × POLR2K − 0.018 × TIMM8B`；
训练集（bulk）**AUC = 0.823**（POLR2K 0.742、TIMM8B 0.777），
外部验证 **GSE49036 AUC = 0.706**。

---

## 二、已核验事实（论文文本 + 公开元数据，非分析结果）

### 2.1 数据构成（★ 论文方法段原文）

> "The single cell datasets GSE140231, GSE157783 and GSE204796 were downloaded from the NCBI GEO database...
> There were **12 healthy samples in GSE140231**, **6 healthy and 5 disease samples in GSE157783**,
> and **9 healthy samples in GSE204796**. Therefore, **a total of 5 disease and 27 healthy samples**
> were included in this study as the subsequent model validation dataset...
> Bulk datasets GSE7621, GSE20164, GSE20163, and GSE8397 ... **41 healthy and 59 disease samples**
> were included in the subsequent model construction analysis..."

| 臂 | 样本数 | 来源数据集 | 每集贡献 |
|---|---|---|---|
| **病例（PD）** | **5** | **仅 GSE157783** | 5 |
| **对照（healthy）** | **27** | GSE140231 + GSE157783 + GSE204796 | 12 / 6 / 9 |
| bulk 模型构建 | 100 | GSE7621、GSE20164、GSE20163、GSE8397 | 41 H + 59 D |

### 2.2 三个 scRNA-seq 来源的身份（按 accession 反查文献确认）

| 数据集 | 原始研究 | 期刊 / 年 | 材料性质 |
|---|---|---|---|
| **GSE140231** | A single-cell atlas of the human substantia nigra… | Nat Commun 2020 | **尸检人脑黑质**（成人脑组织） |
| **GSE157783** | Single-cell sequencing of human midbrain reveals glial activation and a Parkinson-specific neuronal state | Brain 2022 | **尸检人中脑**（PD + 对照，本论文 5 例病例**全部**来自此集） |
| **GSE204796** | Human midbrain dopaminergic neuronal differentiation markers predict cell therapy outcomes in a Parkinson's disease model | J Clin Invest 2022 (PMID 35700056) | ★ **hPSC 诱导分化的中脑多巴胺能神经元培养物**（**体外**细胞分化体系，非人脑组织） |

### 2.3 bulk 队列的批次状态（Gemma 公开策展字段，零下载）

| GSE | Gemma id | n | batchEffect |
|---|---|---|---|
| GSE7621 | 549 | 25 | BATCH_CORRECTED_SUCCESS |
| GSE20164 | 3887 | 11 | BATCH_CORRECTED_SUCCESS |
| GSE20163 | 3886 | 17 | BATCH_CORRECTED_SUCCESS |
| GSE8397 | 873 | 47 | **PROBLEMATIC_BATCH_INFO_FAILURE** |
| GSE49036（验证） | 11617 | 28 | **BATCH_EFFECT_FAILURE** |

### 2.4 论文方法学口径

- 未提及任何**批次校正检验**（论文使用 Harmony 做 sc 整合，但未报告整合前后批次可分性，也未检查 batch-disease 嵌套）；
- 训练集 AUC **在训练数据上同集评估**（无嵌套交叉验证、无独立训练/测试划分）；
- 焦亡相关 DEG 阈值：sc 数据 `p < 0.05 且 |logFC| > 0.1`（**极宽松**），得 **3088** 个「焦亡相关 DEG」；
  bulk 数据另得 **1261** 个焦亡相关基因；
- 特征筛选：LASSO（24 变量）+ SVM-RFE + 随机森林 → 交集 → 6 → 多因素 logistic → **2** 个诊断基因。

---

## 三、候选池与排除原因（计划书 §六原则 7）

检索通道：Europe PMC REST（NCBI 全域封锁，见 §五）。
检索式与全量结果：`results/caseB_candidates.json`、`results/caseB_candidates2.json`、
筛选记录 `results/caseB_target_screen.json`。

**同体裁（疾病 × 面板 × 诊断模型）命中 17 篇**，其中：

| 候选 | 期刊/年 | 入选/排除 | 理由 |
|---|---|---|---|
| **Wang 2024（本目标）** | Sci Rep 2024 | ✅ **入选** | 双模态（sc+bulk）、基因集与 accession 全可追溯、有明确 AUC、OA 全文可获 |
| "Identification of Autophagy-Related Biomarkers and Diagnostic Model in AD" | Genes 2024 | ❌ 排除 | 无单细胞成分 → **无法做 L1 对照来源审计**（案例 B 的承重层） |
| "Identification and validation of pyroptosis-related genes in AD based on multi-transcriptome…" | Front Aging Neurosci 2025 | ⚠️ 备选 | 有 snRNA-seq（GSE242902），但未报告 AUC、面板大小不明 |
| "Using Blood-Based Biomarkers… Pyroptosis-Associated Signature Model"（AD 血液） | Mol Neurobiol 2025 | ❌ 排除 | **非 OA**，方法细节不可核验 |
| "Application of machine learning reveals diagnostic biomarkers related to pyroptosis in AD" | J Alzheimers Dis 2025 | ❌ 排除 | **非 OA** |
| "Identification of Pyroptosis-Related Molecular Subtypes and Diagnostic Model… MDD" | Mol Biotechnol 2025 | ❌ 排除 | **非 OA** |
| "Identification and clinical validation of autophagy-related genes… PD" | Clin Neurol Neurosurg 2026 | ❌ 排除 | **非 OA**；单队列 |
| "Integrated Bioinformatics Analysis of … Mitochondrial Autophagy… AD" | J Neurochem 2026 | ❌ 排除 | **非 OA** |
| "…lncRNA LINC00269 … PD via pyroptosis pathways" | Eur J Med Res 2024 | ❌ 排除 | 研究对象是 lncRNA，非基因集诊断面板；无单细胞 |
| 其余 8 篇 | — | ❌ 排除 | 机制/治疗类，非「面板 × 诊断模型」体裁 |

**排除主因统计**：非 OA（5）、无单细胞（2）、非诊断面板体裁（10）。

---

## 四、预注册审计靶点（跑之前写定）

1. **L1 对照来源**：27 个「healthy」中，**21 个（77.8%）来自 GSE140231 与 GSE204796——这两个集合不含任何 PD 病例**；
   且 GSE204796 为**体外 hPSC 分化培养物**，与尸检脑组织不同源。→ 质疑「健康对照」的生物学同质性。
2. **L2 批次-疾病嵌套**：病例臂 100% 来自 GSE157783，对照臂横跨三集 → **疾病与数据集完全共线**；
   预期 Cramér's V ≈ 1、batch-only AUC 很高。bulk 侧训练集全部经批次校正、验证集 `BATCH_EFFECT_FAILURE`。
3. **L3 技术负载**：诊断评分是 2 基因加权和 → 检验评分与**全局表达水平/检出率**的耦合（案例 A 同法），
   以及 POLR2K/TIMM8B 是否只是「管家/线粒体样」深度代理基因。
4. **L4 随机集校准**：2 基因极短面板 → 与同规模随机基因对比较训练集 AUC（案例 A 同法）。
5. **L5 检验力**：sc「验证」仅 **5 例病例**；训练集同集评 AUC=0.823（无嵌套 CV）；验证 AUC 掉到 0.706。
6. **领域层 L1（第二域泛化）**：CELLxGENE 神经退行全域的对照来源污染率
   （与 Part D 的 UC 结果对照，检验「脑图谱复用是否更严重」）。

---

## 五、通道与环境（冻结时点已知）

- **NCBI 全域封锁**（www/geo/ftp/eutils 全 http=000，`/geo/download/` 返回 reCAPTCHA 页）→
  三个 scRNA-seq 原始数据**无法下载**；GSE 的身份改用**文献反查 + 公开元数据**确认。
- **可行通道**：Gemma REST（bulk 表达矩阵 + 批次策展字段）、CELLxGENE 全库缓存（来源审计）、
  Europe PMC（全文/引文）、mygene.info（symbol→Entrez）、Bioconductor（探针注释）。
- 结论：案例 B 的 **L1（sc 侧）在「设计/元数据」层面完成**，
  **L2/L3/L4/L5 在 bulk 侧（Gemma）完成** —— 与案例 A（全在 bulk）互补且不重复。

---

## 六、下一步（Step 2–6）

- [ ] Step 2：固定基因集（论文未公布 PRG 名单 → 以 MSigDB **REACTOME_PYROPTOSIS** 重建并标注差异）；
      固定 POLR2K(5440) / TIMM8B(26521) Entrez ID。
- [ ] Step 3：用 Gemma 的 5 个 bulk 队列复现诊断评分与 AUC（目标 ≈0.823 训练 / 0.706 验证）。
- [ ] Step 4：执行审计阶梯 L1(领域+设计) / L2 / L3 / L4 / L5。
- [ ] Step 5：外部复现（GSE49036 即论文自己的验证集；另可跨病种）。
- [ ] Step 6：最终判定（四选一）。
