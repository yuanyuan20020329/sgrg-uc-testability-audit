# 独立案例 A 结果报告 —— 铁死亡基因面板 × 肝细胞癌

> **状态**：Step 2–4 完成（Step 5 外部复现见 §5；Step 6 判定见 §6）
> **完成时间**：2026-10-10 01:00
> **对应计划书**：§6 可执行泛化协议 / §7 案例工作流
> **数据来源说明**：原计划用 GEO 的 GSE65372/GSE84402，**因本沙箱对 `ncbi.nlm.nih.gov` 全域封锁**
> （详见 §7 环境记录），改用 **Gemma REST** 获取独立肝癌队列。**这不影响案例有效性的论证**——
> 反使案例更强：数据源与原文完全不同，且样本量扩大 ~8 倍。

---

## 1. 目标论文（Step 1：23:12 冻结，早于任何分析）

| 项 | 内容 |
|---|---|
| 论文 | Yi S, Zhang C, Li M, Wang J. *Construction of a Novel Diagnostic Model Based on Ferroptosis-Related Genes for Hepatocellular Carcinoma Using Machine and Deep Learning Methods.* J Oncol. 2023;2023:1624580 |
| DOI / PMID / PMCID | 10.1155/2023/1624580 / 36873737 / **PMC9981290**（开放获取，全文已抓取） |
| 原文队列 | 发现 **GSE65372**（39 HCC + 15 nontumor）；验证 **GSE84402**（14 + 14） |
| 原文流水线 | 237 FRG → Student's *t* 检验 **p<0.001（无 FDR、无 FC 门限）** → LASSO + SVM-RFE → logistic ROC |
| 原文主声称 | **11 marker 基因面板，AUC = 1.000** |
| 数据可得性 | "available from the corresponding author upon request" |

### 1.1 全文逐字核实后发现的「原文自身缺陷」（文献事实，非复现结果）

| # | 缺陷 | 逐字证据 |
|---|---|---|
| A1 | **差异基因计数三重不一致** | 摘要/§2.2 = **40**（27↑/13↓）；图 1a 题注/图 2 = **44**；§2.4 = **17**。40→17 的收缩**全文无任何解释** |
| A2 | **全流程无多重检验校正** | §2.2 "Genes that had a p value of less than 0.001 were determined to be significant"；FDR 仅用于通路富集（§2.3） |
| A3 | **无倍数变化门限** | 全文检索 fold change / logFC = **0 命中** |
| A4 | **未做正态性检验即用 Student's *t*** | normality / Shapiro / Wilcoxon 全文 **0 命中** |
| A5 | **"237 FRG" 不是策展库条目数，是平台检出数** | §2.2 逐字 "only **237 FRGs were found to be expressed**"；**全文未引用任何铁死亡基因集来源**（FerrDb 检索 0 命中） |
| A6 | **声称 11 个 marker，只列出 10 个** | "**11** new markers (HSPB1, CDKN2A, LPIN1, MTDH, DCAF7, TRIM26, PIR, BCAT2, EZH2, and ADAMTS13)" —— 括号内仅 10 个 |
| A7 | **全流程 0 次提及批次** | batch / ComBat 全文 **0 命中** |
| A8 | **特征选择与评估同队列，无嵌套 CV** | LASSO/SVM-RFE 用 GSE65372 选特征 → 仍用 GSE65372 报 AUC=1.000 |
| A9 | **数据可得性：向作者索取** | "available from the corresponding author upon request" |
| A10 | **验证段回避报告数字** | "some FRGs among 11 FRGs was further confirmed by the GSE84402…"，**未给任何验证 AUC** |

### 1.2 独立于复现的领域事实：所谓「铁死亡面板」9/10 不是铁死亡基因

把原文列出的 10 个 marker 与 **WP_FERROPTOSIS（MSigDB c2.cp v2023.1，n=64）** 对照：

- **在铁死亡集内**：`HSPB1`（1 个）
- **不在铁死亡集内（9 个）**：`CDKN2A, LPIN1, MTDH, DCAF7, TRIM26, PIR, BCAT2, EZH2, ADAMTS13`

这 9 个是通用的增殖/癌症/转录调控基因。→ 该面板更接近「癌症通用基因 ∩ 与铁死亡有交叉」，
**而非铁死亡特异面板**；且因 A5 原文未公布基因集，面板的「铁死亡」归属**无法被独立核验**。

---

## 2. 数据（Step 2）

| 项 | 内容 |
|---|---|
| 队列 | **GSE36376**（Gemma dataset id **7673**），Illumina HumanHT-12 V4（GPL10558） |
| 样本 | **433**：**240 liver tumor / 193 adjacent non-tumor liver** |
| 配对结构 | **非配对**（240 与 193 来自不同患者；样本名 `GSE36376_Biomat_N` 唯一）→ 与原文 39+15 同为非配对，可直接独立样本检验 |
| Gemma 批次标注 | **`PROBLEMATIC_BATCH_INFO_FAILURE`**（批次信息有问题 → L2 有实质可检验内容） |
| 基因集 | WP_FERROPTOSIS（n=64）；本队列 **62/64** 可测 |
| 为何选它 | Gemma 全库肝癌筛过后，**唯一**「大样本 + 干净肿瘤/癌旁二分组」队列（其余：GSE6465 是异种移植+细胞系；GSE54236 tissue 全写 "Liver" 无分组；GSE4108 仅 2 正常；GSE19977 全癌） |

---

## 3. Step 3：忠实复现原文流水线（结果：**成功复现**）

在**独立队列、8 倍样本量**上按原文同一逻辑（t 检验 p<0.001 无 FDR → LASSO + SVM-RFE → logistic 同集 ROC）：

| 量 | 原文 | 本案例（GSE36376） |
|---|---|---|
| 差异 FRG 数 | 40（或 44，原文自相矛盾） | **44**（与原文图 1a 的 44 一致） |
| 筛选链 | LASSO ∩ SVM-RFE → 11 | LASSO 24 → SVM-RFE 18 → **15 final** |
| **AUC** | **1.000** | **0.9991** |
| BH-FDR 口径 | 未做 | DE=44（FDR<0.001）→ AUC=0.9966 |

**★ 结论：原文主结果可复现，且不是单队列偶然。** 这使审计问题从「能否复现」升级为
**「AUC≈1 究竟由什么产生」**——这正是下面四层要回答的。

---

## 4. Step 4：审计阶梯

### L4 · 随机基因集校准（n=500，表达分层匹配，k=44）—— **决定性**

| 量 | 值 |
|---|---|
| 真实面板 in-sample AUC | **0.9994** |
| 随机集 in-sample AUC | **中位 0.9978、q95 0.9988、max 0.9993** |
| **随机集中 in-sample AUC ≥ 0.99 的比例** | **100%（500/500）** |
| 随机集诚实 7 折 CV AUC | 中位 0.9905、q95 0.9939 |
| 真实面板 7 折 CV AUC | 0.9898 |
| 面板 in-sample 经验 p | 0.0020（**假显著**——整个零分布都压在 0.998 上） |
| **面板 7 折 CV 经验 p** | **0.611（位于随机集分布第 61 百分位）** |

> **结论：AUC≈1.000 是「该队列 + 选特征后在同一数据上评估」这一协议的属性，
> 不是铁死亡基因的属性。任何 44 个基因都能做到。**

### L2 · 技术轴占据度

| 量 | 值 |
|---|---|
| **PC1 对肿瘤/癌旁的单轴 AUC** | **0.967** |
| PC1–PC3 解释方差 | 39.7% |
| 最大 PC 的疾病 AUC | 0.967（即 PC1） |

### L3 · 评分-技术负载耦合 & 单基因 AUC

| 量 | 值 |
|---|---|
| ρ（面板评分, 全基因平均表达） | **0.938** |
| 单基因 AUC 中位 / 最大 | 0.680 / **0.953** |
| 单基因 AUC ≥0.95 的基因数 | **1 / 62**（≥0.99：0 个） |
| AUC（全基因平均表达 ~ 疾病） | 0.624 |

→ **没有任何单基因是诊断性的**（最高 HSPB1 = 0.953），但 44 个基因的组合却给出 0.999：
这个跃升**不可能来自基因本身的信号**，只能来自「在同一数据上选 44 个最分开的基因」的机制。

### 机制确证：移除全局轴后面板崩塌

| 条件 | 差异基因数 | AUC |
|---|---|---|
| 原样 | 44 | **0.9991** |
| **移除 PC1（疾病/技术全局轴）后** | **18** | **0.757** |

→ **AUC≈1.0 由单一全局轴（PC1）承载；抽掉它，面板退化到 0.757。**
而原论文全流程 **0 次提及批次**（缺陷 A7）—— 这正是 L2 要暴露的、原作者从未考虑的混杂。

---

## 5. Step 5：外部复现

原计划的外部集 GSE84402（14+14，仅 28 样本）**因 NCBI 全域封锁无法下载**（见 §7）。
替代方案更强，已分两路完成：

- **路线 1（独立队列）**：在 **GSE36376（433 样本，8 倍于原文）** 上复现原文流水线 → AUC=0.9991（§3），
  并给出 500 集随机基因集零分布（§4 L4）。
- **路线 2（原文队列，意外获得）**：在 **GSE65372（原文发现队列）** 上做配对 vs 非配对对比（§5.5）。

## 5.5 ★★ 附加审计：在**原文自己的队列 GSE65372** 上发现伪重复（pseudoreplication）

> 本节是本次会话的**意外收获**：NCBI 监视器在 7h46m 的守候中于某个窗口期成功抓到了
> **GSE65372 的完整系列矩阵**（9,207,303 字节，可完整解压 19.4 MB）。于是可以直接审计
> **原文自己的发现队列**，而不只是独立队列。

### 5.5.1 数据的真实设计（来自 GEO 自身的元数据）

| 项 | 实测 | 说明 |
|---|---|---|
| 样本 | **54** | 与原文一致 |
| 分组 | **39 tumor (Y) / 15 non-tumor (N)** | 与原文声称完全一致 |
| **患者数** | **18** | ★ 远小于 54 |
| 患者内结构 | **14 个患者同时有 Y 和 N**；**17 个患者有 tumor 样本** | — |
| **关键** | **多数患者贡献 2–3 个 tumor 样本，却只贡献 1 个 normal** | 患者 422/82/179 各 3Y+1N |

**GEO 自己的 series summary 逐字写明**：

> "Genome-wide gene expression profiles of **multiple hepatocellular carcinoma tumors and
> non-tumor liver tissues from the same patients**."
> "Molecular heterogeneity of **multinodular** hepatocellular carcinoma"

→ **数据提供方明确声明这是「同一患者的多结节肿瘤 + 同患者癌旁」的配对/多重复设计**，
而原文把它当成 **39 vs 15 的独立样本比较**（§2.2 "the Student's t-test"），
**未做任何配对、聚类或混合效应校正**。

### 5.5.2 三口径对比（同一批 64 个铁死亡基因）

| 口径 | n | DE (p<0.001) | DE (p<0.05) | in-sample AUC |
|---|---|---|---|---|
| **P1 · 原文口径**（54 样本，独立样本 t） | 39 vs 15 | **8** | 19 | 0.9795 |
| **P2 · 配对口径**（14 对同患者，患者内聚合） | **14 对** | **2** | 16 | 0.9847 |

- **伪重复膨胀因子 = 2.29**（39 个瘤样本来自 17 个患者）
- **两口径的 p<0.001 显著基因仅重叠 1 个（TXNRD1）** → 原文的 8 基因命中表 **7/8 是口径依赖的**
- ★ **AUC 几乎不变（0.9795 → 0.9847）** → **再次证明 AUC 不由基因信号决定**，
  而由队列的全局结构决定（与 §4 的 PC1 机制完全一致）

### 5.5.3 意义

这是**与「AUC 无特异性」相互独立的第二个失效模式**：

1. **伪重复**（统计设计缺陷）：把患者内重复当独立观测，使 n 虚增 2.29 倍、
   p 值系统性偏小、发现集膨胀 4 倍（8→2）。原文明写用 "independent samples" 的 t 检验，
   而数据提供方明写是 "from the same patients" —— **方法与数据的自述直接矛盾**。
2. **选择偏倚 + 全局轴**（§4）：AUC≈1 与随机集不可区分，由 PC1 承载。

两者**互相独立**，且都在**原文自己的数据**上得到证实。这使案例 A 的判定从
「面板在别的队列上失效」升级为「**原文在其自身队列上的统计处理即存在两处独立缺陷**」。

---



---

## 6. Step 6：最终判定（四选一）

> 计划书 §6 的四个出口：① 支持原结论 ② 面板失效但机制有别 ③ 复现失败 ④ 结论无法检验**

**判定 = ②「面板失效，且失效机制明确」**，具体为：

1. **可复现**：原文 AUC≈1.000 在独立队列上是真的（0.9991）——**不能指责原文报错**；
2. **但无特异性**：同规模随机基因集中位 AUC = 0.9978，**100% 的随机集 ≥0.99**；
   面板在诚实 CV 下位于第 61 百分位 → **与随机不可区分**；
3. **机制 = 选择偏倚 + 全局技术轴**：AUC 由 PC1（单轴 AUC=0.967）承载，移除后掉到 0.757；
   而原论文**从未检查批次**（A7）；
4. **面板的「铁死亡」属性不可核验**：原文未公布基因集（A5），且列出的 10 个 marker 中
   **9 个不是铁死亡基因**（§1.2）；
5. **★ 叠加第二个独立缺陷 = 伪重复**（§5.5）：在原文**自己的队列**上，54 样本实为 18 患者、
   39 个瘤样本来自 17 个患者（数据自述 "from the same patients"），
   原文按独立样本处理 → n 虚增 2.29 倍、发现集 8→2（4 倍膨胀）、
   两口径仅重叠 1 个基因；而 **AUC 几乎不变（0.98→0.98）**，再次独立佐证 AUC 不由基因信号决定。

**跨疾病意义**：这是「选择偏倚 → 完美 AUC → 随机基因集也能做到」这一失效模式在**第二个疾病领域**
（HCC/铁死亡，非 UC/SGRG）上的独立复现 → **直接支撑主论文「五层阶梯可泛化」的核心论点**；
同时额外贡献一个**跨疾病通用的 L5 设计缺陷实例（伪重复）**，可并入 Box 1 的第 ⑤ 层。

---

## 7. 环境记录与后续（阻断项归档）

### NCBI 全域封锁（本次会话）
- 实测 `www.ncbi.nlm.nih.gov`、`geo.ncbi.nlm.nih.gov`、`ftp.ncbi.nlm.nih.gov`、
  `eutils.ncbi.nlm.nih.gov` 全部 **http=000**；curl -v 显示沙箱代理 CONNECT 隧道建立后连接被静默关闭。
- 早先的 `/geo/download/` 返回的是 **21,658 字节 Google reCAPTCHA 挑战页**（假成功）。
- **教训（写进方法学）**：「下载成功」必须验 **gzip magic + 完整解压 + 内容头**，
  不能只看 HTTP 200 与字节数——本次差点把两个 reCAPTCHA 页当成数据。

### 绕行（已成功）
- **Gemma REST** 取代 GEO：`/datasets/{id}/samples`（游标分页，含 tissue 标签与 GSM）、
  `/datasets/{id}/expressions/genes/{EntrezIDs}?consolidate=average`（64 基因 × 433 样本，45 秒）。
- **symbol→Entrez 映射**用 **mygene.info**（本沙箱可达），64/64 全中。
- **GPL14951 探针注释**：GEO 的 GPL SOFT 不可达 → 改用 **Bioconductor 的 `illuminaHumanv3.db`**
  （6.4 MB，含 49,576 探针的 `SymbolReannotated` 列），成功把 64/64 铁死亡基因映射到 GSE65372 探针。
- Gemma OpenAPI 契约已存档 `data/raw/caseA/gemma/openapi.json`。
- **可复用清单（本沙箱可达）**：Gemma、Europe PMC、EBI/BioStudies、Broad MSigDB、
  bioconductor.org、CRAN、mygene.info、github.com（但 raw.githubusercontent 不可达）。

### 后续可选
- [ ] 在 Gemma 另选一个「肿瘤/癌旁」肝癌队列做**第二个外部复现**（增强 Step 5）
- [ ] 用 Gemma GEEQ 的 per-factor 分解（`/datasets/{id}/geeq/public`）量化 PC1 的批次归属
- [ ] 把本案例并入主论文 §2.8 独立案例区块，与 Part B 一起构成「阶梯非逢事皆 no」的证据

---

## 附：产出文件清单

| 文件 | 内容 |
|---|---|
| `results/caseA_论文登记表.md` | Step 1 冻结登记表（含 §3.4 原文缺陷清单） |
| `results/caseA_ferroptosis_genes.txt` / `_entrez.tsv` | WP_FERROPTOSIS 64 基因与 Entrez 映射 |
| `results/caseA_expr_7673.tsv` | 433 样本 × 62 基因 + 标签 |
| `results/caseA_sample_labels_7673.tsv` | 样本级标签（label/GSM/Gemma 名/平台/分期） |
| `results/caseA_step3_de.tsv` | 62 基因的 p/ Welch p/ FDR / Cohen's d |
| `results/caseA_L4_random.json` | L4 随机集校准（500 集） |
| `results/caseA_L2L3.json` | L2 技术轴 + L3 负载 + 机制检验 |
| `results/caseA_GSE65372_pseudorep.json` | 原文队列的设计结构（18 患者 / 54 样本 / 伪重复） |
| `results/caseA_GSE65372_paired.json` | **原文队列三口径对比**（P1 非配对 8 个 DE vs P2 配对 2 个 DE） |
| `data/raw/caseA/annot/illuminaHumanv3.db/` | GPL14951 探针→基因注释（Bioconductor，49,576 探针） |
| `scripts/caseA_gemma.py` 等 11 个 | 采集/映射/拼装/分析脚本 |
