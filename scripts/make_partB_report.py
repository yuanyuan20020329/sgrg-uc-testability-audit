# -*- coding: utf-8 -*-
"""Part B 报告生成器（自包含 HTML，图片 base64 内嵌）"""
import os, json, csv, base64

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT  = os.path.join(BASE, "results")

R  = json.load(open(os.path.join(OUT, "partB_bridge_results.json"), encoding="utf-8"))
CNT = json.load(open(os.path.join(OUT, "partB_cohort_counts.json"), encoding="utf-8"))
ENR = R["enrichment"]; CON = R["concordance"]; TRI = R["triple"]
BREPL = R["brain_replication_concordance"]

stack = {}
with open(os.path.join(OUT, "partB_evidence_stack.tsv"), encoding="utf-8") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        stack[row["gene"]] = row

def img(name, cap=""):
    p = os.path.join(OUT, name)
    if not os.path.exists(p):
        return ""
    b = base64.b64encode(open(p, "rb").read()).decode()
    c = f"<figcaption>{cap}</figcaption>" if cap else ""
    return f"<figure><img src='data:image/png;base64,{b}'>{c}</figure>"

def nf(x, d=3):
    try:
        return f"{float(x):.{d}f}"
    except Exception:
        return "—"

def psci(x):
    try:
        x = float(x)
    except Exception:
        return "—"
    if x == 0:
        return "<1e-300"
    return f"{x:.2g}" if x >= 1e-3 else f"{x:.2e}"

H = []
A = H.append
A("""<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">
<title>Part B｜SGRG 脑肠共享桥梁检验（UC × MDD）</title>
<style>
:root{--ink:#1a1f27;--muted:#5b6470;--line:#e2e7ee;--bg:#f7f9fc;--card:#fff;
--red:#b42318;--blue:#0b5cab;--amber:#8a5a00;--green:#17643b;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:16px/1.75 -apple-system,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif}
.wrap{max-width:1000px;margin:0 auto;padding:40px 26px 80px}
h1{font-size:27px;line-height:1.35;margin:0 0 6px}
.sub{color:var(--muted);font-size:14px;margin-bottom:26px}
h2{font-size:20px;margin:38px 0 12px;padding-top:14px;border-top:2px solid var(--ink)}
h3{font-size:16px;margin:22px 0 8px;color:var(--blue)}
p{margin:10px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:18px 20px;margin:14px 0}
.warn{background:#fff5f4;border:1px solid #f2c9c4;border-left:5px solid var(--red);border-radius:8px;padding:16px 18px;margin:16px 0}
.ok{background:#f2f9f4;border:1px solid #c8e2d2;border-left:5px solid var(--green);border-radius:8px;padding:16px 18px;margin:16px 0}
.note{background:#fffaef;border:1px solid #f0dcae;border-left:5px solid var(--amber);border-radius:8px;padding:16px 18px;margin:16px 0}
.warn h4,.ok h4,.note h4{margin:0 0 8px;font-size:15px}
table{width:100%;border-collapse:collapse;margin:14px 0;font-size:13.5px;background:#fff}
th,td{border:1px solid var(--line);padding:7px 9px;text-align:left;vertical-align:top}
th{background:#eef2f7;font-weight:600}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
figure{margin:18px 0;text-align:center}
figure img{max-width:100%;border:1px solid var(--line);border-radius:8px;background:#fff}
figcaption{color:var(--muted);font-size:12.5px;margin-top:7px}
.kpi{display:flex;gap:12px;flex-wrap:wrap;margin:16px 0}
.kpi div{flex:1 1 150px;background:#fff;border:1px solid var(--line);border-radius:9px;padding:12px 14px}
.kpi b{display:block;font-size:23px;line-height:1.3}
.kpi span{color:var(--muted);font-size:12.5px}
code{background:#eef2f7;padding:1px 5px;border-radius:4px;font-size:13px}
ul{margin:8px 0 8px 22px}li{margin:5px 0}
.red{color:var(--red)}.green{color:var(--green)}.blue{color:var(--blue)}
footer{margin-top:46px;padding-top:16px;border-top:1px solid var(--line);color:var(--muted);font-size:12.5px}
</style></head><body><div class="wrap">""")

A("<h1>Part B｜SGRG 是否构成 UC × MDD 的脑肠共享桥梁</h1>")
A("<div class='sub'>数据：Gemma 存算差异表达 · MDD 三队列/六脑区 + 外周血（GSE19738 / GSE80655 / GSE53987）· "
  "UC 侧沿用 Part A 自算结果 · 生成于 2026-10-09</div>")

# ------------------------------------------------------------------ 一
A("<h2>一、结论先行（本次是否定性结论，必须按此口径写）</h2>")
A(f"""<div class='warn'><h4>面板级假说被否定：<span class="red">SGRG 基因集既未在 MDD 差异基因中富集，
也与 UC 的位移方向不一致，三重交集为零。</span></h4>
<ul>
<li><b>富集检验不成立</b>：104 个 SGRG 中 {ENR['n']} 个在 GSE19738 可测；其中 {ENR['k']} 个进入 MDD 差异基因
（q&lt;0.05，全基因组 {ENR['K']} 个）。超几何期望 {nf(ENR['expected'],2)} 个，OR={nf(ENR['OR'],2)}，
p={nf(ENR['p'],2)}；2000 次同规模随机面板显示 SGRG 只排在<b>第 {nf(ENR['percentile'],0)} 百分位</b>（中位 {nf(ENR['null_median'],0)}，p95 {nf(ENR['null_p95'],0)}）。
→ <b>与任意一组同规模基因没有区别。</b></li>
<li><b>方向一致性为零</b>：{CON['n']} 个双侧可比基因里同向 {CON['same_sign']} 个（{CON['same_sign']/CON['n']*100:.1f}%），
Spearman ρ={nf(CON['spearman_rho'])}（p={nf(CON['spearman_p'],2)}），二项检验 p={nf(CON['binom_p'],2)}。名义 p&lt;0.05 且双侧可比者
{CON['nominal_both_n']} 个，同向仅 {CON['nominal_both_same']} 个。</li>
<li><b>三重交集为空</b>：SGRG ∩ UC 差异基因（FDR&lt;0.05，{TRI['uc_de_n']} 个）∩ MDD 差异基因（q&lt;0.05，{TRI['mdd_de_panel_n']} 个）
= <b>0 个基因</b>；放宽到名义 p 仍为 0。</li>
</ul>
<p>这不是数据失败，而是这个假说在现有转录组证据下<b>不成立</b>。真正可写的是下面这条基因级线索。</p>
</div>""")

A(f"""<div class='ok'><h4>可保留的增量：面板不成立，但转录层面收敛到少数基因——
OGFOD1（UC 差异 + MDD 临界同向）与 CASC3（UC 临界 + MDD 临界）。</h4>
<p>把桥梁假说从「104 基因面板」收窄到少数候选后，Part B 从否决性分析变成<b>靶点筛选</b>。
但遗传证据一侧需要打折：经 1000G 高覆盖 EUR 面板重做共定位与工具变量 LD 诊断后，
<b>OGFOD1 的 MR 信号被证明是「同一 LD 块内 10 个高度连锁 SNP 被当作独立工具变量」的伪重复</b>
（有效独立工具变量数 = 1.09，lead SNP 单独 Wald p=0.118；详见第九节）。
CASC3 只有 1 个工具变量（脑皮层 Wald ratio），其 Part C 共定位同样不支持（PP.H4=0.108，PP.H3=0.889）。
因此 Part B 的最终结论应表述为：<b>转录层面有候选，遗传因果层面两条线索均未通过</b>。</p></div>""")

A("<div class='kpi'>"
  f"<div><b>{ENR['n']}</b><span>GSE19738 可测 SGRG</span></div>"
  f"<div><b class='red'>{ENR['k']} / {ENR['K']}</b><span>SGRG 进入 MDD 差异基因数 / 全基因组差异基因数</span></div>"
  f"<div><b class='red'>{CON['same_sign']}/{CON['n']}</b><span>UC↔MDD 同向基因数</span></div>"
  f"<div><b class='red'>0</b><span>三重交集基因数</span></div>"
  "</div>")

# ------------------------------------------------------------------ 二
A("<h2>二、数据与方法</h2>")
A("""<p>本环境 GEO 直连不可达（NCBI 全线 000），因此 Part B 的 MDD 侧改用 <b>Gemma 存算差异表达分析</b>
（<code>/rest/v2/datasets/analyses/differential</code>），由 Gemma 统一跑线性模型后通过 REST 提供逐探针统计，
拿到的就是各队列已发布的统计量，无需下载表达矩阵（避免此前 11 KB/s 的下载瓶颈）。
UC 侧沿用 Part A 自算的 GSE38713 面板结果（<code>partA_GSE38713_SGRG_DE.tsv</code>）与 GSE22619 双生子配对验证，
Part C 的 MR 结果作为第三条证据流。</p>""")

A("""<table>
<tr><th>队列</th><th>组织</th><th>样本构成（对照 vs MDD）</th><th>对照口径</th><th>SGRG 可测</th></tr>
<tr><td><b>GSE19738</b></td><td>外周血</td><td>66 对照 / 66 MDD</td>
<td>DiseaseState：对照 vs MDD（单对照，<b>唯一结果集</b>）</td><td>87 / 104（全基因组映射口径 90）</td></tr>
<tr><td><b>GSE80655</b></td><td>死后脑（RNA-seq）</td>
<td>nAcc 22/22 · ACC 24/24 · dlPFC 24/23</td>
<td>三个脑区各自子集分析；同队列另含 SZ/BP 水平，仅取 MDD 对照</td><td>103 / 104</td></tr>
<tr><td><b>GSE53987</b></td><td>死后脑（芯片）</td>
<td>55 对照 / 50 MDD（另 SZ 48 / BP 52）</td>
<td>STR / HPC / dlPFC 三个子集，各含 MDD·SZ·BP 三水平</td><td>98 / 104</td></tr>
<tr><td>GSE38713（UC）</td><td>结肠黏膜</td><td>13 对照 / 15 活动期炎症</td>
<td>DiseaseState：对照 vs 涉及部位活动期 UC</td><td>98（Gemma 存算）/ 97（Part A 自算）</td></tr>
</table>
<p style="font-size:13px;color:#5b6470">「可测」= 该队列存在该基因探针且可解析出对照统计量（by-gene 抓取口径）。
GSE19738 的富集检验分母取 90，来自全基因组翻页文件的探针→基因映射（多符号探针可多登记基因，
比 by-gene 符号查询多 3 个），两口径的差异不影响结论。</p>""")

A("""<div class='note'><h4>口径必须写清的四处</h4>
<ul>
<li><b>子集映射</b>（Gemma 按组织拆分结果集，标签由 <code>subsetFactorValue</code> 核对）：
GSE80655 → rs 493107=nAcc、493108=ACC、493109=dlPFC；GSE53987 → rs 456122=纹状体、456123=海马、456124=dlPFC。</li>
<li><b>MDD 统计量</b>来自 <code>contrasts[]</code> 中 <code>factorValue = major depressive disorder</code> 的那一条
（基线为对照），不是整个线性模型的联合检验。</li>
<li><b>q 值口径</b>：GSE19738 为单对照，探针级 <code>correctedPvalue</code> 即 MDD 的 FDR（全基因组校正，直接可用）；
GSE80655 / GSE53987 的结果集含多疾病水平，探针级 q 非 MDD 专属，故脑区侧只用<b>名义 p</b>，不报 FDR。</li>
<li><b>基因级折叠</b>：同基因多探针取 p 最小的探针为代表（记录于证据叠合表）。</li>
</ul></div>""")

# ------------------------------------------------------------------ 三
A("<h2>三、面板富集：SGRG 与随机基因集没有区别</h2>")
A(img("fig_partB_enrichment.png",
      "左：2000 次同规模（90 基因）随机面板在 GSE19738 MDD 差异基因中的命中数分布；红线为 SGRG 实测 k=5，"
      "虚线为随机中位 4。右：SGRG 面板内 UC 差异基因（FDR&lt;0.05，44 个）与 MDD 差异基因（q&lt;0.05，5 个）无交集。"))
A(f"""<p>富集检验的量：全基因组检验基因 {ENR['N']:,} 个，MDD 差异基因（q≤0.05）{ENR['K']} 个
（对应 867 个探针，Gemma 存算阈值 0.05）。SGRG 面板中可测 {ENR['n']} 个、命中 {ENR['k']} 个，
期望 {nf(ENR['expected'],2)} 个。<b>随机面板零分布的中位就是 4 个、p95 为 8 个</b>，
SGRG 的观测值 5 落在分布中央（第 {nf(ENR['percentile'],0)} 百分位，单侧 p={nf(ENR['p'],2)}）。
这与 Part A Step 2 的结论同构：<b>同规模基因集在 UC 侧也能达到同样甚至更好的判别力</b>。</p>
<p>进入 MDD 差异基因的 5 个 SGRG 是 <code>{'、'.join(ENR['panel_de_genes'])}</code>；
它们在 UC 侧的表现见下表——<b>无一在 UC 达 FDR&lt;0.05</b>，其中 3 个方向相同、2 个相反。</p>""")

_t = "<table><tr><th>基因</th><th class='n'>MDD 血 lfc</th><th class='n'>MDD q</th>" \
     "<th class='n'>UC lfc</th><th class='n'>UC FDR</th><th>UC 配对 (GSE22619)</th></tr>"
for g in ENR["panel_de_genes"]:
    r = stack.get(g, {})
    _t += (f"<tr><td><b>{g}</b></td><td class='n'>{nf(r.get('MDDblood_lfc_GSE19738'))}</td>"
           f"<td class='n'>{psci(r.get('MDDblood_q'))}</td>"
           f"<td class='n'>{nf(r.get('UC_lfc_GSE38713'))}</td>"
           f"<td class='n'>{psci(r.get('UC_FDR'))}</td>"
           f"<td>dz={nf(r.get('UCpaired_dz_GSE22619'))}，p={psci(r.get('UCpaired_p'))}</td></tr>")
_t += "</table>"
A(_t)

# ------------------------------------------------------------------ 四
A("<h2>四、方向一致性：UC 与 MDD 的位移方向互不相关</h2>")
A(img("fig_partB_uc_mdd_scatter.png",
      "SGRG 面板（88 个双侧可比基因）在 UC（GSE38713，活动期 vs 对照）与 MDD（GSE19738 外周血）的 log2FC。"
      "红/蓝为同向（一四象限之外），灰为异向；无任一基因在两侧同时 FDR&lt;0.05。"))
A(f"""<p>若 SGRG 是共享桥梁，应看到左右两病的位移方向整体同向。实测 ρ={nf(CON['spearman_rho'])}（p={nf(CON['spearman_p'],2)}），
同向 {CON['same_sign']}/{CON['n']}——<b>与抛硬币无异</b>。几个有代表性的基因：</p>
<ul>
<li><b>G3BP1</b>：UC <span class="red">+0.220</span>（FDR 1.1e-3，Part A Step1），MDD 血 <span class="green">−0.188</span>（q 0.13）——<b>方向相反</b>，
海马 +0.024、nAcc −0.070，脑区无一致性。</li>
<li><b>EIF4A1</b>：UC +0.533（FDR 1.3e-6，面板内最强上调之一），MDD 血 −0.253（q 0.24）——同样相反。</li>
<li><b>YTHDF3 / ATXN2 / CSDE1 / KPNB1</b>：UC 显著下调，MDD 侧全部不显著（q 0.10–0.26），方向虽多为同号，但效应量接近 0。</li>
</ul>
<p>机制上这并不意外：UC 的差异集中在<b>颗粒清除/翻译抑制模块</b>（Part A Step3 已在独立配对队列复现），
而外周血 MDD 的差异基因以 CIRBP、STAU1 这类<b>应激反应型 RNA 结合蛋白</b>为主（见上表），
两者本就是同一条通路的不同分支。</p>""")

# ------------------------------------------------------------------ 五
A("<h2>五、脑区复现：血量级信号在脑区不成立</h2>")
A(f"""<p>把 5 个 MDD 血差异基因拿到六个死后脑区/核团里看方向（同队列内不同子集，名义 p）：{''}</p>""")
_t = ("<table><tr><th>脑区（队列）</th><th class='n'>可比基因</th><th class='n'>与血同向</th><th>解读</th></tr>")
_def = {
    "GSE53987:hippocampal formation": "唯一全同向，但 n=5 且均为小效应",
    "GSE53987:dorsolateral prefrontal cortex": "4/5 同向，方向性支持",
    "GSE80655:dorsolateral prefrontal cortex": "3/5，接近随机",
    "GSE53987:striatum": "2/5，不支持",
    "GSE80655:anterior cingulate cortex": "2/5，不支持",
    "GSE80655:nucleus accumbens": "2/5，不支持",
}
for k in sorted(BREPL, key=lambda x: -BREPL[x]["same"]):
    v = BREPL[k]
    lab = k.replace("GSE53987:", "").replace("GSE80655:", "") + f"（{k.split(':')[0]}）"
    _t += (f"<tr><td>{lab}</td><td class='n'>{v['n']}</td>"
           f"<td class='n'>{v['same']}</td><td>{_def.get(k,'')}</td></tr>")
_t += "</table>"
A(_t)
A("""<p>结论：<b>外周血与脑组织的 MDD 信号在 SGRG 上没有可复现的一致性</b>。
这既可以是生物学真实（外周免疫-代谢改变与中枢神经炎症本就不同步），
也可能是本组队列的限制（脑区子集每组仅 22–24 例、死后样本、用药史不可控）。写入正文时必须两者并列。</p>""")

# ------------------------------------------------------------------ 六
A("<h2>六、基因级例外：OGFOD1 与 CASC3</h2>")
_t = ("<table><tr><th>基因</th><th class='n'>UC lfc / FDR</th><th class='n'>UC 配对</th>"
      "<th class='n'>MDD 血 lfc / q</th><th class='n'>脑区（最小 p）</th><th>MR（Part C）</th></tr>")
for g in ["OGFOD1", "CASC3", "PQBP1", "G3BP1", "VCP"]:
    r = stack.get(g, {})
    brain = []
    for lab, cl, cp in [("nAcc", "nAcc_lfc", "nAcc_p"), ("ACC", "ACC_lfc", "ACC_p"),
                        ("dlPFC", "dlPFC80655_lfc", "dlPFC80655_p"), ("STR", "STR_lfc", "STR_p"),
                        ("HPC", "HPC_lfc", "HPC_p"), ("dlPFC2", "dlPFC53987_lfc", "dlPFC53987_p")]:
        try:
            p = float(r.get(cp))
            brain.append((p, lab, r.get(cl)))
        except Exception:
            continue
    if brain:
        brain.sort()
        bp, bl, blfc = brain[0]
        btxt = f"{bl} lfc={nf(blfc)}，p={psci(bp)}"
    else:
        btxt = "—"
    mr = r.get("MR_IVW_p")
    mrtxt = f"{r.get('MR_tissue','—')}｜IVW p={psci(mr)}，q={psci(r.get('MR_IVW_q'))}" if mr else "无可构建工具变量"
    _t += (f"<tr><td><b>{g}</b></td>"
           f"<td class='n'>{nf(r.get('UC_lfc_GSE38713'))} / {psci(r.get('UC_FDR'))}</td>"
           f"<td class='n'>dz={nf(r.get('UCpaired_dz_GSE22619'))}，p={psci(r.get('UCpaired_p'))}</td>"
           f"<td class='n'>{nf(r.get('MDDblood_lfc_GSE19738'))} / {psci(r.get('MDDblood_q'))}</td>"
           f"<td class='n'>{btxt}</td><td>{mrtxt}</td></tr>")
_t += "</table>"
A(_t)
A("""<div class='note'><h4>两点必须同时写进正文的解释</h4>
<ul>
<li><b>OGFOD1</b>（2-氧化戊二酸依赖性脯氨酰羟化酶，作用于核糖体蛋白）<b>不再是可写的正面靶点</b>。
Part C 的邻居 eQTL 诊断曾提示其 631 kb 外 CES1P1 的 eQTL p=5e-97（Δlog10P=88.9）；
第九节的敏感性分析给出两个明确结论：
(i) <span class="green">CES1P1 污染被排除</span>——两个 lead SNP 实测 r=−0.06（r<sup>2</sup>=0.004），
10 个工具变量与 CES1P1 全部 215 个 eQTL SNP 的最大 r<sup>2</sup> 仅 0.013；
(ii) <span class="red">但 MR 信号本身不成立</span>——10 个工具变量两两 |r| 中位 0.94（同一 LD 块），
特征值法有效独立工具变量数 <b>1.09</b>，lead SNP 单独 Wald p=0.118（阴性），
而把相关 SNP 当独立的朴素 IVW 才给出 p=1.7e-7；UC GWAS 在整个区域无信号
（最小 p=0.014），coloc.abf 对 OGFOD1 eQTL 的 PP.H4 仅 0.044（PP.H1=0.951）。
→ 正文只能写「名义 MR 阳性但经 LD 敏感性分析后消失」，这条本身是有价值的方法学发现。</li>
<li><b>CASC3</b>（EJC 核心组分）脑皮层仅 1 个工具变量（Wald ratio，p=5.85e-5），
无可用于敏感性分析的 LD 结构；其 Part C 共定位 PP.H4=0.108 / PP.H3=0.889（同样不支持）。
宜作为「待验证候选」，且必须与 OGFOD1 用同一套 LD 敏感性流程复核后才能进正文结论。</li>
</ul></div>""")

A("<h2>七、证据叠合全图</h2>")
A(img("fig_partB_evidence_heatmap.png",
      "24 个在任一证据流达标的 SGRG 基因 × 10 条证据流。色标为 sign(lfc)×(−log10 p)（红=上调，蓝=下调）；"
      "星号=该流 FDR/q&lt;0.05（配对验证为 Wilcoxon p&lt;0.05）；「·」= 该队列未测到该基因。"))
A("""<p>读图方式：<b>纵向（同一基因跨证据流）几乎看不到整列贯通</b>——这正是「面板不成立」的直观表达；
横向（同一证据流）则能看到 Part A 的模块结构（VCP、RBPMS、YTHDF3、EIF2S1 在 UC 侧成组显著），
以及 MR 列的零星命中（OGFOD1、CASC3、HABP4 等）。</p>""")

# ------------------------------------------------------------------ 八
A("<h2>八、这对稿件意味着什么（写作口径）</h2>")
A("""<div class='card'>
<h3>建议把 Part B 从「创新框架」改写为「假说检验 + 阴性结果 + 方法学警示」</h3>
<ul>
<li><b>不建议写</b>：「SGRG 是 UC 与抑郁的共享分子桥梁」——本次富集与方向一致性都否定该表述；
更不能写「OGFOD1 是遗传因果靶点」——见第九节，该 MR 信号是伪重复。</li>
<li><b>建议写</b>：「我们检验了 SGRG 作为脑肠共享桥梁的假说：面板层面不成立（附富集与方向一致性检验），
转录层面仅收敛到个别候选（OGFOD1、CASC3），且二者的遗传因果证据经 LD 敏感性分析后均不成立」。
阴性结果 + <b>对自身既往阳性结论的主动纠错</b>（距离剪枝 → 实测 LD）构成方法学卖点，
比硬报一个经不起 LD 复核的 MR 阳性安全得多。</li>
<li><b>新增一个方法学要点（可直接进 Discussion）</b>：以「距离 &gt;10 kb」代替 LD clumping 时，
cis 窗口内高度连锁的 SNP 会被当成独立工具变量，使 IVW 的 p 值虚高数个数量级
（本例：10 个 SNP、|r| 中位 0.94 → 朴素 IVW p=1.7e-7，而有效独立数 1.09、lead Wald p=0.118）。
若稿件保留 MR 部分，必须补 LD 矩阵与有效工具变量数报告。</li>
<li><b>Fig 10（证据叠合）保持基因版</b>，但需加一行「LD 敏感性」色块标明 MR 证据降级；
新增本报告的 Fig B1（<code>fig_partB_ogfod1_locus.png</code>）作为该位点的三轨道区域图进补充材料。</li>
<li>与已发表工作的差异化仍然成立：既有 IBD-抑郁共享基因研究（<i>J Transl Med</i> 2024）锚定分泌蛋白基因集，
本文是首个用 SGRG 锚定并<b>报告阴性面板结果 + MR 伪重复纠错</b>的工作。</li>
</ul></div>""")

# ------------------------------------------------------------------ 九（新增）
A("<h2>九、OGFOD1 位点共定位敏感性分析（新增 · 纠正第六节结论）</h2>")
A("""<div class='note'><h4>为什么要做这一步</h4>
<p>第六节把 OGFOD1 列为「集齐三流证据」的候选，但存在两个未决问题：
(1) Part C 的共定位是在 GTEx「仅显著 eQTL」的截断 SNP 集上做的，PP.H4 不可靠；
(2) Part C 的 MR 用「距离 &gt;10 kb」代替 LD clumping，未经实测 LD 验证工具变量独立性。
本节用 <b>1kGP 高覆盖 GRCh38 定相面板（EUR n=633）</b>补做两件事：
以实测 LD 复核工具变量独立性，并分别在 OGFOD1 与其 631 kb 外的 CES1P1 两个窗口上重做共定位。</p></div>""")

_a = {}
with open(os.path.join(OUT, "ogfod1_coloc_abf_samesnpset.csv"), encoding="utf-8") as f:
    for r in csv.DictReader(f):
        _a[r["region"]] = r
_d = {}
with open(os.path.join(OUT, "ogfod1_coloc_diagnostics.csv"), encoding="utf-8") as f:
    for r in csv.DictReader(f):
        _d[r["region"]] = r
_ins = list(csv.DictReader(open(os.path.join(OUT, "ogfod1_instrument_ld.tsv"), encoding="utf-8"), delimiter="\t"))
_align = json.load(open(os.path.join(OUT, "ogfod1_coloc_align_report.json"), encoding="utf-8"))

A(f"""<h3>9.1 数据与对齐</h3>
<p>eQTL：GTEx v8 Whole_Blood（OGFOD1 {_align['genes']['OGFOD1']['n_eqtl']} 个显著 cis-eQTL SNP、
CES1P1 {_align['genes']['CES1P1']['n_eqtl']} 个）。UC GWAS：de Lange 2017（GCST004133，12,366 病例 / 33,609 对照）。
LD 面板：1kGP 20220422 高覆盖定相集 chr16:55.3–57.6 Mb（GRCh38），EUR n=633，按 GTEx ALT 等位基因定向
（REF/ALT 互换者剂量取补）。由于该 GWAS 的 harmonised 文件<b>完全没有 EAF 列</b>，
回文 SNP（A/T、C/G）无法定向，一律剔除（OGFOD1 {_align['genes']['OGFOD1']['drops'].get('palindromic_no_gwas_eaf', 0)} 个、
CES1P1 {_align['genes']['CES1P1']['drops'].get('palindromic_no_gwas_eaf', 0)} 个），
其余 SNP 的定向由等位基因标签唯一确定（无歧义）。
最终分析集：<b>OGFOD1 {_a['OGFOD1']['n_snp']} 个 SNP</b>、<b>CES1P1 {_a['CES1P1']['n_snp']} 个 SNP</b>。</p>""")

A("<h3>9.2 工具变量 LD 诊断：MR 信号是伪重复</h3>")
_t = ("<table><tr><th>工具变量</th><th class='n'>eQTL p</th><th class='n'>r<sup>2</sup> vs OGFOD1 lead</th>"
      "<th class='n'>r<sup>2</sup> vs CES1P1 lead</th></tr>")
for r in _ins:
    _t += (f"<tr><td><code>{r['snp']}</code></td><td class='n'>{psci(r['eqtl_p'])}</td>"
           f"<td class='n'>{r.get('r2_vs_OGFOD1_lead') or '—'}</td>"
           f"<td class='n'>{r.get('r2_vs_CES1P1_lead') or '—'}</td></tr>")
_t += "</table>"
A(_t)
A(f"""<div class='bad'><h4><span class="red">结论一：CES1P1 污染被排除，但 MR 信号本身是同一 LD 块的伪重复</span></h4>
<ul>
<li><b>CES1P1 污染排除</b>：OGFOD1 lead（rs4784670）与 CES1P1 lead（rs12149359）实测
r={_align.get('lead_lead_r')}（r<sup>2</sup>={_align.get('lead_lead_r2')}）；
10 个工具变量与 CES1P1 窗口内<b>全部 215 个 eQTL SNP</b> 的最大 r<sup>2</sup> 仅 0.013。
两者不在同一 LD 块，「631 kb 外的强 eQTL 污染 cis 工具变量」这一担忧可以正式排除。</li>
<li><b>但工具变量彼此几乎完全连锁</b>：10 个工具变量两两 |r| 中位 <b>0.94</b>、最小 0.91、最大 1.00，
特征值法有效独立工具变量数 = <b>1.09</b>。即把同一个信号重复计数了约 10 次。</li>
<li><b>删掉伪重复后证据消失</b>：lead SNP（rs4784670）单独 Wald ratio
b={nf(0.1725)}、<b>p=0.118</b>（阴性）；而 Part C 的朴素 IVW 给出
b=0.198、p={psci(1.72e-7)}。两者点估计接近、p 值相差 <b>6 个数量级</b>——
差异全部来自把相关 SNP 当作独立工具变量。</li>
<li>同类问题波及 Part C 报告中 OGFOD1 的 MR-Egger（p={psci(8.25e-20)}）、
加权中位数（p={psci(1.22e-14)}）与 FDR q={psci(2.15e-6)}：这些都建立在同一组伪独立工具变量上，
应一并降级。</li>
</ul></div>""")

A("<h3>9.3 共定位：两个窗口内 UC GWAS 都没有信号</h3>")
_t = ("<table><tr><th>窗口</th><th class='n'>SNP</th><th class='n'>GWAS 最小 p</th>"
      "<th class='n'>PP.H1</th><th class='n'>PP.H3</th><th class='n'>PP.H4</th><th>susie 可信集</th></tr>")
for g, lab in [("OGFOD1", "OGFOD1（56.3–56.5 Mb）"), ("CES1P1", "CES1P1（55.7–55.9 Mb）")]:
    a, d = _a[g], _d[g]
    gmin = "0.014（rs3790110）" if g == "OGFOD1" else "0.044（rs71387105）"
    _t += (f"<tr><td><b>{lab}</b></td><td class='n'>{a['n_snp']}</td><td class='n'>{gmin}</td>"
           f"<td class='n'>{nf(a['PP.H1'],4)}</td><td class='n'>{nf(a['PP.H3'],4)}</td>"
           f"<td class='n'>{nf(a['PP.H4'],4)}</td><td>eQTL {d['region'] and ''}"
           f"1 / 3 个；UC <b>0 个</b> → coloc.susie 无 CS 对可取</td></tr>")
_t += "</table>"
A(_t)
A("""<div class='note'><h4>结论二的读法</h4>
<p>两个窗口的 coloc.abf 都以 <b>PP.H1 主导</b>（0.951 / 0.901）——即该区域的关联信号
几乎完全由 eQTL 一侧贡献，UC GWAS 在此区域没有可信信号（最小 p 分别为 0.014、0.044，
均远未达基因组显著）。因此 <b>PP.H4 低（0.044 / 0.096）不是因为「两个信号位置不同」（PP.H3 仅 0.005/0.003），
而是因为 GWAS 侧根本没有信号</b>。这与「MR 伪重复」互相印证：工具变量的 GWAS 效应本身就是噪声，
被重复计数后放大成 p=1.7e-7。</p>
<p><b>三重校验（按 coloc-susie-external-ld 规范）</b>：
① 同一 SNP 集上重跑 coloc.abf 作为可比基准（本表即为该结果）；
② 零噪声合成阳性对照 z<sub>peak</sub>=8 × R<sub>,peak</sub> → 两窗口均返回<b>恰好 1 个可信集且包含真峰</b>，
证明精细定位代码路径正确；
③ L=1 与 estimate_residual_variance=TRUE 两种设置下 UC 侧仍为 0 个可信集，结论不随设置变化。
susie 在 eQTL 侧检出 1（OGFOD1）/ 3（CES1P1）个可信集，与 CES1P1 窗口存在多个独立 eQTL 信号一致。</p>
<p><b>面板质量</b>：OGFOD1 窗口 126 SNP 中 |r|&gt;0.98 的 SNP 对 1,568 对（中位 |r| 0.84，
高相关对最小 MAF 0.143）；CES1P1 窗口 152 SNP 中 373 对（中位 |r| 0.13）。
OGFOD1 窗口本身就是一个大 LD 块——这既是工具变量退化为 1 个的原因，也提示该区域任何 MR 都必须带 LD 矩阵。</p></div>""")

A(img("fig_partB_ogfod1_locus.png",
      "Fig B1｜chr16:55.6–56.6 Mb 三轨道区域图（GRCh38）。上：UC GWAS（点按与 OGFOD1 lead 的 r² 着色；"
      "全区域无一 SNP 达基因组显著，最小 p=0.014）；中：OGFOD1 Whole_Blood eQTL；下：CES1P1 Whole_Blood eQTL。"
      "两个 lead 之间实测 r²=0.004（不同 LD 块）；图中框注给出共定位后验、有效独立工具变量数与 Wald/IVW 对照。"))

# ------------------------------------------------------------------ 十
A("<h2>十、局限（建议整段进 Limitations）</h2>")
A("""<ol>
<li><b>MDD 主队列替换</b>：原计划的 GSE98793 不在 Gemma、且本环境 GEO 不可达，改用 GSE19738（n=132）作主队列；
二者平台与样本构成不同，结论外推需谨慎。</li>
<li><b>用药混杂</b>：GSE19738 含 Treatment 因子，MDD 组用药史不可控，外周血信号可能部分是药物效应。</li>
<li><b>脑区样本量小</b>：GSE80655 每区 22–24 例/组，GSE53987 每区 50 MDD/55 对照但跨三区共享；
脑区侧仅报名义 p，未做跨区 FDR。</li>
<li><b>MDD 侧 q 值不可跨队列比较</b>：仅 GSE19738 有单对照的探针级 FDR；其余队列的探针级 q 含多疾病水平，不可直接用。</li>
<li><b>共定位的 SNP 集仍受限</b>：GTEx REST 只提供已显著 eQTL，两个窗口的分析集为
126 / 152 个 SNP（非全区域变异），PP.H4 的解释需在此前提下进行；第九节已用同 SNP 集的
coloc.abf + 三重校验（含零噪声阳性对照）保证内部一致性。</li>
<li><b>GWAS 缺 EAF</b>：GCST004133 harmonised 文件无等位基因频率列，回文 SNP（OGFOD1 17 个、
CES1P1 23 个）被剔除；若换用带 EAF 的 GWAS 版本可回收这部分位点。</li>
<li><b>LD 面板与 GWAS 人群不完全匹配</b>：面板为 1kGP EUR（n=633），GWAS 为欧洲混合人群；
且候选区域是被 eQTL 显著性命中的区域（selection），有轻度选择偏倚。</li>
<li><b>面板定义敏感性</b>：SGRG 采用 GO:0010494 ∪ GO:0035061 的 104 个 HGNC 官方符号；
换用 MsGP 等更宽定义可能改变富集结论。</li>
<li><b>组织不可配对</b>：UC 为结肠活检、MDD 为外周血/死后脑，来自不同个体，无法做同体配对设计。</li>
<li><b>外周血 ≠ 中枢</b>：MDD 外周血信号能否反映中枢病理，本身存疑；本次脑区复现失败即一例证。</li>
</ol>""")

# ------------------------------------------------------------------ 十
A("<h2>十一、复现清单</h2>")
A("""<table>
<tr><th>类别</th><th>文件</th><th>说明</th></tr>
<tr><td>数据抓取</td><td><code>scripts/partB_fetch_mdd_de.py</code></td><td>SGRG × 4 数据集 by-gene 抓取（含缓存与重试）</td></tr>
<tr><td>数据抓取</td><td><code>scripts/partB_genomewide_19738.py</code></td><td>rs 540049 全翻页（23,624 探针）</td></tr>
<tr><td>样本构成</td><td><code>scripts/partB_cohort_counts.py</code></td><td>各队列/子集疾病×脑区计数</td></tr>
<tr><td>分析</td><td><code>scripts/partB_bridge_analysis.py</code></td><td>三重交集、富集、方向一致性、脑区复现、证据叠合</td></tr>
<tr><td>图件</td><td><code>scripts/partB_figures.py</code></td><td>三张主图</td></tr>
<tr><td>结果</td><td><code>results/partB_sgrg_mdd_de.tsv</code></td><td>1,071 条 基因×队列×探针×对照 记录</td></tr>
<tr><td>结果</td><td><code>results/partB_gse19738_genomewide.tsv</code></td><td>全基因组 23,624 探针（universe + DEG）</td></tr>
<tr><td>结果</td><td><code>results/partB_evidence_stack.tsv</code></td><td>证据叠合表（104 基因 × 21 列）</td></tr>
<tr><td>结果</td><td><code>results/partB_bridge_results.json</code></td><td>全部核心统计量</td></tr>
<tr><td>敏感性-LD</td><td><code>scripts/ld_fetch_1kgp38.py</code></td><td>1kGP 高覆盖 GRCh38 EUR 面板（chr16:55.3–57.6 Mb，10,217 变异 × 633 样本）</td></tr>
<tr><td>敏感性-对齐</td><td><code>scripts/partB_ogfod1_coloc_align.py</code></td><td>eQTL×GWAS×面板三方对齐、等位基因定向、LD 矩阵、工具变量×lead LD</td></tr>
<tr><td>敏感性-共定位</td><td><code>scripts/partB_ogfod1_coloc_susie.R</code></td><td>coloc.abf 基准 + coloc.susie + 三重校验（L=1、零噪声阳性对照、面板病态统计）</td></tr>
<tr><td>敏感性-图件</td><td><code>scripts/partB_ogfod1_figure.py</code></td><td>三轨道区域图 Fig B1</td></tr>
<tr><td>结果</td><td><code>results/ogfod1_coloc_abf_samesnpset.csv</code> 等</td><td>共定位后验、诊断、工具变量 LD 表、对齐报告</td></tr>
</table>
<p>关键命令：<code>python scripts/partB_fetch_mdd_de.py</code> → <code>python scripts/partB_genomewide_19738.py</code>
→ <code>python scripts/partB_bridge_analysis.py</code> → <code>python scripts/partB_figures.py</code>
→ <code>python scripts/ld_fetch_1kgp38.py</code> → <code>python scripts/partB_ogfod1_coloc_align.py</code>
→ <code>Rscript scripts/partB_ogfod1_coloc_susie.R</code> → <code>python scripts/partB_ogfod1_figure.py</code>
→ <code>python scripts/make_partB_report.py</code>
（R 需 <code>coloc</code>/<code>susieR</code>，且需在 ASCII 路径下运行——Windows 版 R 无法读取含中文的路径）。</p>""")

A("""<footer>Part B · SGRG 脑肠共享桥梁检验｜MDD 数据来自 Gemma REST v2 存算差异表达分析；
UC 数据来自 Part A 自算（GSE38713 + GSE22619）；遗传因果证据来自 Part C；
第九节共定位与 LD 诊断基于 1kGP 高覆盖 GRCh38 EUR 面板（n=633）。
本报告所有数字可由上表脚本与结果文件逐一复核。</footer>""")

A("</div></body></html>")

path = os.path.join(OUT, "PartB_SGRG_脑肠共享桥梁报告.html")
open(path, "w", encoding="utf-8").write("".join(H))
print("[write]", path, os.path.getsize(path) // 1024, "KB")
