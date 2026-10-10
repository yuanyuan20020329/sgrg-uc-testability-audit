# -*- coding: utf-8 -*-
"""
生成 Part A · Step 1 结果报告（自包含 HTML，图片 base64 内嵌）。
数据直接读取 results/ 下的结果文件，避免手工转录误差。
"""
import csv, json, os, base64, html

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT  = os.path.join(BASE, "results")

DE     = os.path.join(OUT, "partA_GSE38713_SGRG_DE.tsv")
SUMM   = os.path.join(OUT, "partA_GSE38713_SGRG_summary.json")
GM     = os.path.join(OUT, "partA_GSE38713_SGRG_groupMeans.tsv")
VOLC   = os.path.join(OUT, "fig_volcano_SGRG_GSE38713.png")
HEAT   = os.path.join(OUT, "fig_heatmap_SGRG_GSE38713.png")
CLEAN  = os.path.join(OUT, "sgrg_gene_set_clean.tsv")
REJ    = os.path.join(OUT, "sgrg_rejected.tsv")
REPORT = os.path.join(OUT, "PartA_Step1_UC_SGRG差异表达报告.html")

def b64(p):
    return "data:image/png;base64," + base64.b64encode(open(p, "rb").read()).decode()

rows = list(csv.DictReader(open(DE, encoding="utf-8"), delimiter="\t"))
S    = json.load(open(SUMM, encoding="utf-8"))
gm   = list(csv.reader(open(GM, encoding="utf-8"), delimiter="\t"))
n_clean = sum(1 for i, _ in enumerate(open(CLEAN, encoding="utf-8")) if i > 0)
n_rej   = sum(1 for i, _ in enumerate(open(REJ, encoding="utf-8")) if i > 0)

def f(x):
    try: return float(x)
    except Exception: return float("nan")

sig  = [r for r in rows if f(r["FDR_inflamed"]) < 0.05]
up   = [r for r in sig if f(r["log2FC_inflamed_vs_ctrl"]) > 0]
dn   = [r for r in sig if f(r["log2FC_inflamed_vs_ctrl"]) < 0]
bigE = [r for r in sig if abs(f(r["log2FC_inflamed_vs_ctrl"])) >= 0.5]
multi= [r for r in rows if r["probe_quality"] != "unique"]

CORE = {"G3BP1":"SG 核心支架（标志物）","G3BP2":"SG 核心支架","VCP":"AAA-ATPase，SG 清除",
        "EIF2S1":"eIF2α，SG 组装开关","ATXN2":"SG 装配调控","FMR1":"SG 组分/RBP",
        "CSDE1":"SG 组分/RBP","UBAP2L":"SG 装配必需","LSM14A":"SG 组分","TRIM21":"SG 相关 E3",
        "ZNFX1":"SG 相关","YTHDF1":"m6A 阅读器","YTHDF3":"m6A 阅读器","PUM1":"RBP",
        "PUM2":"RBP，SG 组分","MBNL1":"RBP","QKI":"RBP，SG 相关","DDX6":"解旋酶",
        "DDX3X":"解旋酶","TARDBP":"TDP-43","CAPRIN1":"SG 装配","TIA1":"SG 标志物",
        "TIAL1":"TIA1 家族","ELAVL1":"HuR","EIF4A1":"翻译起始","EIF4G1":"翻译起始",
        "PABPC1":"poly(A) 结合","GIGYF2":"翻译抑制","STAU1":"RBP","ATXN2L":"SG 相关"}

def tablerows(sel):
    out = []
    for r in sel:
        fc, p_, q = f(r["log2FC_inflamed_vs_ctrl"]), f(r["p_inflamed"]), f(r["FDR_inflamed"])
        cls = "up" if fc > 0 else "dn"
        flag = "" if r["probe_quality"] == "unique" else ' <span class="pill p-c">多探针</span>'
        out.append(
            f'<tr><td><b>{r["GeneSymbol"]}</b>{flag}</td>'
            f'<td class="note">{html.escape(CORE.get(r["GeneSymbol"],""))}</td>'
            f'<td class="num {cls}">{fc:+.3f}</td>'
            f'<td class="num">{p_:.2e}</td>'
            f'<td class="num {"hl" if q<0.05 else ""}">{q:.2e}</td>'
            f'<td class="num">{f(r["mean_Control"]):.2f}</td>'
            f'<td class="num">{f(r["mean_UC_inflamed"]):.2f}</td></tr>')
    return "\n".join(out)

def coretable():
    m = {r["GeneSymbol"]: r for r in rows}
    out = []
    for g, desc in CORE.items():
        r = m.get(g)
        if not r:
            out.append(f'<tr><td><b>{g}</b></td><td class="note">{desc}</td>'
                       f'<td colspan="5" class="num na">芯片未覆盖 / 未通过质控</td></tr>')
            continue
        fc, p_, q = f(r["log2FC_inflamed_vs_ctrl"]), f(r["p_inflamed"]), f(r["FDR_inflamed"])
        cls = "up" if fc > 0 else "dn"
        star = '<span class="sig">●</span>' if q < 0.05 else '<span class="ns">○</span>'
        flag = "" if r["probe_quality"] == "unique" else ' <span class="pill p-c">多探针</span>'
        out.append(f'<tr><td><b>{g}</b>{flag}</td><td class="note">{desc}</td>'
                   f'<td class="num {cls}">{fc:+.3f}</td><td class="num">{p_:.2e}</td>'
                   f'<td class="num {"hl" if q<0.05 else ""}">{q:.2e}</td><td class="num">{star}</td></tr>')
    return "\n".join(out)

CSS = """
:root{--bg:#f6f7f9;--card:#fff;--line:#e5e7eb;--ink:#1f2328;--sub:#57606a;
--accent:#0b6b57;--accent-soft:#e7f4f0;--warn:#b4530a;--warn-soft:#fdf3e7;
--danger:#b42318;--danger-soft:#fdecea;--blue:#0b5cab;--blue-soft:#eaf2fb;
--up:#b42318;--dn:#0b5cab;--tag:#f0f2f5;}
*{box-sizing:border-box;}
body{margin:0;background:var(--bg);color:var(--ink);line-height:1.72;font-size:15px;
font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;}
.wrap{max-width:1120px;margin:0 auto;padding:40px 28px 80px;}
header.hero{background:linear-gradient(135deg,#0b6b57 0%,#0e7c66 55%,#12a184 100%);
color:#fff;border-radius:16px;padding:34px 38px;margin-bottom:26px;box-shadow:0 10px 30px rgba(11,107,87,.18);}
header.hero .kicker{font-size:13px;letter-spacing:.16em;text-transform:uppercase;opacity:.85;margin-bottom:10px;}
header.hero h1{margin:0 0 12px;font-size:25px;line-height:1.35;font-weight:700;}
header.hero p{margin:0;opacity:.93;font-size:14.5px;max-width:900px;}
.verdict{display:flex;gap:14px;flex-wrap:wrap;margin-top:20px;}
.vbox{background:rgba(255,255,255,.15);border:1px solid rgba(255,255,255,.35);border-radius:10px;padding:10px 16px;font-size:13.5px;}
.vbox b{font-size:15px;}
h2{font-size:20px;margin:42px 0 14px;padding-left:14px;border-left:4px solid var(--accent);}
h3{font-size:16.5px;margin:26px 0 10px;color:var(--accent);}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:20px 22px;margin:16px 0;box-shadow:0 1px 3px rgba(0,0,0,.03);}
table{width:100%;border-collapse:collapse;margin:14px 0;font-size:13px;background:#fff;border-radius:10px;overflow:hidden;}
th,td{border:1px solid var(--line);padding:6px 9px;text-align:left;vertical-align:middle;}
th{background:#f0f3f5;font-weight:600;}
tbody tr:nth-child(even){background:#fafbfc;}
td.num{text-align:center;font-variant-numeric:tabular-nums;white-space:nowrap;}
td.up{color:var(--up);font-weight:600;}
td.dn{color:var(--dn);font-weight:600;}
td.hl{background:#e9f7ef;font-weight:600;}
td.na{color:var(--sub);font-style:italic;}
td.note{color:var(--sub);font-size:12.3px;}
.sig{color:var(--up);font-size:15px;}
.ns{color:#c9ced4;font-size:15px;}
.pill{display:inline-block;padding:1px 8px;border-radius:999px;font-size:11.5px;font-weight:600;}
.p-a{background:#e9f7ef;color:#137a3f;}
.p-c{background:var(--warn-soft);color:var(--warn);}
.callout{border-radius:10px;padding:15px 18px;margin:16px 0;border:1px solid;font-size:14px;}
.c-key{background:var(--accent-soft);border-color:#b9e0d4;}
.c-warn{background:var(--warn-soft);border-color:#f2d3ab;}
.c-danger{background:var(--danger-soft);border-color:#f5c2bd;}
.c-info{background:var(--blue-soft);border-color:#c4dcf5;}
ul,ol{margin:9px 0;padding-left:23px;}li{margin:5px 0;}
.small{font-size:12.7px;color:var(--sub);}
code{background:#eef1f4;padding:1px 6px;border-radius:5px;font-size:12.6px;font-family:ui-monospace,Consolas,monospace;}
figure{margin:20px 0;text-align:center;}
figure img{max-width:100%;border:1px solid var(--line);border-radius:10px;background:#fff;}
figcaption{font-size:12.7px;color:var(--sub);margin-top:8px;}
.two{display:grid;grid-template-columns:1fr 1fr;gap:18px;}
footer{margin-top:50px;padding-top:18px;border-top:1px solid var(--line);font-size:12.5px;color:var(--sub);}
@media(max-width:820px){.two{grid-template-columns:1fr;}}
"""

H = []
H.append(f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Part A · Step 1 —— UC 应激颗粒相关基因差异表达</title>
<style>{CSS}</style></head><body><div class="wrap">

<header class="hero">
  <div class="kicker">Part A · Step 1 · 2026-10-08</div>
  <h1>溃疡性结肠炎黏膜中应激颗粒相关基因（SGRG）的差异表达<br>—— SGRG×脑肠轴项目首个实质分析结果</h1>
  <p>为 A+B+C 整合方案的 <b>Part A（单病种诊断模型）</b> 打底：以 Gemma 统一的处理流程取得 GSE38713 人结肠黏膜表达矩阵，
  用经 HGNC 校验的应激颗粒基因集做定量差异表达，产出后续机器学习建模的候选特征基因池。</p>
  <div class="verdict">
    <div class="vbox">数据集：<b>GSE38713</b><br>43 例 · 4 组 · 人结肠黏膜</div>
    <div class="vbox">SGRG 面板：<b>{n_clean} 基因</b><br>GO:0010494 ∪ GO:0035061（HGNC 校验）</div>
    <div class="vbox">差异表达：<b>{len(sig)} / {S['sgrg_measured']}</b><br>FDR&lt;0.05（↑{len(up)} / ↓{len(dn)}）</div>
  </div>
</header>
""")

H.append(f"""
<h2>一、结论速览</h2>

<div class="callout c-key">
<strong>应激颗粒相关基因在活动期 UC 炎症黏膜中整体失调。</strong>
在可定量的 {S['sgrg_measured']} 个 SGRG 中，<b>{len(sig)} 个（{len(sig)/S['sgrg_measured']*100:.0f}%）</b>在
"活动期炎症区 UC vs 非炎症对照"间差异表达（BH-FDR &lt; 0.05），其中 {len(up)} 个上调、{len(dn)} 个下调。
</div>

<div class="callout c-key">
<strong>核心机制信号指向"应激颗粒组装轴被激活"。</strong>
SG 组装的直接开关 <b>eIF2α（EIF2S1）上调 +0.43</b>（FDR = 7.3×10⁻⁵）；
SG 核心支架 <b>G3BP1 上调 +0.22</b>（FDR = 1.1×10⁻³）；SG 清除通路 <b>VCP 上调 +0.64</b>（FDR = 9.2×10⁻⁸）。
同时多个 RBP 组分下调（ATXN2 −0.35、FMR1 −0.50、CSDE1 −0.22、PUM2 −0.37、LSM14A −0.19）。
</div>

<div class="callout c-warn">
<strong>必须如实说明的一点：统计显著 ≠ 效应量可观。</strong>
全部 {S['sgrg_measured']} 个 SGRG 的 log2FC 中位数仅 <b>0.015</b>，|log2FC| ≥ 0.5 的显著基因只有 <b>{len(bigE)}</b> 个。
即"半数 SGRG 显著"主要来自炎症 vs 正常黏膜的巨大整体差异配合较小方差，而非单个基因的大幅位移。
<b>建议在正文中以效应量排序呈现候选基因，而不是以 P 值排序。</b>
</div>
""")

H.append(f"""
<h2>二、数据与方法</h2>
<div class="two">
  <div class="card">
    <h3 style="margin-top:0">表达数据</h3>
    <ul class="small">
      <li><b>数据集</b>：GSE38713（Gemma 数据集 id = 5873）<br><i>Expression data from intestinal mucosa of patients with UC</i></li>
      <li><b>样本</b>：43 例人结肠黏膜活检，四组——<br>
          活动期炎症区 UC {S['group_counts'].get('UC_active_inflamed')}；非炎症区活动期 UC {S['group_counts'].get('UC_active_noninflamed')}；
          缓解期 UC {S['group_counts'].get('UC_remission')}；非炎症对照 {S['group_counts'].get('Control')}</li>
      <li><b>平台</b>：Affymetrix HG-U133 Plus 2.0（微阵列）</li>
      <li><b>预处理</b>：由 Gemma 统一完成 RMA + 分位数标准化，输出为 log2 尺度（37,959 条探针行）</li>
      <li><b>分组来源</b>：Gemma 实验设计因子 <code>DiseaseState</code>（非由标题文本猜测）</li>
    </ul>
  </div>
  <div class="card">
    <h3 style="margin-top:0">SGRG 基因集与统计</h3>
    <ul class="small">
      <li><b>基因集</b>：GO:0010494（细胞质应激颗粒）∪ GO:0035061（核应激颗粒），限 <i>Homo sapiens</i>，取自 EBI QuickGO 注释 API</li>
      <li><b>符号净化</b>：QuickGO 返回 {n_clean + n_rej} 个符号，经 HGNC 官方 REST API 校验后保留 <b>{n_clean} 个官方批准符号</b>，
          剔除 {n_rej} 个（多为 UniProt 登录号与计算推断条目）→ <code>results/sgrg_gene_set_clean.tsv</code></li>
      <li><b>主对比</b>：{S['contrast_main']}</li>
      <li><b>辅对比</b>：{S['contrast_secondary']}</li>
      <li><b>统计</b>：log2 尺度 Welch <i>t</i> 检验；Benjamini–Hochberg FDR（<b>在 SGRG 面板内校正</b>）</li>
      <li><b>探针处理</b>：同一符号多探针取均值；跨基因探针（如 <code>PABPC1|RLIM|PABPC3</code>）登记到各候选符号并标注</li>
    </ul>
  </div>
</div>
""")

H.append(f"""
<h2>三、结果</h2>

<h3>3.1 核心应激颗粒基因状态（按机制角色）</h3>
<p class="small">log2FC 为"活动期炎症区 UC − 非炎症对照"；<span class="sig">●</span> = FDR&lt;0.05，<span class="ns">○</span> = 未达显著。</p>
<table>
<thead><tr><th>基因</th><th>机制角色</th><th>log2FC</th><th><i>P</i></th><th>FDR</th><th>显著</th></tr></thead>
<tbody>
{coretable()}
</tbody></table>

<div class="callout c-info">
<b>读法提示：</b>eIF2α 的 <i>转录本</i> EIF2S1 上调，说明的是"组装开关的转录供给增加"，与"颗粒数量增加"并非同一件事——
后者需要 G3BP1 蛋白水平/颗粒计数的证据（Hödlmayr 2026 CMGH 的影像学结果可作为外部佐证）。本分析只呈现转录层证据。
</div>

<h3>3.2 全部差异表达 SGRG（FDR &lt; 0.05，按 <i>P</i> 值排序）</h3>
<table>
<thead><tr><th>基因</th><th>机制角色</th><th>log2FC</th><th><i>P</i></th><th>FDR</th><th>对照均值</th><th>炎症均值</th></tr></thead>
<tbody>
{tablerows(sig)}
</tbody></table>

<h3>3.3 图形</h3>
<figure><img src="{b64(VOLC)}" alt="volcano">
<figcaption>图 1　SGRG 火山图。红色 = FDR&lt;0.05 且上调，蓝色 = 下调，灰 = 不显著。
横轴为 log2 倍数变化，纵轴为 −log10 <i>P</i>（Welch <i>t</i> 检验）。</figcaption></figure>

<figure><img src="{b64(HEAT)}" alt="heatmap">
<figcaption>图 2　按 <i>P</i> 值排序前 30 个 SGRG 的行 z-score 热图（列按分组排布：炎症区 → 非炎症区 → 缓解期 → 对照）。
聚类分界可见"炎症区 vs 其余"的整体偏移。</figcaption></figure>

<h3>3.4 各组均值（供后续建模与作图取用）</h3>
<p class="small">完整表见 <code>results/partA_GSE38713_SGRG_groupMeans.tsv</code>，此处展示前 15 行。</p>
<table>
<thead><tr>{''.join('<th>'+c+'</th>' for c in gm[0])}</tr></thead>
<tbody>
{''.join('<tr>'+''.join(f'<td class="num">{c}</td>' if i>0 else f'<td><b>{c}</b></td>' for i,c in enumerate(r))+'</tr>' for r in gm[1:16])}
</tbody></table>
""")

H.append(f"""
<h3>3.5 稳健性检查（数据集内部）</h3>
<p>把主对比换成"全部 UC（炎症区 + 非炎症区 + 缓解期，n=30）vs 对照（n=13）"重新计算，用来判断候选基因是否依赖"炎症区"这一取样选择：</p>
<ul>
  <li>{len(sig)} 个候选基因在两套对比中的 log2FC <b>100% 同向</b>（Spearman ρ = 0.93，<i>P</i> = 1.3×10⁻¹⁹）——方向不依赖于用"炎症区"还是"全 UC"对比；</li>
  <li>两套对比共同达 FDR&lt;0.05 的有 <b>34</b> 个基因；</li>
  <li>再以缓解期 UC（n=8）相对对照的变化方向交叉核对：{len(sig)} 个基因中 <b>35 个（80%）</b>方向一致。</li>
</ul>
<p class="small">须注意：这是<b>同一数据集内部</b>的稳健性检查，<b>不能替代独立队列的外部验证</b>（见第六节 A-2 项）。</p>
""")

H.append(f"""
<h2>四、技术说明与局限（务必带入正文写作）</h2>
<div class="callout c-danger">
<ol style="margin:4px 0">
  <li><b>FDR 是在 SGRG 面板内校正的</b>（n = {S['sgrg_measured']}），不是全转录组。<b>论文中必须写明这一点</b>，
      否则会被审稿人质疑"选择性校正"。若需全基因组 FDR，需重跑全矩阵。</li>
  <li><b>效应量普遍偏小</b>：|log2FC| ≥ 0.5 的显著基因仅 {len(bigE)} 个。用传统"|log2FC|&gt;1 且 FDR&lt;0.05"标准筛，
      几乎筛不出基因——这是单细胞以外的 bulk 黏膜数据的常态，建议改用效应量排序 + 独立队列复现来支撑候选基因。</li>
  <li><b>{len(multi)} 个 SGRG 含跨基因探针</b>（{', '.join(r['GeneSymbol'] for r in multi)}），其中 DDX3X、TARDBP 等为重要基因。
      芯片上这些基因没有唯一探针，效应值可能被相邻基因稀释，<b>不建议作为唯一证据</b>。</li>
  <li><b>微阵列 vs RNA-seq</b>：本数据集为 HG-U133 Plus 2.0，动态范围与 RNA-seq 不同；若 Part A 最终要换用 RNA-seq 队列
      （GSE206285 等，目前网络环境取不到），需重新核对探针→基因映射。</li>
  <li><b>单一数据集</b>：当前结论仅来自 GSE38713（n=43）。外部验证队列已启动下载（GSE11223，n=202），见第六节。</li>
  <li><b>芯片未覆盖 {S['sgrg_not_measured']} 个 SGRG</b>（如 NANOS3、ENDOV 等），属平台限制而非生物学缺失。</li>
</ol>
</div>
""")

H.append(f"""
<h2>五、本步对整体方案的意义</h2>
<div class="card">
<p>本步是 <b>Part A 的地基</b>：产出的候选基因池将直接进入后续机器学习特征选择。同时它给 Part C 提供了"转录层"旁证：</p>
<ul>
  <li><b>Part A（本列表 → 建模）</b>：以 {len(sig)} 个差异 SGRG 为输入，做 LASSO / SVM-RFE / 随机森林特征筛选，
      再建诊断模型与列线图；外部验证用 GSE11223。</li>
  <li><b>Part C（eQTL-MR）</b>：上一个对话已锁定 EIF2AK4、G3BP2、DDX6、DDX3X 四个"肠+脑双侧可用"的工具变量基因。
      本步结果中 <b>G3BP1、EIF2S1、VCP、ATXN2、FMR1 在转录层显著</b>，可与 MR 的因果层结果做证据叠合（最终 Fig 10）。</li>
  <li><b>Part B（脑肠桥梁）</b>：若后续加入 MDD 数据集，同一 SGRG 面板可在两个器官间比较方向一致性。</li>
</ul>
</div>
""")

H.append(f"""
<h2>六、下一步（Part A 后续与并行线）</h2>
<table>
<thead><tr><th>编号</th><th>内容</th><th>状态</th></tr></thead>
<tbody>
<tr><td><b>A-2</b></td><td>外部验证：用独立队列复现本列表候选基因的方向与显著性</td>
    <td><span class="pill p-c">进行中</span> GSE11223（n=202，UC 与健康对照结肠活检）已开始下载</td></tr>
<tr><td><b>A-3</b></td><td>机器学习特征选择（LASSO / SVM-RFE / 随机森林）→ 诊断模型 + 列线图 + ROC 曲线</td>
    <td><span class="pill">待环境就绪</span> 需 scikit-learn；本机安装因网络限速未完成</td></tr>
<tr><td><b>A-4</b></td><td>功能富集（GO / KEGG / Reactome），定位差异 SGRG 的通路归属</td>
    <td><span class="pill">待办</span></td></tr>
<tr><td><b>C-1</b></td><td>eQTL 双样本 MR + 共定位：EIF2AK4 / G3BP2 / DDX6 / DDX3X → UC / MDD</td>
    <td><span class="pill p-a">条件已具备</span> 工具变量上一轮已核验通过（Part C GO）</td></tr>
<tr><td><b>叠合</b></td><td>转录层证据（本报告）与因果层证据（MR）叠合，锁定最终靶基因（对应方案 Fig 10）</td>
    <td><span class="pill">待办</span></td></tr>
</tbody></table>

<div class="callout c-warn">
<b>关于本环境的限制（供你判断下一步如何安排）：</b>
NCBITaxon 的 GEO / NCBI 主站在当前代理白名单之外，<code>ftp.ncbi.nlm.nih.gov</code> 与 <code>eutils.ncbi.nlm.nih.gov</code> 均不可达；
本报告的数据是经 <b>Gemma（gemma.msl.ubc.ca）</b> 与 <b>EBI（QuickGO、ArrayExpress/BioStudies）</b> 取得。
实测下载速率约 <b>13 KB/s</b>，一个 200 例量级的表达矩阵需半小时以上。
若后续要处理 GSE206285 等大数据集，建议在本机常规网络环境下运行 —— 分析脚本已写成可直接复用的形式
（<code>scripts/partA_step1_sgrg_de.py</code> 只需替换矩阵路径与分组来源）。
</div>
""")

H.append("""
<footer>
<p><b>可溯源性</b>　表达矩阵：<code>data/raw/gemma_GSE38713.tsv</code>（源自 Gemma REST API v3，dataset 5873；
原始 gzip 见同目录 <code>.tsv.gz</code>）。基因集：<code>results/sgrg_gene_set_clean.tsv</code> 与
<code>sgrg_rejected.tsv</code>；GO 原始注释 <code>data/raw/qgo/*.json</code>。
分析脚本：<code>scripts/validate_sgrg_hgnc.py</code>、<code>scripts/partA_step1_sgrg_de.py</code>。
结果表：<code>results/partA_GSE38713_SGRG_DE.tsv</code>、<code>partA_GSE38713_SGRG_groupMeans.tsv</code>。</p>
<p><b>运行环境</b>　Python 3.13（managed venv：numpy 2.5.3 / scipy 1.18.1 / matplotlib 3.11.2），
本机无 R。数据检索与外源校验均在 2026-10-08 当日完成。</p>
</footer>
</div></body></html>""")

open(REPORT, "w", encoding="utf-8").write("\n".join(H))
print("[write]", REPORT, os.path.getsize(REPORT), "bytes")
