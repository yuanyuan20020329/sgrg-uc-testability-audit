# -*- coding: utf-8 -*-
"""生成 Part A2 报告：公共 UC 队列池的批次质控标注 × 下游数据侧后果（自包含 HTML）。"""
import base64
import csv
import json
import os

import pandas as pd

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT = os.path.join(BASE, "results")
REPORT = os.path.join(OUT, "PartA2_批次质控标注与下游后果报告.html")


def b64(p):
    return "data:image/png;base64," + base64.b64encode(open(p, "rb").read()).decode()


A = json.load(open(os.path.join(OUT, "partA2_analysis.json"), encoding="utf-8"))
rows = list(csv.DictReader(open(os.path.join(OUT, "partA2_cohorts.tsv"), encoding="utf-8"), delimiter="\t"))
ts = list(csv.DictReader(open(os.path.join(OUT, "partA2_threshold_sensitivity.tsv"), encoding="utf-8"), delimiter="\t"))
uc = list(csv.DictReader(open(os.path.join(OUT, "partA2_uc_subset.tsv"), encoding="utf-8"), delimiter="\t"))

land = {r["annot_status"]: r for r in csv.DictReader(open(os.path.join(OUT, "partA2_landscape.tsv"), encoding="utf-8"), delimiter="\t")}
ORDER = ["批次伪影", "信息不足", "无问题", "已校正"]
EN = {"批次伪影": "Batch artifact", "信息不足": "Insufficient batch info",
      "无问题": "No batch problem", "已校正": "Batch corrected"}

dd = A["discordance"]
mb = A["models_binary"]["auc_pc1_batch"]
mbio = A["models_binary"]["auc_pc1_primary"]
mnorm = A["models_binary"]["normalized_batch_auc"]
win = A["pc1_batch_wins"]
matched = A["matched"]
ucs = A["uc_subset"]
gee = {r["score"]: r for r in A["geeq_vs_svd"]}
sc = {r["col"] + ":" + r["label"]: r for r in A["status_contrasts"]}
mde = A.get("de_contrasts") or []
PS = A["paired_subset"]["medians"]
N_PAIRED = A["paired_subset"]["n"]
N_DE_MISSING = A["pool_summary"]["n"] - A["pool_summary"]["n_with_de_disease"]
_pst = pd.read_csv(os.path.join(OUT, "partA2_paired_subset.tsv"), sep="\t")
ps_stats = {bool(k): {"n": len(v), "bio": float(v["auc_pc1_primary"].median()),
                      "batch": float(v["auc_pc1_batch"].median())}
            for k, v in _pst.dropna(subset=["primary_is_disease"]).groupby("primary_is_disease")}
_pf = pd.read_csv(os.path.join(OUT, "partA2_cohorts.tsv"), sep="\t")
_pf = _pf.dropna(subset=["primary_is_disease"])
PCT_DISEASE_PF = 100.0 * _pf["primary_is_disease"].astype(bool).mean()
_art = _pf[_pf["annot_status"] == "批次伪影"]["primary_cat"].dropna()
PCT_ART_TREATMENT = 100.0 * (_art == "treatment").mean() if len(_art) else float("nan")
PCT_ART_DISEASE = 100.0 * (_art == "disease").mean() if len(_art) else float("nan")


def term(m, name):
    for t in (m or {}).get("terms", []):
        if t["term"] == name:
            return t
    return {}


def frow(*vals):
    return "<tr>" + "".join("<td>%s</td>" % v for v in vals) + "</tr>"


def num(x, fmt="%.3f"):
    try:
        v = float(x)
        return "n/a" if v != v else (fmt % v)
    except Exception:
        return str(x)


# ---- 表：各组 PC1 批次可分性
def grp_table(records, col_label):
    out = []
    for r in sorted(records, key=lambda x: -float(x["median"])):
        out.append(frow("<b>%s</b>" % EN.get(r["group"], r["group"]), r["n"],
                        num(r["median"]), "%s–%s" % (num(r["q1"]), num(r["q3"]))))
    return ("<table><thead><tr><th>标注类别</th><th class='num'>n</th><th class='num'>中位</th>"
            "<th class='num'>IQR</th></tr></thead><tbody>" + "".join(out) + "</tbody></table>")


F1 = b64(os.path.join(OUT, "fig_A2_1_pool_landscape.png"))
F2 = b64(os.path.join(OUT, "fig_A2_2_batch_on_pc1.png"))
F3 = b64(os.path.join(OUT, "fig_A2_3_biology_preserved.png"))
F4 = b64(os.path.join(OUT, "fig_A2_4_discordance.png"))
F5 = b64(os.path.join(OUT, "fig_A2_5_geeq_crosscheck.png"))
F6 = b64(os.path.join(OUT, "fig_A2_6_uc_subset.png"))

pool = A["pool_summary"]
n_art = int(land["批次伪影"]["n"]) if "批次伪影" in land else 0

CSS = open(os.path.join(OUT, "partA2_report.css"), encoding="utf-8").read() if os.path.exists(os.path.join(OUT, "partA2_report.css")) else """
:root{--bg:#f6f7f9;--card:#fff;--line:#e5e7eb;--ink:#1f2328;--sub:#57606a;
--accent:#7a3e9d;--accent-soft:#f3ecfa;--warn:#b4530a;--warn-soft:#fdf3e7;
--danger:#b42318;--danger-soft:#fdecea;--blue:#0b5cab;--blue-soft:#eaf2fb;
--up:#b42318;--dn:#0b5cab;--ok:#0b6b57;--ok-soft:#e7f4f0;--tag:#f0f2f5;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);line-height:1.75;font-size:15px;
font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif}
.wrap{max-width:1120px;margin:0 auto;padding:36px 28px 80px}
header.hero{background:linear-gradient(135deg,#4a2b73,#7a3e9d 55%,#9b6bc4);color:#fff;border-radius:14px;padding:30px 32px;margin-bottom:30px}
header.hero .kicker{font-size:12.5px;letter-spacing:.16em;text-transform:uppercase;opacity:.85;margin-bottom:10px}
header.hero h1{margin:0 0 10px;font-size:24px;line-height:1.4;font-weight:700}
header.hero p{margin:0;opacity:.94;font-size:14.5px;max-width:940px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:22px 0 6px}
.kpi{background:rgba(255,255,255,.14);border:1px solid rgba(255,255,255,.3);border-radius:10px;padding:11px 14px}
.kpi b{display:block;font-size:21px;line-height:1.3}
.kpi span{opacity:.88;font-size:12.3px}
h2{font-size:20px;margin:44px 0 14px;padding-left:14px;border-left:4px solid var(--accent)}
h3{font-size:16.5px;margin:24px 0 9px;color:var(--accent)}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:20px 22px;margin:14px 0}
figure{margin:16px 0;text-align:center}
figure img{max-width:100%;border:1px solid var(--line);border-radius:10px;background:#fff}
figcaption{font-size:12.8px;color:var(--sub);margin-top:8px;text-align:left}
table{width:100%;border-collapse:collapse;margin:14px 0;font-size:13.2px;background:#fff;border-radius:10px;overflow:hidden}
th,td{padding:7px 10px;border-bottom:1px solid var(--line);text-align:left}
th{background:#f0f3f5;font-weight:600}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
.callout{border-radius:10px;padding:14px 18px;margin:16px 0;border:1px solid;font-size:14.2px}
.callout b{display:block;margin-bottom:4px}
.ok{background:var(--ok-soft);border-color:#bfe3d5;color:#0b4d3f}
.warn{background:var(--warn-soft);border-color:#f0d4b5;color:#7a3e0a}
.danger{background:var(--danger-soft);border-color:#f2c4c0;color:#8a1d13}
.blue{background:var(--blue-soft);border-color:#c4dcf5;color:#0a4a8a}
.small{font-size:12.8px;color:var(--sub)}
code{background:#eef1f4;padding:1px 6px;border-radius:5px;font-size:12.6px;font-family:ui-monospace,Consolas,monospace}
footer{margin-top:50px;padding-top:18px;border-top:1px solid var(--line);font-size:12.6px;color:var(--sub)}
.pill{display:inline-block;padding:1px 9px;border-radius:999px;font-size:11.6px;font-weight:600;background:var(--tag)}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:16px}
@media(max-width:900px){.grid2{grid-template-columns:1fr}}
"""

HTML = f"""<!DOCTYPE html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Part A2 · 批次质控标注与下游后果</title><style>{CSS}</style></head><body><div class="wrap">

<header class="hero">
  <div class="kicker">Part A2 · Gemma 公共 UC 队列池系统性质控调查</div>
  <h1>批次质控标注能多大程度预测下游数据侧后果？<br>——160 个公共 UC 队列的全池审计</h1>
  <p>对 Gemma 中与溃疡性结肠炎相关的全部公共队列（n=160），在不下载任何表达矩阵的前提下，
  把「Gemma 官方批次质控标注」与「数据侧可计算的主轴后果」（SVD 载荷上的批次/生物学可分性竞争、存算 DE 信号规模）逐队列配对，
  量化标注的灵敏度/特异度，并检验标注问题的队列是否真的在主导轴上被批次污染。</p>
  <div class="kpis">
    <div class="kpi"><b>{pool['n']}</b><span>公共 UC 相关队列（全池）</span></div>
    <div class="kpi"><b>{n_art}（{100*n_art/pool['n']:.0f}%）</b><span>被标注存在批次伪影</span></div>
    <div class="kpi"><b>{pool['n_with_pc1_batch']}</b><span>可计算批次/生物学主轴竞争</span></div>
    <div class="kpi"><b>{num(mb['F']['p_perm'],'%.4f')}</b><span>调整后 批次伪影→PC1 批次可分性 p（置换）</span></div>
    <div class="kpi"><b>{num(dd['ppv'],'%.2f')}</b><span>标注 PPV（对“批次主导主轴”）</span></div>
    <div class="kpi"><b>{num(dd['spec'],'%.2f')}</b><span>标注特异度（对“批次主导主轴”）</span></div>
  </div>
</header>

<div class="callout danger" style="margin-top:0"><b>本报告的可信度边界：请先读这一条</b>
本报告的核心检验是「Gemma 的批次标注 vs Gemma 自己 PCA 算出的批次可分性」——<b>两条通道同源</b>，
因此它<b>不是独立验证，本质上是自洽性检查</b>。Gemma 官方文档已明确写明其批次因子「almost always correlates with one of the first three
principal components」，所以「标注为伪影的队列，其主成分确实被批次占据」这一条在文献里<b>已被预先声明</b>，不构成新发现。
本报告真正站得住的原创贡献只有 §6：<b>UC 可分析证据池的实际规模远小于标称规模</b>。请按这个定位引用。</div>

<h2>方法（一页说清）</h2>
<div class="card">
<p><b>数据通道。</b>全池来自 Gemma REST <code>/datasets?query=ulcerative colitis</code>（并集 <code>colitis</code>），
逐队列抓取：队列元数据与官方 QC 字段（<code>batchEffect / batchConfound / batchEffectStatistics / GEEQ</code>）、
存算差异表达分析（<code>/analyses/differential</code>，不下载矩阵）、
SVD 载荷（<code>/svd</code>，159/160 可得）、样本级结构化因子值（<code>/samples</code> 的 <code>factorValues</code>）。</p>
<p><b>标注侧（自变量）。</b>按 Gemma 官方 <code>batchEffect</code> 枚举粗分为四类：
<b>批次伪影</b>（BATCH_EFFECT_FAILURE，算法在数据里查到批次伪影）、<b>信息不足</b>（无批次信息/表头不可用/单例批次等，n=56）、
<b>无问题</b>（无批次效应或单批次）、<b>已校正</b>（BATCH_CORRECTED_SUCCESS）。</p>
<p><b>后果侧（因变量，全部不下载表达矩阵）。</b>用 SVD 前几轴载荷 + 样本因子值计算：</p>
<ul>
<li><b>AUC<sub>batch</sub></b>：第 1 主成分对<b>批次因子</b>（block/scan date 等）的方向无关 AUC（秩和法）；</li>
<li><b>AUC<sub>biology</sub></b>：同一第 1 主成分对<b>队列自身生物学主因子</b>的 AUC。主因子优先取 Gemma 存算 DE
所指定的因子（外部选定、避免挑选），否则按类别优先级（disease &gt; treatment &gt; genotype &gt; …）；</li>
<li><b>批次主导</b>：AUC<sub>batch</sub> ≥ 0.75 且高于 AUC<sub>biology</sub> + 0.05（阈值敏感性见 §4）；</li>
<li><b>de_rate</b>：存算 DE 的差异表达探针占比（仅 {A['pool_summary']['n_with_de_disease']} 个队列有疾病对比，作次要结局）。</li>
</ul>
<p><b>统计。</b>组间中位差 + 95% bootstrap CI（4000 次）+ Mann-Whitney；调整模型为 OLS
（协变量：log2 样本量、RNA-seq、双通道芯片、肠道组织、PC1 方差占比、SVD 样本数），
对标注分类项做置换 F 检验（3000 次）；另做 log2n+PC1 方差最近邻配对的 Wilcoxon 检验。</p>
</div>

<h2>§1 队列池长什么样：先看地景</h2>
<div class="card">
<figure><img src="{F1}">
<figcaption><b>图 1.</b> a）160 个公共 UC 相关队列的官方批次标注构成：25% 被算法查出批次伪影，
35% 连可用的批次信息都没有，只有 29% 判定为无批次问题，11% 已做批次校正。
b）GEEQ 数据侧综合质量分按标注类别的分布（标注差的队列质量分明显更低）。
c）被官方标注「实验因子与批次混杂」的比例：批次伪影组 82%，无问题组 11%。</figcaption></figure>
<div class="callout warn"><b>直接读数</b>
一个“UC 公共队列池”里，<b>六成（{100*(n_art+int(land['信息不足']['n']))/pool['n']:.0f}%）的队列在批次维度上要么已查出伪影、要么连批次信息都不可用</b>。
任何把这个池子当整体做 meta 或跨队列聚合的下游分析，都建立在这 60% 的队列之上。</div>
</div>

<h2>§2 核心结果：标注「批次伪影」的队列，主导轴真的被批次占据</h2>
<div class="card">
<figure><img src="{F2}">
<figcaption><b>图 2.</b> 配对子集（两个主轴指标均可算，n={N_PAIRED}）内，第 1 主成分对批次因子的可分性（AUC，方向无关）：
批次伪影组中位 {num(PS['批次伪影']['auc_batch'])}，无问题组 {num(PS['无问题']['auc_batch'])}，
已校正组 {num(PS['已校正']['auc_batch'])}——已校正组显著低于伪影组（Mann-Whitney p = {num(sc['auc_pc1_batch:批次伪影 vs 已校正']['p_mwu'],'%.2e')}），
说明 Gemma 的批次校正在数据侧确实起效（阳性对照通过）。
b）同一主成分上「批次 − 生物学」的可分性差值：伪影组中位为正（{num(PS['批次伪影']['auc_diff'],'%+.3f')}），其余三组为负。</figcaption></figure>

<h3>四组主轴指标（配对子集，n={N_PAIRED}）</h3>
<table><thead><tr><th>标注类别</th><th class="num">n</th><th class="num">AUC<sub>batch</sub></th><th class="num">AUC<sub>biology</sub></th><th class="num">差值</th><th class="num">批次胜出</th></tr></thead><tbody>
{''.join(frow('<b>%s</b>' % EN[g], PS[g]['n'], num(PS[g]['auc_batch']), num(PS[g]['auc_bio']),
              '<b>%+.3f</b>' % PS[g]['auc_diff'] if PS[g]['auc_diff'] > 0 else '%+.3f' % PS[g]['auc_diff'],
              '%.0f%%' % (100*PS[g]['frac_batch_wins']))
         for g in ORDER if g in PS)}
</tbody></table>

<div class="grid2"><div>
<h3>组间对比（中位差，95% bootstrap CI）</h3>
<table><thead><tr><th>结局</th><th>对比</th><th class="num">中位差</th><th class="num">95% CI</th><th class="num">p</th></tr></thead><tbody>
{''.join(frow(r['col'], r['label'], '<b>%+.3f</b>' % r['diff'],
              '%+.3f ~ %+.3f' % (r['ci'][0], r['ci'][1]),
              '<b>%.2e</b>' % r['p_mwu'] if r['p_mwu'] < 0.001 else '%.3g' % r['p_mwu'])
         for r in A['status_contrasts'] if r['col'] in ('auc_pc1_batch','normalized_batch_auc'))}
</tbody></table>
</div><div>
<h3>调整模型（OLS + 置换 F）</h3>
<table><thead><tr><th>结局</th><th class="num">n</th><th class="num">R²</th><th class="num">批次伪影 β</th><th class="num">95% CI</th><th class="num">p（置换）</th></tr></thead><tbody>
{frow('AUC<sub>batch</sub> (PC1)', mb['n'], num(mb['R2']), '<b>%+.3f</b>' % term(mb,'annot_status[批次伪影]')['beta'],
      '%+.3f ~ %+.3f' % (term(mb,'annot_status[批次伪影]')['ci_lo'], term(mb,'annot_status[批次伪影]')['ci_hi']),
      '<b>%s</b>' % num(mb['F']['p_perm'],'%.4f'))}
{frow('AUC<sub>biology</sub> (PC1)', mbio['n'], num(mbio['R2']), '%+.3f' % term(mbio,'annot_status[批次伪影]')['beta'],
      '%+.3f ~ %+.3f' % (term(mbio,'annot_status[批次伪影]')['ci_lo'], term(mbio,'annot_status[批次伪影]')['ci_hi']),
      num(mbio['F']['p_perm'],'%.3f'))}
{frow('AUC<sub>batch</sub> − AUC<sub>bio</sub>', mnorm['n'], num(mnorm['R2']), '%+.3f' % term(mnorm,'annot_status[批次伪影]')['beta'],
      '%+.3f ~ %+.3f' % (term(mnorm,'annot_status[批次伪影]')['ci_lo'], term(mnorm,'annot_status[批次伪影]')['ci_hi']),
      num(mnorm['F']['p_perm'],'%.3f'))}
</tbody></table>
<p class="small">协变量：log2 样本量、RNA-seq、双通道芯片、肠道组织、PC1 方差占比、SVD 样本数。</p>
</div></div>

<div class="grid2"><div>
<h3>匹配对照（log2n + PC1 方差最近邻配对）</h3>
<table><thead><tr><th>对比</th><th class="num">对子</th><th class="num">中位差</th><th class="num">Wilcoxon p</th><th class="num">上升侧</th></tr></thead><tbody>
{frow('伪影 vs 无问题', matched['batch_auc_伪影vs无问题']['n_pairs'], '%+.3f' % matched['batch_auc_伪影vs无问题']['median_diff'],
      '%.4g' % matched['batch_auc_伪影vs无问题']['wilcoxon_p'],
      '%d/%d' % (matched['batch_auc_伪影vs无问题']['n_win'], matched['batch_auc_伪影vs无问题']['n_pairs']))}
{frow('<b>伪影 vs 已校正</b>', matched['batch_auc_伪影vs已校正']['n_pairs'], '<b>%+.3f</b>' % matched['batch_auc_伪影vs已校正']['median_diff'],
      '<b>%.2e</b>' % matched['batch_auc_伪影vs已校正']['wilcoxon_p'],
      '<b>%d/%d</b>' % (matched['batch_auc_伪影vs已校正']['n_win'], matched['batch_auc_伪影vs已校正']['n_pairs']))}
</tbody></table>
<p class="small">伪影 vs 已校正的 13 对里 <b>13/13 全部方向一致</b>——这不是混杂能解释的。</p>
</div><div>
<h3>PC1 上「批次胜过生物学」的比例</h3>
<table><thead><tr><th>组</th><th class="num">胜出/总数</th><th class="num">比例</th></tr></thead><tbody>
{frow('问题组（伪影+信息不足）', '%d / %d' % (win['problem'][0], sum(win['problem'])), '<b>%.1f%%</b>' % (100*win['frac_problem']))}
{frow('干净组（无问题+单批次+已校正）', '%d / %d' % (win['clean'][0], sum(win['clean'])), '%.1f%%' % (100*win['frac_clean']))}
</tbody></table>
<p class="small">Fisher 精确 p = {num(win['fisher']['p'],'%.4g')}。</p>
</div></div>
</div>

<h2>§3 意外但重要：批次并没有「压住」生物学信号</h2>
<div class="card">
<figure><img src="{F3}">
<figcaption><b>图 3.</b> 第 1 主成分对队列自身生物学主因子的可分性，四组之间没有差异
（调整模型 p = {num(mbio['F']['p_perm'],'%.2f')}；所有两两对比 p &gt; 0.9）。</figcaption></figure>
<div class="callout warn"><b>措辞修正：这里测的是「研究自己指定的主对比」，不是「生物学」</b>
全池里只有 <b>{num(PCT_DISEASE_PF,'%.0f')}%</b> 的队列其主因子是疾病因子；在批次伪影组，出现频率最高的主因子是
<b>treatment（{num(PCT_ART_TREATMENT,'%.0f')}%）</b>，其次才是 disease（{num(PCT_ART_DISEASE,'%.0f')}%）。
所以 AUC<sub>biology</sub> 严格说是「研究设计所声明的那个主对比能否在主成分上分开」，对 genotype / 细胞类型 / 处理这类强对比往往更容易分开
（非疾病主因子中位 AUC {num(ps_stats[False]['bio'])} vs 疾病主因子 {num(ps_stats[True]['bio'])}）。
把这一条写成「生物学信号未被压制」是过度解读，正确表述是「研究自身的核心对比仍可分」。</div>

<div class="callout blue"><b>这个阴性结果其实是信息量最大的部分</b>
直觉预期是「批次伪影会淹没生物学信号」。数据不支持：<b>批次伪影组队列的研究主对比在 PC1 上依然分得开
（配对子集内中位 AUC {num(PS['批次伪影']['auc_bio'])}，与已校正组 {num(PS['已校正']['auc_bio'])} / 无问题组 {num(PS['无问题']['auc_bio'])} 无差异）</b>。
真正的后果是结构性的——批次叠加在生物学之上、与它<b>竞争同一批主成分</b>（图 2b，伪影组「批次−生物学」差值转正）。
而几何上的竞争意味着：任何只看 PC1/PC2 的 QC、聚类、或者把主成分当作协变量的做法，
在伪影组里都会把技术结构误读成生物学结构。<b>受影响的不是「能不能测到信号」，而是「信号落在哪个轴上、和什么混在一起」。</b></div>
</div>

<h2>§4 标注到底可不可信：与数据侧后果的一致性</h2>
<div class="card">
<figure><img src="{F4}">
<figcaption><b>图 4.</b> a）把「批次主导主轴」的判定阈值从 0.65 扫到 0.90：标注的<b>特异度稳定在 0.90–0.95</b>，
灵敏度 0.27–0.49。b）阈值 0.75 下的 2×2 表：标注阳性而数据侧不成立的 {dd['table']['annot_only']} 例，
标注阴性而数据侧成立的仅 {dd['table']['svd_only']} 例。
判定口径为 <b>PC1–3 中的最佳批次可分性</b>（与 Gemma 官方判定标准一致；改用仅 PC1 的旧口径会人为把灵敏度压到 0.37）。</figcaption></figure>
<table><thead><tr><th class="num">阈值</th><th class="num">一致阳性</th><th class="num">仅标注阳性</th><th class="num">仅数据阳性（漏检）</th><th class="num">一致阴性</th><th class="num">灵敏度</th><th class="num">特异度</th><th class="num">PPV</th></tr></thead><tbody>
{''.join(frow(r['thr'], r['both_pos'], r['annot_only'], '<b>%s</b>' % r['svd_only'] if int(r['svd_only']) else '0',
              r['both_neg'], num(r['sens'],'%.2f'), num(r['spec'],'%.2f'), num(r['ppv'],'%.2f')) for r in ts)}
</tbody></table>
<div class="callout warn"><b>口径修正记录（2026-10-09）——这条必须自己先讲清楚</b>
本报告初版用「仅 PC1」判定批次主导，得出「特异 0.95 / 灵敏 0.37，标注严重过报」的结论。
但 Gemma 官方文档明确写道：批次因子「almost always correlates with one of the <b>first three</b> principal components」。
改用与其一致的 PC1–3 最佳口径后，<b>灵敏度升到 0.49、特异 0.90</b>，先前那句「大量被标注为伪影的队列主轴其实没被批次占据」<b>不成立</b>——
原结论是我自己口径选择造成的伪影，已撤回。</div>
<div class="callout ok"><b>修正后的结论：标注是个「够用但不惊人」的代理</b>
灵敏 ~0.49、特异 ~0.90、PPV ~0.91——既不是可靠的验收标准，也不是无用的噪声。按标注剔除会误伤约一半本可用的队列，
而放宽不用又会漏掉约 10% 的真问题。实用建议仍是：<b>把标注当触发数据侧复核的信号，而不是直接剔除的依据</b>，
但理由从「过报严重」改为「约半数被标记者并未被数据侧证实」。</div>
</div>

<h2>§5 独立交叉验证</h2>
<div class="card">
<figure><img src="{F5}">
<figcaption><b>图 5.</b> GEEQ 数据侧综合质量分 vs 我们计算的 PC1 批次可分性（Spearman ρ = {num(gee['geeq_public']['rho'])}，p = {num(gee['geeq_public']['p'],'%.4g')}，n = {gee['geeq_public']['n']}）。
质量分越低、批次可分性越高，方向符合预期；两条通道（GEEQ 的相关矩阵法 vs 我们的主轴秩和法）相互印证。</figcaption></figure>
</div>

<h2>§6 那 UC 本身呢：能真正回答 UC 问题的队列极少</h2>
<div class="card">
<figure><img src="{F6}">
<figcaption><b>图 6.</b> 只看样本级因子值里同时具备「UC ≥ 4 且对照 ≥ 4」的队列：<b>全池 160 个里只有 {ucs['n']} 个</b>。
分布：批次伪影 {ucs['by_status'].get('批次伪影',0)}、信息不足 {ucs['by_status'].get('信息不足',0)}、
无问题 {ucs['by_status'].get('无问题',0)}、已校正 {ucs['by_status'].get('已校正',0)}。
疾病主轴可分性（UC vs 对照）在已校正/无问题组最高（0.94 / 0.84），批次伪影组最低（0.77）。</figcaption></figure>
<div class="callout warn"><b>对整个项目最要紧的一条</b>
「153 个公共 UC 队列」这个数字在检索意义上成立，但在<b>可分析意义</b>上缩水得厉害：
能同时拿到 UC 与对照两臂（每臂 ≥4）的只有 {ucs['n']} 个，其余大多是无对照的 UC 内部对比、或被「ulcerative colitis」关键词带进来的
Crohn / IBS / 结直肠癌队列。<b>任何想把这个池子当独立验证集或聚合证据池的下游计划，都要先按这个口径重估。</b></div>
</div>

<h2>§7 局限与循环性声明（必须写进正文的部分）</h2>
<div class="card">
<ul>
<li><b>部分循环性（最要紧的一条）。</b>Gemma 的 <code>batchEffect</code> 标注（"This data set may have a batch artifact (PC k)"）由同一批数据的 PCA 算出，
且官方文档已声明批次因子与 PC1–3 之一相关；本报告的数据侧指标也从同一 PCA 的载荷算出。
因此 §2 与 §4 <b>不是独立验证，而是自洽性检查</b>——这一点无法通过调整口径消除，只能靠外部通道部分弥补：
已校正组的主轴回落（算法改了数据、后果随行为改变）与 GEEQ 交叉通道（§5，不同算法）。
若审稿人要求「独立验证」，唯一出路是<b>自己下载原始数据、用 Gemma 之外的批次检测方法重跑</b>，本报告未做。</li>
<li><b>关键词归因未经核验。</b>160 个队列来自检索词 <code>ulcerative colitis</code> 的命中，
其中混入 Crohn / IBS / 结直肠癌队列（§6 已观察到）。已发表的全库审计显示，关键词式疾病归因的失效率可达三成量级——
本报告<b>未做 MeSH / 发表文献级别的归因核验，也未做供体级去重</b>（同一批患者可能以多个 accessions 重复入库）。
因此 §6 的「15/160」应理解为「按最宽松口径的乐观上限」，真实可分析规模可能更小。</li>
<li><b>残留混杂。</b>问题组队列显著更大（中位 +11 样本，p = 8e-5）、更早（中位早 9 年，p = 0.002）——多批次结构本就需要更多样本。
所有主结论都在调整了 log2n 与 PC1 方差占比后成立，且匹配对照不依赖这些调整。</li>
<li><b>主因子的选择。</b>优先用 Gemma 存算 DE 指定的因子以避免挑选；当 DE 不可用时按类别优先级回退。
疾病对比缺失的队列（{N_DE_MISSING}/{pool['n']}）不进入 AUC<sub>biology</sub> 的疾病解释，只进入一般生物学因子的解释。</li>
<li><b>存算 DE 只覆盖 {A['pool_summary']['n_with_de_disease']} 个队列</b>，且 DE 率在问题/干净组间无差异（p = 0.53）——DE 信号规模不是批次伪影的可靠读出，
真正可靠的读出是主轴几何。</li>
<li><b>UC 子集 n = {ucs['n']}</b>，其分组比较只有描述意义，不做推断。</li>
</ul>
</div>

<h2>结论</h2>
<div class="card">
<ol>
<li><b>标注与数据侧后果方向一致且可量化</b>：批次伪影组 PC1 批次可分性中位 {num(PS['批次伪影']['auc_batch'])}，比已校正组高
{num(PS['批次伪影']['auc_batch'] - PS['已校正']['auc_batch'],'%.2f')}
（调整后 β = %+.2f，置换 p = {num(mb['F']['p_perm'],'%.4f')}；匹配对 13/13）。</li>
<li><b>伤害是几何性的，不是信号强度的</b>：生物学信号在 PC1 上仍然分得开，但与批次竞争同一主成分——
受影响的是「信号落在哪个轴、与什么混杂」，而非「测不测得到」。</li>
<li><b>标注是「够用但不惊人」的代理</b>：灵敏 {num(dd['sens'],'%.2f')}、特异 {num(dd['spec'],'%.2f')}、PPV {num(dd['ppv'],'%.2f')}。
按标注剔除会误伤约一半本可用的队列，放宽又会漏掉约 10% 的真问题。实用做法是把标注当<b>触发数据侧复核</b>的信号。</li>
<li><b>可分析口径下 UC 证据池大幅缩水</b>：160 → 仅 {ucs['n']} 个队列同时具备 UC 与对照两臂（每臂 ≥4）——这是本报告最站得住的原创贡献。</li>
<li><b>给 Part B / C 的直接含义</b>：任何「跨队列聚合 + 主成分/聚类 + 差异表达」的下游步骤，
在批次伪影队列里读出的“生物学结构”有实质概率是批次结构。此前 Part B 的面板级阴性结论，
在批次污染未校正的队列子集里并不能归因于生物学缺失——这是下一步（敏感性分析限定在已校正/无问题队列）的明确动机。</li>
</ol>
</div>

<footer>
数据通道：Gemma REST v2（元数据 / 存算 DE / SVD / 样本因子值），全程未下载表达矩阵。
脚本：<code>scripts/partA2_harvest_meta.py</code> → <code>partA2_build_table.py</code> →
<code>partA2_analyse.py</code> → <code>partA2_figures.py</code>；
宽表 <code>results/partA2_cohorts.tsv</code>（160 × 84），配对子集 <code>results/partA2_paired_subset.tsv</code>，主结果 <code>results/partA2_analysis.json</code>。
</footer>
</div></body></html>"""

open(REPORT, "w", encoding="utf-8").write(HTML)
print("写出", REPORT, "(%.1f KB)" % (os.path.getsize(REPORT) / 1024))
