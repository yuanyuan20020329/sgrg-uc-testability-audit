# -*- coding: utf-8 -*-
"""Part A · Step 2 报告生成器（综合 2/2b/2c/2d 四项分析）"""
import os, json, base64, datetime

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT = os.path.join(BASE, "results")
log = print

def J(name, default=None):
    p = os.path.join(OUT, name)
    if not os.path.exists(p): return default
    return json.load(open(p, encoding="utf-8"))

def b64(p):
    if not os.path.exists(p): return None
    return "data:image/png;base64," + base64.b64encode(open(p, "rb").read()).decode()

def fig(p, cap, maxw="88%"):
    d = b64(p)
    if not d: return f"<p class='miss'>[图缺失: {os.path.basename(p)}]</p>"
    return f"<figure><img src='{d}' style='width:{maxw}'/><figcaption>{cap}</figcaption></figure>"

def n3(x): return "—" if x is None else f"{x:.3f}"
def pc(x): return "—" if x is None else f"{x*100:.1f}%"

M = J("partA_step2_model_metrics.json", {})
S = J("partA_step2b_panel_specificity.json", {})
G = J("partA_step2c_global_structure.json", {})
B = J("partA_step2d_batch_corrected.json", {})

g13 = G.get("GSE38713", {})
s13 = S.get("GSE38713", {})
m13 = M.get("GSE38713", {})

html = []
html.append("""<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">
<title>Part A · Step 2 — UC 应激颗粒相关基因诊断模型（含批次混杂自检）</title>
<style>
:root{--ink:#1a1d21;--sub:#5b6470;--line:#e3e7ec;--accent:#0b5cab;--warn:#b42318;--ok:#1a7f4b;--bg:#fbfcfd;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font-family:-apple-system,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;font-size:15px;line-height:1.75}
.wrap{max-width:960px;margin:0 auto;padding:44px 30px 80px}
h1{font-size:25px;margin:0 0 6px;letter-spacing:-.3px}
h2{font-size:19px;margin:40px 0 12px;padding-bottom:7px;border-bottom:2px solid var(--line)}
h3{font-size:15.5px;margin:24px 0 8px;color:var(--accent)}
p{margin:9px 0} small,.sub{color:var(--sub)}
.lede{color:var(--sub);font-size:14px;margin-top:2px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(158px,1fr));gap:12px;margin:20px 0}
.kpi{background:#fff;border:1px solid var(--line);border-radius:10px;padding:13px 15px}
.kpi .v{font-size:23px;font-weight:700;letter-spacing:-.5px}
.kpi .k{font-size:12px;color:var(--sub);margin-top:2px}
.kpi.bad .v{color:var(--warn)} .kpi.good .v{color:var(--ok)}
table{border-collapse:collapse;width:100%;margin:14px 0;font-size:13.5px;background:#fff}
th,td{border:1px solid var(--line);padding:8px 10px;text-align:left;vertical-align:top}
th{background:#f2f5f8;font-weight:600;font-size:12.5px}
td.n{text-align:right;font-variant-numeric:tabular-nums}
figure{margin:18px 0;background:#fff;border:1px solid var(--line);border-radius:10px;padding:12px;text-align:center}
figcaption{font-size:12.5px;color:var(--sub);margin-top:8px;text-align:left}
.box{background:#fff;border:1px solid var(--line);border-left:4px solid var(--accent);border-radius:8px;padding:13px 16px;margin:16px 0}
.box.warn{border-left-color:var(--warn);background:#fff8f7}
.box.ok{border-left-color:var(--ok);background:#f6fbf8}
.box h4{margin:0 0 6px;font-size:14px}
ul{margin:8px 0 8px 20px;padding:0} li{margin:5px 0}
code{background:#eef2f6;padding:1px 5px;border-radius:4px;font-size:12.5px}
footer{margin-top:44px;padding-top:14px;border-top:1px solid var(--line);color:var(--sub);font-size:12.5px}
.miss{color:var(--warn);font-size:12.5px}
</style></head><body><div class="wrap">""")

html.append("<h1>Part A · Step 2 — 溃疡性结肠炎（UC）应激颗粒相关基因诊断模型</h1>")
html.append("<p class='lede'>SGRG 面板判别力 · 面板特异性自检 · 批次混杂诊断 · 去批次稳健性</p>")
html.append(f"<p><small>生成时间 {datetime.datetime.now():%Y-%m-%d %H:%M} ｜ 队列 GSE38713（人结肠黏膜活检，n=43）"
            f"｜ 环境 Python 3.13 · scikit-learn 1.9</small></p>")

# ---------- 结论先行 ----------
html.append(f"""
<div class='box warn'><h4>结论先行（请先读这一段）</h4>
<p>模型本身"看起来"很好：SGRG 面板在 GSE38713 上的留一交叉验证 AUC = <b>{n3(m13.get('cv_auc'))}</b>，
置换检验 p = {m13.get('permutation_p', float('nan')):.4f}。但三项自检一致表明，
<b>这个判别力基本不能归因于应激颗粒基因本身</b>：</p>
<ol>
<li>同等规模的<b>随机基因面板</b>就能达到中位 AUC {n3(s13.get('rand_median'))}（SGRG 只排在第 {s13.get('sgrg_percentile','—')} 百分位）；</li>
<li>该数据集的<b>批次与疾病状态高度混杂</b>（Cramér's V = {n3(g13.get('batch_cramers_v'))}，p = {g13.get('batch_chi2_p', float('nan')):.4g}），
仅用批次标签预测疾病的 AUC 就达 {n3(g13.get('batch_only_auc'))}；</li>
<li><b>去掉批次效应后</b>，SGRG 面板 AUC 塌到 <b>{n3((B.get('batch_centered') or {}).get('auc_sgrg'))}</b>，与随机面板无实质差别。</li>
</ol>
<p>因此本步的正确结论是：<b>GSE38713 单队列无法支撑"SGRG 具有 UC 诊断价值"的论断</b>；
该数据集的分组差异主要由处理批次驱动。需要独立队列（且批次平衡）才能做可靠判断。</p></div>
""")

html.append("<div class='kpis'>")
html.append(f"<div class='kpi good'><div class='v'>{n3(m13.get('cv_auc'))}</div><div class='k'>SGRG 面板 LOO-CV AUC<br>(原始，未去批次)</div></div>")
html.append(f"<div class='kpi'><div class='v'>{s13.get('sgrg_percentile','—')}<span style='font-size:14px'>%</span></div><div class='k'>在随机 97 基因面板中的百分位<br>(即不特异)</div></div>")
html.append(f"<div class='kpi bad'><div class='v'>{n3(g13.get('batch_cramers_v'))}</div><div class='k'>疾病×批次 Cramér's V<br>(1=完全混杂)</div></div>")
html.append(f"<div class='kpi bad'><div class='v'>{n3((B.get('batch_centered') or {}).get('auc_sgrg'))}</div><div class='k'>去批次后 AUC<br>(≈随机水平)</div></div>")
html.append("</div>")

# ---------- 一、设计与数据 ----------
html.append("<h2>一、数据与设计</h2>")
html.append("""
<table>
<tr><th>项</th><th>内容</th></tr>
<tr><td>队列</td><td>GSE38713（Gemma dataset 5873）— 人结肠黏膜活检，43 例</td></tr>
<tr><td>分组（取自 Gemma 实验设计因子 DiseaseState，非按标题猜测）</td>
<td>非炎症对照 13 ｜ 缓解期 UC 8 ｜ 活动期 UC 非炎症区 7 ｜ 活动期 UC 炎症区 15</td></tr>
<tr><td>批次因子</td><td>7 个处理批次（Batch_01 7/24/07 … Batch_07 5/13/09），由 Gemma 从原始 GEO 的 "batch" 特征导入</td></tr>
<tr><td>特征</td><td>SGRG 面板（GO:0010494 ∪ GO:0035061，经 HGNC 校验的 104 个官方符号）中本芯片可定量的 97 个</td></tr>
<tr><td>平台</td><td>Affymetrix HG-U133 Plus 2.0；Gemma 已完成 RMA + 分位数标准化</td></tr>
</table>
<p><small>原计划用 GSE11223（n=202）作为独立验证队列，但本机到 Gemma 的连接在传输约 5.5 MB 后即被中断
（完整文件 20.8 MB，实测约 10 KB/s），多次重试未取得完整文件，故本轮未能纳入外部验证。</small></p>
""")

# ---------- 二、方法 ----------
html.append("<h2>二、方法</h2>")
html.append("""
<ul>
<li><b>模型</b>：L1 / L2 惩罚逻辑回归（liblinear）。正则强度 C 在训练集内部网格搜索、以交叉验证 AUC 择优。</li>
<li><b>交叉验证</b>：n ≤ 60 用留一法（LOO）；准确率与敏感度/特异度取 Youden 最优阈值。</li>
<li><b>置换检验</b>：随机打乱标签 200 次，重算交叉验证 AUC，给出经验 p 值。</li>
<li><b>面板特异性</b>：从同一芯片、同一预处理后的基因池中随机抽 200 个与 SGRG 等规模的基因面板，用完全相同的模型配置计算 AUC 分布。</li>
<li><b>批次混杂</b>：疾病 × 批次的列联表 Cramér's V；以及仅用批次 one-hot 作为特征预测疾病的交叉验证 AUC。</li>
<li><b>去批次</b>：逐基因减去其所在批次的均值（最简单的批次中心化），再重算面板 AUC 与随机面板分布。</li>
</ul>
""")

# ---------- 三、模型性能 ----------
html.append("<h2>三、模型的表观性能</h2>")
if m13:
    html.append("<table><tr><th>指标</th><th class='n'>值</th><th>说明</th></tr>")
    rows = [("交叉验证 AUC", n3(m13.get("cv_auc")), f"最优配置 {m13.get('penalty','').upper()}，C={m13.get('best_C')}"),
            ("L1 / L2 各自最优 AUC", f"{n3(m13.get('cv_auc_l1'))} / {n3(m13.get('cv_auc_l2'))}", "两种惩罚下的最优"),
            ("准确率", pc(m13.get("cv_acc")), "Youden 阈值"),
            ("敏感度", pc(m13.get("cv_sens")), "UC 召回率"),
            ("特异度", pc(m13.get("cv_spec")), "对照正确排除率"),
            ("置换检验 p", f"{m13.get('permutation_p'):.4f}", "200 次标签置换"),
            ("样本", f"{m13.get('n_train')}（UC {m13.get('n_pos')} / 对照 {m13.get('n_train',0)-m13.get('n_pos',0)}）", "")]
    for a, b, c in rows:
        html.append(f"<tr><td>{a}</td><td class='n'>{b}</td><td>{c}</td></tr>")
    html.append("</table>")
    html.append(fig(os.path.join(OUT, "fig_ROC_SGRG_diagnostic_model.png"),
                    "ROC：仅训练队列交叉验证曲线（本轮无外部队列）。AUC 看似很高，但需结合第四节的特异性检验解读。", "60%"))

# ---------- 四、签名 ----------
html.append("<h2>四、L1 稀疏签名</h2>")
if m13.get("signature"):
    html.append(f"<p>用 L1（LASSO）在其自身最优 C={m13.get('signature_C')} 上全量拟合，"
                f"97 个特征中保留 {m13.get('signature_n')} 个。系数为特征标准化后每 1 个标准差的 log-odds，"
                f"正号偏向 UC。</p>")
    html.append(fig(os.path.join(OUT, "fig_signature_GSE38713.png"), "L1 签名基因与系数。", "62%"))
    html.append("<table><tr><th>#</th><th>基因</th><th class='n'>系数</th><th>方向</th></tr>")
    for i, s in enumerate(m13["signature"], 1):
        html.append(f"<tr><td>{i}</td><td>{s['gene']}</td><td class='n'>{s['coef']:+.3f}</td>"
                    f"<td>{'↑ UC' if s['coef']>0 else '↓ UC'}</td></tr>")
    html.append("</table>")
    html.append("""<div class='box warn'><h4>不要把这个签名当作生物标志物</h4>
<p>n=43、p=97 时 L1 的选基因极不稳定，换一次重抽样就可能换一批基因；而且下一节的检验说明这些基因的
判别力并不来自应激颗粒生物学本身。签名只能作为"候选"，不能作为结论。</p></div>""")

# ---------- 五、面板特异性 ----------
html.append("<h2>五、面板特异性检验：这个 AUC 是 SGRG 特有的吗？</h2>")
if s13:
    html.append(f"""
<p>在<b>完全相同的基因池、相同的预处理、相同的模型配置</b>下，随机抽取 200 个与 SGRG 等规模（{s13.get('n_features')} 个）的基因面板：</p>
<table><tr><th>面板</th><th class='n'>交叉验证 AUC</th></tr>
<tr><td>SGRG 面板（97 基因）</td><td class='n'><b>{n3(s13.get('auc_sgrg'))}</b></td></tr>
<tr><td>随机面板 — 中位数</td><td class='n'>{n3(s13.get('rand_median'))}</td></tr>
<tr><td>随机面板 — 95 百分位</td><td class='n'>{n3(s13.get('rand_p95'))}</td></tr>
<tr><td>随机面板 — 最大值</td><td class='n'>{n3(s13.get('rand_max'))}</td></tr>
<tr><td>SGRG 所处百分位</td><td class='n'>{s13.get('sgrg_percentile')}%</td></tr>
</table>
<p>SGRG 面板确实优于随机（位于第 {s13.get('sgrg_percentile')} 百分位），但<b>优势很小</b>——
随机面板的中位数 {n3(s13.get('rand_median'))} 与它的 {n3(s13.get('auc_sgrg'))} 只差 {abs(s13.get('auc_sgrg',0)-s13.get('rand_median',0)):.3f}，
而随机面板的上分位甚至更高。这说明该数据集中"任意 97 个基因"都能把 UC 与对照分开，
提示存在<b>非特异性的全局结构</b>（见下一节）。</p>
""")
    html.append(fig(os.path.join(OUT, "fig_panel_specificity_GSE38713.png"),
                    "SGRG 面板（红线）与 200 个随机等规模面板的 AUC 分布。红线并不明显偏右。", "70%"))

# ---------- 六、批次混杂 ----------
html.append("<h2>六、全局结构与批次混杂诊断</h2>")
if g13:
    bc = g13.get("batch_cats", {})
    html.append(f"""
<h3>6.1 存在一条与疾病对齐的全局主轴</h3>
<p>对 43 个样本的全转录组做逐基因标准化后取主成分：<b>PC1 解释 {g13.get('pc1_var',0)*100:.1f}% 的方差，
单独一个 PC1 对"UC vs 对照"的 AUC 就达到 {n3((g13.get('pc_auc_vs_disease') or [None])[0])}</b>，
而 PC2–PC5 分别只有 {g13.get('pc_auc_vs_disease',[0,0,0,0,0])[1:]} 。也就是说，
数据里第一主成分几乎就是"疾病标签"。在真实生物学信号通常分散的情况下，这种"单一主轴=疾病"的形态，
更常见的解释是<b>技术因素</b>。</p>
<p>与之呼应：PC1 与样本在文件中的排列顺序（处理顺序的粗糙代理）存在中等相关（Spearman ρ = {n3(g13.get('pc1_order_spearman'))}）。</p>
""")
    html.append(fig(os.path.join(OUT, "fig_global_structure_GSE38713.png"),
                    "左：主成分方差贡献（PC1 独占鳌头）。右：PC1 取值按疾病状态分组。", "86%"))

    html.append(f"""
<h3>6.2 批次与疾病状态高度混杂</h3>
<p>GSE38713 的 7 个处理批次与疾病状态显著关联：<b>Cramér's V = {n3(g13.get('batch_cramers_v'))}
（χ² = {g13.get('batch_chi2',0):.1f}，df = {g13.get('batch_dof',0)}，p = {g13.get('batch_chi2_p',0):.4g}）</b>。</p>
<table><tr><th>批次</th><th class='n'>UC</th><th class='n'>对照</th><th>读法</th></tr>
""")
    for b in sorted(bc):
        r = bc[b]
        note = "纯 UC，无对照" if r["Control"] == 0 else ("几乎全是对照" if r["UC"] <= 1 else "混合")
        cls = " style='color:#b42318;font-weight:600'" if r["Control"] == 0 else ""
        html.append(f"<tr><td>{b}</td><td class='n'>{r['UC']}</td><td class='n'>{r['Control']}</td><td{cls}>{note}</td></tr>")
    html.append(f"""</table>
<p>7 个批次中有 <b>4 个批次完全没有对照样本</b>。在这种结构下，任何能捕捉批次差异的特征
都能顺带"预测"疾病——<b>仅用批次标签本身</b>（one-hot，不含任何基因表达）预测疾病的交叉验证
AUC 就达到 <b>{n3(g13.get('batch_only_auc'))}</b>。</p>
""")
    html.append(fig(os.path.join(OUT, "fig_batch_confound_GSE38713.png"),
                    "左：每个批次的 UC/对照构成（多个批次为纯 UC）。右：PC1 取值随批次分层。", "92%"))

# ---------- 七、去批次 ----------
html.append("<h2>七、去掉批次之后还剩什么</h2>")
if B:
    raw = B.get("raw", {}); bc2 = B.get("batch_centered", {})
    html.append(f"""
<p>逐基因减去其所在批次的均值后，用<b>完全相同的流程</b>重算：</p>
<table><tr><th>版本</th><th class='n'>SGRG 面板 AUC</th><th class='n'>随机面板中位数</th><th class='n'>随机面板 95 分位</th><th class='n'>SGRG 百分位</th></tr>
<tr><td>原始（未去批次）</td><td class='n'>{n3(raw.get('auc_sgrg'))}</td><td class='n'>{n3(raw.get('rand_median'))}</td><td class='n'>{n3(raw.get('rand_p95'))}</td><td class='n'>{raw.get('percentile')}%</td></tr>
<tr><td><b>去批次后</b></td><td class='n'><b style='color:#b42318'>{n3(bc2.get('auc_sgrg'))}</b></td><td class='n'>{n3(bc2.get('rand_median'))}</td><td class='n'>{n3(bc2.get('rand_p95'))}</td><td class='n'>{bc2.get('percentile')}%</td></tr>
</table>
<p>去批次后 SGRG 面板的 AUC 从 {n3(raw.get('auc_sgrg'))} 塌到 {n3(bc2.get('auc_sgrg'))}，几乎回到随机水平。
它仍略高于随机中位数（位于第 {bc2.get('percentile')} 百分位，说明残留了一点非零信号），
但 <b>AUC ≈ 0.58 在临床上没有任何判别价值</b>。</p>
<div class='box warn'><h4>这条检验的边界条件</h4>
<p>当疾病与批次混杂时，去批次会连同"被混杂的那部分疾病信号"一起移除，所以
{n3(bc2.get('auc_sgrg'))} 是判别力的<b>下界</b>而非真值。它不能证明"SGRG 与 UC 无关"，
只能证明<b>用这份数据无法把两者分开</b>。</p></div>
""")
    html.append(fig(os.path.join(OUT, "fig_batch_corrected_specificity_GSE38713.png"),
                    "左：原始数据——SGRG 红线与随机分布几乎重合。右：去批次后——SGRG 进一步贴近随机中位数。", "92%"))

# ---------- 八、结论与写作建议 ----------
html.append("<h2>八、结论与正文写作建议</h2>")
html.append("""
<div class='box'><h4>可以直接写的</h4>
<ul>
<li>在 GSE38713 中，活动期 UC 与对照之间存在明显的全转录组结构差异，SGRG 面板可达到较高交叉验证 AUC。</li>
<li>但该差异的相当一部分可归因于处理批次：批次与疾病状态显著混杂，且去批次后 SGRG 判别力大幅衰减。</li>
</ul></div>
<div class='box warn'><h4>不能写的</h4>
<ul>
<li><b>不要说</b>"SGRG 可作为 UC 诊断标志物"。当前证据不支持——同等规模的随机基因集表现相当。</li>
<li><b>不要说</b>签名里的 HECA / DAZAP2 / EIF4A1 等是"关键标志基因"。这些是 n=43 下 L1 的不稳定选择。</li>
<li><b>不要</b>把交叉验证 AUC 当作"模型性能"上报而不提批次问题，这在评审中很容易被质疑。</li>
</ul></div>
<div class='box ok'><h4>建议的下一步（按性价比排序）</h4>
<ol>
<li><b>换一个批次平衡的独立队列做验证</b>（如 GSE11223，n=202），在同预处理的 Gemma 矩阵上重跑本节全部脚本；
只有在"批次不混杂 + 外部队列仍成立"时，SGRG 的诊断价值才站得住。</li>
<li>若坚持用 GSE38713，应在正文中<b>显式报告批次因子</b>，并给出<b>去批次后的结果</b>作为敏感性分析，而不是只报原始 AUC。</li>
<li>把重心从"诊断模型"转到<b>机制/因果</b>（Part C 的孟德尔随机化），转录组部分改为"表达层面的支持性证据"，
并明确其混杂局限。</li>
</ol></div>
""")

html.append("<footer>Part A · Step 2 ｜ 脚本：<code>partA_step2_sgrg_model.py</code>、"
            "<code>partA_step2b_panel_specificity.py</code>、<code>partA_step2c_global_structure.py</code>、"
            "<code>partA_step2d_batch_corrected.py</code>、共享模块 <code>sgrg_common.py</code> ｜ "
            "结果 JSON：<code>partA_step2*_*.json</code></footer>")
html.append("</div></body></html>")

p = os.path.join(OUT, "PartA_Step2_SGRG诊断模型报告.html")
open(p, "w", encoding="utf-8").write("\n".join(html))
print("written:", p, os.path.getsize(p), "bytes")
