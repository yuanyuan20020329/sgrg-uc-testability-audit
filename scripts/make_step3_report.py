# -*- coding: utf-8 -*-
"""Part A · Step 3 报告生成器（自包含 HTML，图片内嵌 base64）"""
import os, sys, json, csv, base64, html
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sgrg_common as C

OUT = C.OUT
S = json.load(open(os.path.join(OUT, "partA_step3_summary.json"), encoding="utf-8"))
PS = json.load(open(os.path.join(OUT, "partA_step3_panel_specificity.json"), encoding="utf-8"))
CC = json.load(open(os.path.join(OUT, "partA_step3_cross_cohort.json"), encoding="utf-8"))
DE = list(csv.DictReader(open(os.path.join(OUT, "partA_step3_paired_DE.tsv"), encoding="utf-8"),
                         delimiter="\t"))

P = S["paired_DE"]; R = S["replication"]; BD = S["batch_diag"]
c226, c387 = S["cohort_22619"], S["cohort_38713"]


def img(name, cap=""):
    p = os.path.join(OUT, name)
    if not os.path.exists(p):
        return ""
    b = base64.b64encode(open(p, "rb").read()).decode()
    c = f"<figcaption>{cap}</figcaption>" if cap else ""
    return f"<figure><img src='data:image/png;base64,{b}'>{c}</figure>"


def n(x, d=3):
    try:
        return f"{float(x):.{d}f}"
    except Exception:
        return "—"


def pct(x):
    try:
        return f"{float(x)*100:.1f}%"
    except Exception:
        return "—"


H = []
A = H.append
A("""<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">
<title>Part A · Step 3｜换批次验证（GSE22619 双生子配对队列）</title>
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

A(f"<h1>Part A · Step 3｜换批次验证：UC 不一致型双生子队列</h1>")
A(f"<div class='sub'>GSE22619（10 对同卵双生子 × 结肠黏膜）作为独立验证队列 · "
  f"生成于 2026-10-08 · 与 GSE38713 双向交叉队列验证</div>")

# ---------------- 结论先行
A("<h2>一、结论先行</h2>")
A(f"""<div class='warn'><h4>本次换批次把两件事分开来了：<span class="green">基因层面的位移方向真能复现</span>，
但<span class="red">「SGRG 面板可用于诊断」依然不成立</span>。</h4>
<p>换到批次与疾病完全正交的队列后：</p>
<ul>
<li><b class="green">方向复现成立</b>：GSE38713 中 FDR&lt;0.05 的 {R['n_GSE38713_sig']} 个基因，在配对设计中
{R['n_same_direction']} 个（{R['pct']}%）方向一致，二项检验 p={R['binom_p']:.2g}。全部 {c226['sgrg_measured']} 个可比基因整体同向率
{S['concordance']['pct']}%（ρ={n(S['concordance']['rho'])}，p={S['concordance']['spearman_p']:.2g}）。</li>
<li><b class="red">诊断/预测特异性仍不成立</b>：把模型配置固定后做跨队列迁移，
随机 82 基因面板的中位 AUC 反而高于 SGRG 面板（0.820 vs 0.740，SGRG 位于第 4 百分位）。
即所谓跨队列 AUC≈0.85 的水平，<b>任意同规模基因集都能达到</b>。</li>
<li><b>这不是又一次"什么都做不出来"</b>：复现是有<b>模块选择性</b>的——颗粒清除/翻译抑制模块
（VCP↑、ZFAND1↓、ATXN2↓、PUM1/2↓、CSDE1↓、TRIM21↑）稳定复现，
而颗粒<b>组装核心</b>（G3BP1、G3BP2、TIA1、CAPRIN1）不复现。这本身是可写进正文的结论。</li>
</ul></div>""")

A("<div class='kpi'>"
  f"<div><b>{c226['n_pairs']}</b><span>可用双生子对</span></div>"
  f"<div><b class='green'>{n(BD['GSE22619']['cramers_v'],4)}</b><span>GSE22619 批次×疾病 Cramér's V（0＝完全正交）</span></div>"
  f"<div><b>{c226['sgrg_measured']}</b><span>GSE22619 可测 SGRG 基因</span></div>"
  f"<div><b class='blue'>{S['concordance']['n_same']}/{S['concordance']['n_total']}</b><span>整体同向基因数</span></div>"
  "</div>")

# ---------------- 为什么换
A("<h2>二、为什么要换批次</h2>")
A(f"""<p>Step 2 在 GSE38713 上得到一个"看起来很好"的模型（留一 CV AUC 0.936），但三项自检一致指向批次假象：
该队列疾病状态与处理批次高度混杂（Cramér's V={n(BD['GSE38713']['cramers_v'],3)}，7 个批次中 4 个是纯 UC），
<b>只用批次标签</b>预测疾病的 AUC 就有 {n(BD['GSE38713']['batch_only_auc'],3)}；去批次后 SGRG 面板 AUC 从 0.936 塌到 0.58。</p>
<p>所以需要一个<b>批次与疾病正交</b>的独立队列。在 Gemma 收录的 153 个 UC 数据集中，我按两个硬条件筛：
① 人体结肠黏膜活检（与 GSE38713 同组织）；② 单队列压缩体积 ≤ 5 MB（本机下载通道约 11 KB/s 且约 8 分钟被掐断，单次预算≈5 MB）。
逐个体检查批次结构后，只有 <b>GSE22619</b> 满足"每批 UC 与对照严格 1:1"。</p>""")

A("<h2>三、新队列结构与批次正交证据</h2>")
A(f"""<p>GSE22619 是 <b>UC 不一致型同卵双生子</b>队列：10 对双生子，每对 1 例 UC + 1 例健康同胞，
全结肠黏膜活检。三个处理批次 <b>每一批都包含完整配对</b>，因此疾病与批次在数学上完全正交。</p>""")
A("<table><tr><th>队列</th><th class='n'>样本</th><th class='n'>UC / 对照</th>"
  "<th class='n'>Cramér's V</th><th class='n'>χ²</th><th class='n'>仅用批次预测疾病 AUC</th><th>判断</th></tr>"
  f"<tr><td>GSE22619（换批次后）</td><td class='n'>{c226['n_samples']}</td>"
  f"<td class='n'>{c226['n_uc']} / {c226['n_control']}</td>"
  f"<td class='n green'>{n(BD['GSE22619']['cramers_v'],4)}</td>"
  f"<td class='n'>{n(BD['GSE22619']['chi2'],2)} (p={BD['GSE22619']['chi2_p']:.3g})</td>"
  f"<td class='n green'>{n(BD['GSE22619']['batch_only_auc'],3)}</td><td>批次与疾病正交 ✔</td></tr>"
  f"<tr><td>GSE38713（原队列）</td><td class='n'>{c387['n_samples']}</td>"
  f"<td class='n'>{c387['n_uc_active']} / {c387['n_control']}</td>"
  f"<td class='n red'>{n(BD['GSE38713']['cramers_v'],3)}</td>"
  f"<td class='n'>{n(BD['GSE38713']['chi2'],2)} (p={BD['GSE38713']['chi2_p']:.3g})</td>"
  f"<td class='n red'>{n(BD['GSE38713']['batch_only_auc'],3)}</td><td>严重混杂 ✘</td></tr></table>")
A(img("fig_GSE22619_batch_balance.png",
      "图 1　两个队列的批次×疾病构成。左：GSE22619 每批 UC 与对照数量严格相等，Cramér's V=0；"
      "右：GSE38713 多批为纯 UC，混杂严重。"))
A("""<div class='ok'><h4>换批次这一步达成了既定目的</h4>
<p>GSE22619 的批次与疾病完全正交（V=0.0000，仅用批次标签预测疾病 AUC=0.260，低于随机水平），
因此任何在该队列中观察到的信号<b>不可能由处理批次解释</b>。这正是 GSE38713 做不到的。</p></div>""")

# ---------------- 配对差异表达
A("<h2>四、配对差异表达：唯一同时消掉批次与遗传背景的设计</h2>")
n_pair_ = c226["n_pairs"]
A(f"""<p>对每一对双生子取差值 <code>Δ = 表达(UC) − 表达(健康同胞)</code>，再做单样本 Wilcoxon 符号秩检验。
由于每对样本处在同一批次内，且为同卵双生（遗传背景一致），<b>配对差把批次效应与个体遗传背景一起消掉了</b>。</p>""")
A(f"""<div class='note'><h4>必须交代：本设计的 FDR 在数学上不可达</h4>
<p>{n_pair_} 对样本的 Wilcoxon 双侧检验最小可能 p = 2/2<sup>{n_pair_}</sup> = <b>{P['min_possible_p']:.5f}</b>；
而 {P['n_tests']} 个检验经 BH 校正后最小阈值 = 0.05/{P['n_tests']} = <b>{P['bh_min_threshold']:.5f}</b>。
二者相比，<b>最小可达 p 仍然大于 BH 阈值</b>，因此"FDR&lt;0.05 为 0 个"是统计功效的结构性上限，
<b>不能读作"没有差异"</b>。本文因此改用名义 p + 效应量 dz + 方向复现作为证据链。</p></div>""")
A(f"<p>结果：{len(DE)} 个 SGRG 基因中，名义 p&lt;0.05 共 <b>{P['n_nominal_p05']}</b> 个，"
  f"且全部 |dz| ≥ 0.3（即效应量不小）。效应量最大的几个：</p>")
A("<table><tr><th>基因</th><th class='n'>Δ(UC−同胞)</th><th class='n'>Cohen's dz</th>"
  "<th class='n'>Wilcoxon p</th><th class='n'>FDR</th><th class='n'>GSE38713 log2FC</th><th>方向</th></tr>")
for r in P["top"][:12]:
    d = next((x for x in DE if x["GeneSymbol"] == r["gene"]), {})
    fl = d.get("GSE38713_log2FC", "")
    dm = d.get("direction_match", "")
    mark = "<span class='green'>一致</span>" if dm == "+" else ("<span class='red'>相反</span>" if dm == "-" else "—")
    A(f"<tr><td><b>{r['gene']}</b></td><td class='n'>{r['mean_diff']:+.3f}</td>"
      f"<td class='n'>{r['dz']:+.3f}</td><td class='n'>{r['wilcoxon_p']:.4f}</td>"
      f"<td class='n'>{r['wilcoxon_FDR']:.3f}</td>"
      f"<td class='n'>{n(fl,3) if fl else '—'}</td><td>{mark}</td></tr>")
A("</table>")
A(img("fig_GSE22619_paired_volcano.png",
      "图 2　GSE22619 配对差异表达。横轴为同胞对内差值，纵轴为 Wilcoxon p。"
      "所有点都落在名义显著之上有限，正是 10 对样本的功效上限所致。"))

# ---------------- 方向复现
A("<h2>五、方向复现：本次换批次最扎实的正向结果</h2>")
A(f"""<p>把 GSE38713 里差异表达的基因（FDR&lt;0.05）拿到 GSE22619 的配对设计中检验<b>符号是否一致</b>：
{R['n_same_direction']}/{R['n_GSE38713_sig']} = <b>{R['pct']}%</b> 同向，二项检验
p={R['binom_p']:.2g}。扩大到全部可比基因是 {S['concordance']['n_same']}/{S['concordance']['n_total']}
= {S['concordance']['pct']}%，Spearman ρ={n(S['concordance']['rho'])}（p={S['concordance']['spearman_p']:.2g}）。</p>
<p>这是<b>在批次被完全消除、且遗传背景匹配的条件下</b>得到的复现，因此不能再用"批次假象"解释。</p>""")
A(img("fig_step3_direction_concordance.png",
      "图 3　两队列位移方向的一致性。红点为配对检验 p<0.05 的基因。绝大多数点落在第一、三象限（同向）。"))

A("<h3>复现是有模块选择性的——这才是真正的信息量</h3>")
core = ["VCP", "ZFAND1", "ATXN2", "PUM1", "PUM2", "CSDE1", "TRIM21", "RBPMS",
        "G3BP1", "G3BP2", "TIA1", "CAPRIN1", "ELAVL1", "DDX3X", "EIF2S1", "QKI", "YTHDF3", "MBNL1", "FMR1"]
mm = {r["GeneSymbol"]: r for r in DE}
A("<table><tr><th>模块</th><th>基因</th><th class='n'>Δ 配对差</th><th class='n'>dz</th>"
  "<th class='n'>Wilcoxon p</th><th class='n'>GSE38713</th><th>复现</th></tr>")
grp = {"清除 / 翻译抑制 / 泛素-蛋白酶体":
       ["VCP", "ZFAND1", "ATXN2", "PUM1", "PUM2", "CSDE1", "TRIM21", "RBPMS"],
       "颗粒组装核心": ["G3BP1", "G3BP2", "TIA1", "CAPRIN1"],
       "其他核心 SG 因子": ["ELAVL1", "DDX3X", "EIF2S1", "QKI", "YTHDF3", "MBNL1", "FMR1"]}
for gname, gl in grp.items():
    for i, g in enumerate(gl):
        r = mm.get(g)
        if not r:
            A(f"<tr><td>{gname if i == 0 else ''}</td><td>{g}</td><td colspan='5'>未测到</td></tr>"); continue
        p_ = float(r["wilcoxon_p"]); s_ = p_ < 0.05
        dm = r["direction_match"]
        rp = ("<b class='green'>同向且显著</b>" if (s_ and dm == "+") else
              ("<span class='blue'>同向未显著</span>" if dm == "+" else
               ("<span class='red'>方向相反</span>" if dm == "-" else "—")))
        A(f"<tr><td>{gname if i == 0 else ''}</td><td><b>{g}</b></td>"
          f"<td class='n'>{float(r['mean_diff_UC_minus_twin']):+.3f}</td>"
          f"<td class='n'>{float(r['cohen_dz']):+.3f}</td>"
          f"<td class='n'>{p_:.4f}</td>"
          f"<td class='n'>{n(r['GSE38713_log2FC'],3) if r['GSE38713_log2FC'] else '—'}</td>"
          f"<td>{rp}</td></tr>")
A("</table>")
A("""<div class='ok'><h4>可写进正文的结论（有方向、有模块）</h4>
<p>在批次正交、遗传背景匹配的双生子队列中，<b>复现的是"应激颗粒清除与翻译抑制"这条轴</b>：
<b>VCP 上调</b>（Δ=+0.348，dz=+1.08，p=0.0098）、<b>ZFAND1 下调</b>（p=0.0020）、
<b>ATXN2 下调</b>（p=0.027）、<b>PUM1/PUM2 下调</b>（p=0.014 / 0.049）、<b>CSDE1 下调</b>（p=0.049）、
<b>TRIM21 上调</b>（p=0.020）、<b>RBPMS 上调</b>（p=0.027）。</p>
<p>而<b>颗粒组装核心 G3BP1、G3BP2、TIA1、CAPRIN1 均未复现</b>（G3BP1 Δ=−0.087, p=0.28；
G3BP2 方向相反）。VCP↑ 与 ZFAND1↓ 同时出现，恰好对应"含异常蛋白的应激颗粒清除通路"被激活——
这比笼统地说"SGRG 整体上调"精确得多，也更抗审稿质疑。</p></div>""")

# ---------------- 诊断模型
A("<h2>六、诊断模型：换批次之后，特异性问题反而更清楚了</h2>")
A("<h3>6.1 跨队列外部验证</h3>")
A("<table><tr><th>训练 → 验证</th><th>最优配置</th><th class='n'>非零基因</th>"
  "<th class='n'>训练内 CV AUC</th><th class='n'>外部 AUC [95%CI]</th>"
  "<th class='n'>灵敏度 / 特异度</th></tr>")
labels = {"38713_allUC→22619": "GSE38713(全部 UC) → GSE22619",
          "38713_inflamed→22619": "GSE38713(活动期炎症区) → GSE22619",
          "22619→38713_inflamed": "GSE22619 → GSE38713(活动期炎症区)"}
for k, v in S["cross_cohort"].items():
    e = v["ext"]
    A(f"<tr><td>{labels.get(k,k)}</td><td>{v['penalty'].upper()} · C={v['C']}</td>"
      f"<td class='n'>{v['n_nonzero']}</td><td class='n'>{n(v['cv_auc'])}</td>"
      f"<td class='n'><b>{n(e['auc'])}</b> [{n(e['auc_lo'])}, {n(e['auc_hi'])}]</td>"
      f"<td class='n'>{pct(e['sens'])} / {pct(e['spec'])}</td></tr>")
A("</table>")
A(img("fig_step3_crosscohort_ROC.png", "图 4　跨队列 ROC。"))
A(f"<p>表面上「炎症区训练 → 双生子队列验证」拿到外部 AUC "
  f"{n(S['cross_cohort']['38713_inflamed→22619']['ext']['auc'])}，看起来可迁移。但下面的检验否掉了它。</p>")

A("<h3>6.2 面板特异性：<span class='red'>关键的反转</span></h3>")
A("""<div class='warn'><h4>上一版的"面板特异性检验"在本场景下是无效的</h4>
<p>n=20、特征 82 个时，任何基因集的<b>同队列交叉验证 AUC 都会逼近 1</b>——
本次实测随机 82 基因面板的中位 CV AUC 就是 0.955，SGRG 面板 0.920 反而低于随机中位数。
因此"同队列 CV AUC"完全不能用来论证 SGRG 的特异性。真正有判别力的口径是<b>跨队列迁移</b>。</p></div>""")
A("<table><tr><th>检验</th><th>模型配置</th><th class='n'>SGRG 面板迁移 AUC</th>"
  "<th class='n'>随机面板中位</th><th class='n'>随机面板 p95</th><th class='n'>SGRG 百分位</th>"
  "<th class='n'>经验 p</th></tr>")
for k, v in S["cross_cohort_panel_specificity"].items():
    A(f"<tr><td>{k.split('·')[1].strip() if '·' in k else k}</td><td>{v['model']}</td>"
      f"<td class='n'><b>{n(v['sgrg_transfer_auc'])}</b></td>"
      f"<td class='n'>{n(v['rand_median'])}</td><td class='n'>{n(v['rand_p95'])}</td>"
      f"<td class='n'>{v['percentile']:.0f}</td><td class='n'>{n(v['empirical_p'],3)}</td></tr>")
A("</table>")
A(img("fig_step3_cross_panel_specificity.png",
      "图 5　跨队列面板特异性。红竖线为 SGRG 面板，灰虚线为随机面板中位数。"))
A("""<div class='warn'><h4>结论：跨队列 AUC 不是 SGRG 特有的</h4>
<ul>
<li>固定配置（L2, C=0.01）下，SGRG 面板迁移 AUC <b>低于</b>随机面板中位数（0.740 vs 0.820，第 4 百分位）。</li>
<li>即使改用对 SGRG 最有利的配置（L1, C=0.3），SGRG 也只是到第 72 百分位，经验 p=0.28，<b>未达显著</b>。</li>
<li>换句话说：<b>任意 82 个基因就能达到 AUC≈0.85 的跨队列分离水平</b>，因为 UC 黏膜炎症本身是一个
全转录组规模的广泛位移，任何足够大的基因集都会捕获这条主轴。</li>
</ul>
<p>因此"用 SGRG 面板做 UC 诊断模型"这一主张，<b>不能由本数据支持</b>。</p></div>""")

# ---------------- 合并判断
A("<h2>七、与前一步结论合并后的最终判断</h2>")
A("<table><tr><th>主张</th><th>证据</th><th>结论</th></tr>"
  "<tr><td>SGRG 基因在 UC 黏膜中差异表达</td><td>GSE38713 44 个 FDR&lt;0.05；GSE22619 19 个名义 p&lt;0.05</td>"
  "<td class='green'>支持（但效应量普遍偏小）</td></tr>"
  f"<tr><td>差异方向可复现</td>"
  f"<td>{R['pct']}% 同向，二项 p={R['binom_p']:.2g}；整体 ρ={n(S['concordance']['rho'])}</td>"
  "<td class='green'><b>支持</b>（且在批次正交队列中）</td></tr>"
  "<tr><td>复现集中于「清除/翻译抑制」模块而非组装核心</td>"
  "<td>VCP↑ ZFAND1↓ ATXN2↓ PUM1/2↓ CSDE1↓ TRIM21↑ 复现；G3BP1/2、TIA1、CAPRIN1 未复现</td>"
  "<td class='green'>支持（新的、更精确的结论）</td></tr>"
  "<tr><td>SGRG 面板可作为 UC 诊断/分型模型</td>"
  f"<td>固定配置下 SGRG 迁移 AUC 仅第 {S['cross_cohort_panel_specificity']['fixed_L2C0.01 · 38713(inflamed)→22619']['percentile']:.0f} 百分位；"
  "最优配置下 p=0.28 不显著；GSE38713 内 0.936 由批次驱动</td>"
  "<td class='red'><b>不支持</b></td></tr>"
  "<tr><td>GSE38713 单队列 AUC 0.936 反映生物学</td><td>批次×疾病 V=0.712，去批次后塌至 0.58</td>"
  "<td class='red'>不支持</td></tr></table>")

A("""<h3>正文写作建议</h3>
<p><b>可以写</b>：SGRG 在 UC 黏膜中存在方向一致、可跨队列复现的表达改变，且改变集中于
应激颗粒清除与翻译抑制通路（VCP 上调、ZFAND1 下调为最稳健的一对）。建议以这 7 个基因
（VCP、ZFAND1、ATXN2、PUM1、PUM2、CSDE1、TRIM21）作为候选面板进入后续功能验证。</p>
<p><b>不能写</b>：SGRG 面板具有 UC 诊断价值；SGRG 整体上调；任何单独依赖 GSE38713 的 AUC 数字。</p>
<p><b>建议降级处理</b>：把转录组部分定位为"支持性证据"，用配对双生子队列的方向复现作为主要卖点，
把因果推断的权重交给 Part C 的孟德尔随机化（不受批次与横断面设计限制）。</p>""")

# ---------------- 局限
A("<h2>八、局限与必须声明的技术问题</h2>")
A(f"""<ul>
<li><b>矩阵被截断（重要）</b>：GSE22619 完整矩阵 4.3 MB，本机下载通道在 2.39 MB（53%）处被中断，
仅保留<b>按行截断的前缀</b>（28,742 行基因 × 20 样本列结构完整）。因此该队列只覆盖
{c226['sgrg_measured']}/104 个 SGRG 基因，且缺失基因集中在阵列后段探针（前缀非随机），
可能带来弱选择偏倚。头部基因（G3BP1、G3BP2、EIF2S1、VCP、ATXN2、FMR1、CSDE1、QKI 等）均在覆盖范围内。</li>
<li><b>配对设计功效上限</b>：仅 10 对样本，Wilcoxon 最小可达 p={P['min_possible_p']:.5f}，
大于 BH 阈值 {P['bh_min_threshold']:.5f}，所以 FDR&lt;0.05 在数学上不可达。所有配对结论只能用名义 p + 效应量。</li>
<li><b>样本量小</b>：20 例样本使同队列 CV AUC 失真（见 6.2），任何在该队列内做的交叉验证数字都不可用。</li>
<li><b>平台差异</b>：两队列均基于 Affymetrix 平台但批次/年份不同，跨队列评分为"方向级"比较，
不构成绝对阈值可迁移的临床模型。</li>
<li><b>组织类型差异</b>：GSE22619 为全黏膜活检（含全层黏膜细胞），GSE38713 为活检黏膜，细胞组成差异未做反卷积校正。</li>
</ul>""")

A(f"""<footer>数据来源：Gemma（gemma.msl.ubc.ca）镜像的 GSE38713（id=5873）、GSE22619（id=6731）；
SGRG 基因集为 GO:0010494 ∪ GO:0035061（人类，经 HGNC 官方接口校验，104 个批准符号）。
分组与批次一律取自 Gemma 实验设计的权威因子（bioMaterialAssignments），非文本推断。<br>
复现脚本：<code>scripts/partA_step3_paired_validation.py</code>　共享模块：<code>scripts/sgrg_common.py</code><br>
矩阵完整性状态：<b>{S['matrix_status_of_GSE22619']}</b>（截断前缀，已在全文标注）。</footer>""")

A("</div></body></html>")

p = os.path.join(OUT, "PartA_Step3_换批次验证报告.html")
open(p, "w", encoding="utf-8").write("\n".join(H))
print("[write]", p, os.path.getsize(p), "bytes")
