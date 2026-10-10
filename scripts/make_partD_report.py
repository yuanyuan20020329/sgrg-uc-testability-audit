# -*- coding: utf-8 -*-
"""Part D 报告生成器（自包含 HTML，图片 base64 内嵌）"""
import base64
import json
import os
import csv

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT = os.path.join(BASE, "results")


def jload(n, d=None):
    p = os.path.join(OUT, n)
    if not os.path.exists(p):
        return d
    return json.load(open(p, encoding="utf-8"))


def tload(n):
    p = os.path.join(OUT, n)
    if not os.path.exists(p):
        return []
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def img(name, cap=""):
    if not os.path.splitext(name)[1]:
        name = name + ".png"
    p = os.path.join(OUT, name)
    if not os.path.exists(p):
        print("[warn] 缺图", p)
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
        return "&lt;1e-300"
    return f"{x:.3f}" if x >= 1e-3 else f"{x:.2e}"


def tbl(headers, rows, cls=""):
    h = "".join(f"<th>{x}</th>" for x in headers)
    b = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f"<table class='{cls}'><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table>"


ov = jload("partD_ds1_overview.json", {})
det = ov.get("sgrg_detection", {})
det_sorted = sorted(det.items(), key=lambda kv: -kv[1]) if det else []
n_ge1 = sum(1 for _, v in det_sorted if v >= 0.01)
n_ge5 = sum(1 for _, v in det_sorted if v >= 0.05)

pb = tload("partD_ds1_pseudobulk_score.tsv")
tests = tload("partD_ds1_pseudobulk_tests.tsv")
adj = {r["cell_type"]: r for r in tload("partD_ds1_depth_adjusted_tests.tsv")}
pg = tload("partD_ds1_sgrg_gene_pseudobulk.tsv")
diag = tload("partD_ds1_depth_diagnostics.tsv")
coup = tload("partD_ds1_state_coupling.tsv")
q5 = {r["cell_type"] + "|" + r["quintile"]: r for r in tload("partD_ds1_state_quintile.tsv")}
ds2 = jload("partD_ds2_myeloid_summary.json", {})
ds3 = jload("partD_ds3_largeintestine_summary.json", {})
t2 = tload("partD_ds2_myeloid_uc_vs_hc_tests.tsv")
t3 = tload("partD_ds3_largeintestine_uc_vs_hc_tests.tsv")
ds4 = jload("partD_ds4_tearle_summary.json", None)
t4 = tload("partD_ds4_tearle_uc_vs_hc_tests.tsv")
t4r = tload("partD_ds4_tearle_ds1k_uc_vs_hc_tests.tsv")
# 负对照（随机基因集，两个口径）
rnd_raw = jload("partD_ds4_tearle_rand_summary.json", None)
rnd_ds = jload("partD_ds4_tearle_ds1k_rand_summary.json", None)

baudit = tload("partD_batch_audit_summary.tsv")
prov = jload("partD_ds4_prov_summary.json", None)
tpv = tload("partD_ds4_prov_contrast_tests.tsv")
tin = tload("partD_ds4_prov_inflamed_tests.tsv")
spv = tload("partD_ds4_prov_geneset_summary.tsv")

SHORT = {"colon epithelial cell": "上皮", "myeloid cell": "髓系", "T cell of anorectum": "T 细胞",
         "plasma cell": "浆细胞", "stromal cell of lamina propria of colon": "基质"}
ORDER = ["colon epithelial cell", "myeloid cell", "T cell of anorectum", "plasma cell",
         "stromal cell of lamina propria of colon"]

# ds4 展示顺序：按「原始 → 抽稀」的收缩幅度降序，取前 8
SHORT_D4 = {
    "CD4-positive, alpha-beta T cell": "CD4 T", "CD8-positive, alpha-beta T cell": "CD8 T",
    "CD8-positive, alpha-beta memory T cell": "CD8 记忆 T", "T follicular helper cell": "Tfh",
    "T-helper 17 cell": "Th17", "capillary endothelial cell": "毛细血管内皮",
    "conventional dendritic cell": "常规 DC", "endothelial cell of artery": "动脉内皮",
    "epithelial cell of large intestine": "大肠上皮", "fibroblast": "成纤维细胞",
    "gamma-delta T cell": "γδ T", "group 3 innate lymphoid cell": "ILC3",
    "intestinal crypt stem cell of large intestine": "隐窝干细胞",
    "intestine goblet cell": "杯状细胞", "macrophage": "巨噬细胞", "mast cell": "肥大细胞",
    "memory B cell": "记忆 B", "myofibroblast cell": "肌成纤维细胞", "naive B cell": "初始 B",
    "pericyte": "周细胞", "plasma cell": "浆细胞", "regulatory T cell": "Treg",
    "tuft cell of colon": "tuft 细胞", "type I enteric glial cell": "肠胶质细胞",
    "vein endothelial cell": "静脉内皮"}
_d4o = {r["cell_type"]: float(r["diff"]) for r in (t4 or []) if r["score"] == "SGRG"}
_d4d = {r["cell_type"]: float(r["diff"]) for r in (t4r or []) if r["score"] == "SGRG"}
ORDER_D4 = sorted(_d4o, key=lambda c: -abs(_d4o[c] - _d4d.get(c, 0.0)))

H = []
A = H.append
A("""<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">
<title>Part D｜单细胞方案启动：SGRG 在 UC 黏膜中的细胞定位与炎症状态</title>
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
table{border-collapse:collapse;width:100%;font-size:13.5px;margin:12px 0;background:#fff}
th,td{border:1px solid var(--line);padding:7px 9px;text-align:left;vertical-align:top}
th{background:#eef3f9;font-weight:600}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
figure{margin:18px 0;background:#fff;border:1px solid var(--line);border-radius:10px;padding:12px}
img{width:100%;display:block}
figcaption{color:var(--muted);font-size:13px;margin-top:8px;line-height:1.6}
code{background:#eef1f5;padding:1px 5px;border-radius:4px;font-size:13px}
.kv{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:10px;margin:14px 0}
.kv div{background:#fff;border:1px solid var(--line);border-radius:8px;padding:10px 12px;font-size:13.5px}
.kv b{display:block;font-size:19px;color:var(--blue);font-variant-numeric:tabular-nums}
ul,ol{margin:10px 0 10px 22px}li{margin:5px 0}
.tag{display:inline-block;font-size:12px;padding:1px 8px;border-radius:20px;background:#eef3f9;color:var(--blue);margin-right:6px}
.small{font-size:13px;color:var(--muted)}
</style></head><body><div class="wrap">""")

# ================= 标题 =================
A(f"""
<h1>Part D｜单细胞方案启动：应激颗粒基因程序在溃疡性结肠炎黏膜中的细胞定位与炎症状态</h1>
<div class="sub">可行性核验 · 数据通道 · 首批结果 · 分析设计 ｜ 2026-10-09 ｜ 承接《项目结果审计_结论与出路》的「出路 C」<br>
本报告全部数字可由 <code>results/partD_*.tsv|json</code> 与 <code>scripts/partD_*.py</code> 复现</div>
""")

# ================= 0 一句话 =================
A(f"""
<div class="card">
<h2 style="border:0;margin-top:0;padding-top:0">一句话结论</h2>
<p><b>单细胞这条路本机走得通</b>：CELLxGENE Discover 与 NCBI/GEO 两条通道实测可达，已取回并跑完
4 套 UC 黏膜 scRNA-seq（46,700 / 11,720 / 87,621 / <b>320,979</b> 细胞，合计约 46 万细胞）。</p>
<p><b>四套数据集、三种独立的深度控制方法，结论完全一致</b>：SGRG 面板在 UC 黏膜中
<b>没有任何一种细胞类型</b>出现程序性升高；200 个表达分层匹配的随机基因集走完全相同的流程，
未校正口径下与 SGRG 无法区分（经验 p=0.49–0.86），抽稀后伪影同步归零。</p>
<p><b>但本轮把这条链又往上游推了一层，并且推翻了上一轮结论的解读方式</b>：审计四套数据的
<em>对照来源</em>发现——</p>
<ul>
<li><b>ds4 的「正常」组里 94.4% 的细胞（177,982 / 188,501）根本不是这项研究自己的对照</b>
（<code>is_primary_data=False</code> 的复用切除标本，批次、化学版本、深度完全不同）。
按数据来源重算，深度差从 <b>4.3× 塌缩到 1.04×</b>——
上一轮归因于「疾病」的深度混杂，其实是<b>对照组成</b>造成的；</li>
<li><b>ds2/ds3 的对照池 98.5% / 80.7% 来自完全没有病例的研究</b>（UC 只出现在 2–3/25 与 2–3/11 个研究里），
我们此前的 ds2/ds3「UC vs 正常」检验同样是跨研究对照；</li>
<li>换成<b>设计上受控的口径</b>重算：UC vs 本研究自己的对照只剩 <b>2 例对照</b>
（连 NF-κB / IFN-γ 通路都测不出差异 → <b>零检验力</b>）；
唯一真正干净的「同 donor 炎症 vs 非炎症」配对只有 <b>4 例</b>
（配对设计的最小可达 p=0.125 → 同样无检验力）。</li>
</ul>
<p><b>所以 Part D 的最终判定变了</b>：不是「面板在 UC 未被激活」，而是
<b>现有公共单细胞资源在结构上无法对这个问题做出有检验力的、批次受控的检验</b>——
四种口径要么混杂不可识别、要么检验力为零。这个结构性结论本身，比任何一个「阴性 p 值」都更值得写。</p>
</div>
""")

# ================= 1 定位 =================
A("""
<h2>一、这条线要回答什么，为什么不是 Part A 的重做</h2>
<div class="card">
<p>审计报告给出的「出路 C」把 bulk 面板的两个死结（<b>批次混杂</b>、<b>面板级富集不成立</b>）作为前提，
换一个能被现有数据回答的问题：</p>
<p style="font-size:17px;border-left:4px solid #0b5cab;padding-left:14px;margin:16px 0">
<b>应激颗粒（SGRG）基因程序在 UC 黏膜里由哪一类细胞承担、在什么样的炎症状态下被激活？</b></p>
<p>它与 Part A 的区别在于提问的层次：Part A 问「整个组织里面板的平均表达是否升高」（被批次与细胞组成主导），
Part D 问「<b>在细胞类型内部</b>，面板程序的水平与状态如何变化」——后者天然把细胞组成这个混杂因素拿掉了，
是 bulk 数据做不到的事。</p>
<div class="warn"><b>先立好判据（避免重演 Part A）。</b>本方案在跑数之前先固定三条失败条件，出现任意一条即判定「面板级阳性不成立」：
<ol>
<li>若「UC vs 正常」的差异在<b>按 donor 聚合（pseudobulk）</b>后消失，则只算细胞数伪重复；</li>
<li>若差异在<b>控制测序深度</b>后消失，则算技术混杂；</li>
<li>若只在某一个数据集显著、在独立队列不重复，则算数据集特异。</li>
</ol></div>
</div>
""")

# ================= 2 通道 =================
A("""
<h2>二、沙箱数据通道核验（这一步决定方案能不能落地）</h2>
<div class="card">
<p>本轮实测的可用性（与项目此前记录不同的一处已高亮）：</p>
""")
rows = [
    ["CELLxGENE Discover API", "<span class='tag'>可达</span>", "api.cellxgene.cziscience.com", "397 个公开集合；含标准化的 cell_type / disease / tissue 注释与 h5ad 资产直链"],
    ["CELLxGENE 数据资产 S3", "<span class='tag'>可达</span>", "datasets.cellxgene.cziscience.com", "h5ad 直下，支持 HTTP Range 断点续传"],
    ["NCBI eutils + GEO FTP", "<span class='tag'>本次可达</span>", "eutils.ncbi.nlm.nih.gov / ftp.ncbi.nlm.nih.gov", "<b>与此前记忆相反</b>：现在通。UC 单细胞 GEO 系列检索命中 101 条（含 GSE282124/312415 等空间转录组）"],
    ["EBI（SCEA / BioStudies / HCA）", "<span class='tag'>可达</span>", "www.ebi.ac.uk / ftp.ebi.ac.uk", "作为备份通道"],
    ["MSigDB Hallmark GMT", "<span class='tag'>可达</span>", "data.broadinstitute.org", "直接下 GMT，免注册"],
    ["QuickGO 注释", "<span class='tag'>可达</span>", "www.ebi.ac.uk/QuickGO", "把 SGRG 面板拆成胞质/核两亚组"],
    ["Enrichr 基因集库", "<span class='tag'>可达</span>", "maayanlab.cloud", "备用"],
    ["Ensembl REST / UCSC / Zenodo", "<span class='tag'>不可达</span>", "—", "无 liftOver，坐标一律按 GRCh38 原样使用"],
]
A(tbl(["通道", "状态", "端点", "说明"], rows))
A("""<div class="note"><b>对方案的直接影响：</b>「NCBI 反而不通」这个长期假设被推翻后，GEO 直接补数据成为可能——
后面若需要某个特定队列（例如 Smillie 2019 或 Kinchen 2018），不必再绕道 Gemma。</div></div>""")

# ================= 3 数据集 =================
A("""
<h2>三、已取回的数据集</h2>
<div class="card">
<p>选择标准：①结肠/直肠黏膜；②同时含 UC 与正常对照；③有标准化的细胞类型注释；④原始 counts 可得（可自行归一化）。</p>
""")
ds_rows = [
    ["ds1（主）", "Garrido-Trigo 2023 <i>Nat Commun</i>", "GSE／CELLxGENE", "46,700", "32,354", "5", "UC 6 / CD 6 / 正常 6", "183.7 MB", "<b>已完成</b>"],
    ["ds2", "Oliver 2024 <i>Nature</i>（髓系分库）", "CELLxGENE", "50,570（结肠 11,720）", "35,461", "12", "223 donors（UC 少）", "376 MB", "<b>已完成</b>"],
    ["ds3", "Oliver 2024 <i>Nature</i>（大肠分库）", "CELLxGENE", "96,675（结肠 87,621）", "35,461", "9（上皮亚型）", "UC 3–4 donors 有效", "861 MB", "<b>已完成</b>"],
    ["ds4", "Tearle 2026 <i>Nat Commun</i>（PSC + UC 黏膜）", "CELLxGENE", "320,979", "37,490", "34", "UC 8 / 正常 7 / PSC 7（donor）", "2.7 GB", "<b>已完成</b>（含深度匹配重跑）"],
]
A(tbl(["编号", "来源", "通道", "细胞数", "基因数", "细胞类型", "分组", "大小", "状态"], ds_rows))
A(f"""<div class="note"><b>矩阵状态必须先判别</b>：CELLxGENE 的 <code>X</code> 有的存原始 counts、有的存已 log 标准化矩阵
（Oliver 系列即如此，<code>raw/X</code> 里才是 counts）。脚本用「整数占比」自动判定并优先取 <code>raw/X</code>；
若不做这步，会把已标准化的矩阵再归一化一次，结果整体失真。</div>
</div>""")

# ================= 4 方法 =================
A(f"""
<h2>四、方法</h2>
<div class="card">
<ol>
<li><b>归一化</b>：统一取原始 counts → 每细胞 CP10K → log1p。</li>
<li><b>面板</b>：GO:0010494（胞质应激颗粒）+ GO:0035061（核应激颗粒）的 HGNC 批准符号共 {ov.get('sgrg_panel', 104)} 个；
另按 GO 注释拆出<b>胞质亚组</b>与<b>核亚组</b>分别打分（二者结论一致，见结果表）。</li>
<li><b>打分</b>：scanpy <code>score_genes</code> 口径——观测集均值减去按表达量分箱（25 箱）抽样的对照基因集均值；
对照集重复 100–200 次给出零分布，用于换算 z 值。<b>目的是扣掉「表达量本身」带来的偏置</b>。</li>
<li><b>统计单元 = donor，不是细胞</b>：所有组间比较先按 <code>donor × cell_type</code> 聚合成 pseudobulk
（每组 ≥20 细胞才纳入），再做置换检验（20,000 次），BH 校正。<b>这是防止「细胞数伪重复」的关键一步</b>。</li>
<li><b>深度校正</b>：两种独立口径——
  ① <b>残差法</b>：在细胞类型内、pseudobulk 水平上对 log10(nUMI) 做线性回归取残差再检验；
  ② <b>深度匹配抽稀</b>（ds4）：把每个细胞的 counts 用二项抽稀（binomial thinning）降到同一深度
  （1,000 UMI），不足 1,000 的细胞弃去，然后在完全相同的深度上重新打分与检验。
  ②不依赖任何模型假设，是①在「深度与疾病共线」时的唯一可靠替代。</li>
<li><b>状态分析</b>：细胞类型内部按炎症复合评分（TNFA-NFKB / INFLAMMATORY / IFN-α / IFN-γ / IL6-JAK-STAT3
五个 Hallmark 的 z 均值）分五分位，比较深度校正后的面板评分。</li>
</ol>
</div>
""")

# ================= 5 结果 =================
A("<h2>五、结果</h2>")

# 5.1 覆盖
A(f"""
<div class="card"><h3>5.1 面板在单细胞里的可测性远好于预期</h3>
<div class="kv">
<div><b>{len(ov.get('sgrg_present', []))}/{ov.get('sgrg_panel', 104)}</b>面板基因可检出（仅缺 L1RE1）</div>
<div><b>{n_ge1}</b>个基因在 ≥1% 细胞中检出</div>
<div><b>{n_ge5}</b>个基因在 ≥5% 细胞中检出</div>
<div><b>64.3%</b>最高检出率（PABPC1）</div>
</div>
<p class="small">对照组（Part A 的 bulk 芯片）里可测到 97/104。单细胞虽然每个细胞只测到约 1,280 个基因，
但面板基因多为高表达 RNA 结合蛋白，检出率并不吃亏：PABPC1 64.3%、ZFP36 56.7%、CIRBP 54.7%、HNRNPK 54.4%。
低检出端的 PABPC1L2A / LIN28A / DDX25 等（&lt;0.1%）与生物学一致（生殖细胞/胎儿特异），不构成障碍。</p>
</div>""")

# 5.2 组成
comp = ov.get("composition", {})
A("""
<div class="card"><h3>5.2 UC 黏膜的细胞组成本身发生了剧烈变化</h3>
""")
A(img("fig_D1_composition", "图 D1｜发现队列（ds1，46,700 细胞）中三种疾病状态的细胞组成。UC 黏膜中浆细胞占比由正常的 22% 升至 57%，上皮细胞由 44% 降至 7%。"))
A("""<p>这是本方案最重要的一条「前置事实」：<b>UC 与正常之间最大的差异是细胞组成，不是某个细胞类型内部的表达</b>。
任何在整块组织上做的「面板平均表达」比较，都会被这个组成位移主导——Part A 里看到的「SGRG 在 UC 中变化」
至少有相当一部分来自这里，而不是细胞内的应激颗粒程序本身。把这层剥掉，正是单细胞方案的价值所在。</p>
</div>""")

# 5.3 基线
base = ov.get("baseline_by_celltype", [])
if base:
    rows = [[SHORT.get(r["cell_type"], r["cell_type"]), nf(r["mean"], 4), nf(r["median"], 4), f"{int(r['count']):,}"]
            for r in base]
    A("<div class='card'><h3>5.3 面板程序在哪种细胞里最强（基线定位）</h3>" +
      tbl(["细胞类型", "面板评分均值", "中位数", "细胞数"], rows) +
      """<p>基质细胞与 T 细胞的基线下限最高、浆细胞最低。但<b>不能直接读这个排名</b>——下一节说明原因。</p></div>""")

# 5.4 深度混杂
A("""
<div class="card"><h3>5.4 关键方法学发现：面板评分与测序深度系统性相关</h3>
""")
A(img("fig_D3_depth_confounding", "图 D3｜每个细胞类型内，面板评分（纵轴）对测序深度（横轴）的关系；红点为线性拟合。相关系数 ρ 均在 0.33–0.93 之间。"))
if diag:
    rows = [[SHORT.get(r["cell_type"], r["cell_type"]), nf(r["rho_score_nUMI"], 3),
             nf(r["rho_z_nUMI"], 3), f"{float(r['median_nUMI_UC']):,.0f}", f"{float(r['median_nUMI_normal']):,.0f}"]
            for r in diag]
    A(tbl(["细胞类型", "ρ(评分, 深度)", "ρ(z 值, 深度)", "UC 深度中位", "正常深度中位"], rows))
A("""<p>两件事同时成立：①面板评分随深度单调上升（髓系 ρ=0.77、基质 ρ=0.61）；②UC 与正常的细胞深度本身不同
（髓系 UC 3,645 vs 正常 2,214；T 细胞 UC 2,864 vs 正常 3,440）。<b>两者叠加，足以在没有任何生物学差异的情况下
凭空造出「显著差异」。</b></p>
<div class="warn"><b>这条要写进论文的方法学部分。</b>基因集打分（score_genes / AUCell / 模块评分）在单细胞里被大量使用，
但很少有人报告它与 nUMI 的相关性，更少有人在跨条件比较时校正它。本项目的三个 Part 各暴露了一种同类问题：
Part A 的批次-疾病混杂、Part C 的距离剪枝伪重复、Part D 的深度-基因集评分混杂。</div>
</div>""")

# 5.5 UC vs 正常
A("<div class='card'><h3>5.5 UC vs 正常：控制深度后，没有任何细胞类型升高</h3>")
A(img("fig_D2_pseudobulk_score", "图 D2｜pseudobulk（donor 级，6 例/组）面板评分。左：原始评分；右：控制测序深度后的残差。每个点为一个 donor。"))
if tests:
    rows = []
    for c in ORDER:
        t = [r for r in tests if r["cell_type"] == c and r["score"] == "SGRG" and r["contrast"] == "UC_vs_normal"]
        a = adj.get(c)
        if not t:
            continue
        t = t[0]
        p0 = float(t["p_perm"])
        mark0 = "<b>%s</b>" % psci(p0) if p0 < 0.05 else psci(p0)
        p1 = float(a["p_UC"]) if a else float("nan")
        rows.append([SHORT.get(c, c), f"{float(t['mean2']):.4f}", f"{float(t['mean1']):.4f}",
                     nf(t["cohen_d"], 2), mark0, nf(a["diff_UC"], 4) if a else "—",
                     psci(p1) if a else "—"])
    A(tbl(["细胞类型", "正常", "UC", "Cohen's d", "p（原始）", "UC−正常（残差）", "p（深度校正）"], rows,
          cls=""))
A("""<p>原始口径下 T 细胞显示 UC 下降（p=0.020，d=−1.46），看似是一个「UC 中应激颗粒程序下调」的发现。
但控制深度后差异从 −0.030 缩到 −0.010、p 升到 0.189 —— <b>它来自 UC 的 T 细胞测序更浅</b>，不是生物学。</p>
<div class="ok"><b>判定：面板级阳性不成立（三条失败条件中的第 2 条被触发）。</b>
在 UC 黏膜中，没有任何一类细胞的 SGRG 面板程序出现升高；"下调"的说法同样不成立。</div>
</div>""")

# 5.6 逐基因
A("<div class='card'><h3>5.6 逐基因：趋势存在，但无一通过多重校正</h3>")
A(img("fig_D5_gene_effect_heatmap", "图 D5｜各细胞类型内 SGRG 基因的 pseudobulk 效应量（Cohen's d，UC vs 正常）。红=UC 升高，蓝=UC 降低。"))
if pg:
    rows = []
    for c in ORDER:
        sub = [r for r in pg if r["cell_type"] == c and r.get("d_UC_vs_normal") not in (None, "")]
        try:
            sub = sorted(sub, key=lambda r: -abs(float(r["d_UC_vs_normal"])))[:5]
        except Exception:
            continue
        for r in sub:
            rows.append([SHORT.get(c, c), r["gene"], nf(r["hc"], 3), nf(r["uc"], 3),
                         nf(r["d_UC_vs_normal"], 2), psci(r["p_UC_vs_normal"]), psci(r["q_UC_vs_normal"])])
    A(tbl(["细胞类型", "基因", "正常", "UC", "Cohen's d", "p", "q(BH)"], rows))
A("""<p>可读出的模式有三条，但都<b>未通过 FDR</b>（最小 q≈0.069）：</p>
<ul>
<li><b>CIRBP 在髓系、浆细胞、基质三类细胞中一致下降</b>（d=−3.2 / −2.5 / −3.2）：CIRBP 是冷诱导的 RNA 结合蛋白，
也是应激颗粒组分，这条一致性值得在更大队列里追。</li>
<li><b>ZFP36（TTP）在 T 细胞与浆细胞中下降</b>（d=−2.4 / −3.2）：与炎症下 mRNA 稳定性调控的改变方向一致。</li>
<li><b>髓系中 C9orf72、TRIM25、MBNL1 上升</b>（d=3.2 / 3.8 / 2.8）：这是唯一一处「UC 中升高」的信号，
且集中在固有免疫细胞，与「髓系细胞在炎症中被激活」的先验一致。</li>
</ul>
<p class="small">注：以上为发现队列（n=6/组）的名义趋势，<b>不作为结论</b>；在 ds4（320,979 细胞）中的复核见 5.9——
鉴于 ds4 的面板级信号在深度匹配后整体消失，这些基因级趋势也应当作深度混杂的同源产物对待。</p>
</div>""")

# 5.7 状态耦合
A("<div class='card'><h3>5.7 在哪种炎症状态下形成：耦合存在，但在 UC 中被削弱</h3>")
A(img("fig_D4_inflammation_quintile", "图 D4｜细胞类型内，按炎症复合评分五分位分箱后的深度校正面板评分（点=均值，误差棒=SEM）。"))
if coup:
    ct_order = ORDER
    rows = []
    for c in ct_order:
        r1 = next((r for r in coup if r["cell_type"] == c and r["disease"] == "normal"), None)
        r2 = next((r for r in coup if r["cell_type"] == c and r["disease"] == "ulcerative colitis"), None)
        g = lambda q, k: float(q5.get("%s|%d" % (c, q), {}).get(k, float("nan")))
        rows.append([SHORT.get(c, c),
                     nf(r1["rho_depth_adj"], 3) if r1 else "—",
                     nf(r2["rho_depth_adj"], 3) if r2 else "—",
                     nf(g(5, "sgrg_mean") - g(1, "sgrg_mean"), 4)])
    A(tbl(["细胞类型", "ρ（正常）", "ρ（UC）", "Q5−Q1"], rows))
A("""<p>把深度差异扣掉后，面板程序与炎症状态之间<b>确实存在正相关</b>，但强度依细胞类型而变，且
<b>在 UC 中普遍低于正常</b>（髓系 0.37→0.12，基质 0.32→0.18）。五分位趋势单调上升，跨度 0.026–0.131。</p>
<div class="note">这一条需要谨慎解读：UC 样本整体已处于高炎症区，炎症评分的组内动态范围被压缩，
耦合系数下降可能只是「天花板效应」而非调控关系改变。<b>当前证据不足以区分两者</b>，需要在
更大的队列里用状态匹配（state-matched）的方式重做。</div>
</div>""")

# 5.8 两个补充数据集
A(f"""
<div class="card"><h3>5.8 两个补充数据集：细粒度细胞类型可用，但 UC 样本太薄</h3>
<p>为了把「哪种细胞」问到亚型级别，另外分析了两个数据集：</p>
<ul>
<li><b>ds3（上皮亚型，96,675 细胞）</b>给出清晰的上皮内梯度：肠内分泌细胞 0.183 ≈ 隐窝干细胞 0.171 ＞
增殖性过渡放大细胞 0.164 ＞ 杯状细胞 0.132 ＞ BEST4+ 肠上皮 0.123 ＞ 成熟结肠细胞 0.104。
即<b>越处于增殖/未分化状态的细胞，SGRG 程序越强</b>——与应激颗粒参与翻译重编程的生物学一致。</li>
<li><b>ds2（髓系亚型，12 类）</b>：浆细胞样树突状细胞 0.211 ＞ 常规 DC 0.191 ＞ 巨噬细胞 0.182 ＞
单核细胞 0.181 ＞ TREM2+ 巨噬 0.172 ＞ 肥大细胞 0.150。</li>
</ul>
<p>但两者的 UC 细胞数都太少（ds3 仅 3–4 个 UC donor 达到每类型 ≥20 细胞；ds2 的 UC 仅约 800 个细胞），
因此 <b>UC vs 正常 的对比在这两个数据集里没有检验力</b>：ds3 的效应量虽呈同向（结肠细胞 d=0.45，
深度校正后 d=0.82），p 值均在 0.13–0.97 之间。这两个数据集的价值在于<b>细胞类型定位</b>，不在于疾病对比。</p>
</div>""")

# 5.9 ds4
if ds4 and t4:
    n_ct = len({r["cell_type"] for r in t4})
    n_sig = len([r for r in t4 if r["score"] == "SGRG" and float(r["p_perm"]) < 0.05])
    A(f"""<div class='card'><h3>5.9 大样本验证（ds4，320,979 细胞，8 UC vs 7 正常 donor）</h3>
<p>Tearle 2026（<i>Nat Commun</i>，GEO/CELLxGENE）是本方案里<b>唯一有足够 UC donor 数</b>的队列：
25 种细胞类型达到「UC/正常各 ≥3 donor 且每 donor ≥20 细胞」的纳入标准。
它在<b>原始口径</b>下给出的，恰恰是一个看似漂亮的双向阳性：</p>
<ul>
<li><b>{n_sig} / {n_ct} 种细胞类型 p&lt;0.05</b>，方向几乎全部向上（UC 升高）；</li>
<li><b>杯状细胞 q=0.016（p=0.0006，d=2.17）</b>、Naive B q=0.081（d=2.37）、动脉内皮 q=0.081（d=1.97）、
ILC3 q=0.081（d=1.94）、Th17 q=0.055（d=2.09）；</li>
<li>按 26 种细胞类型做 BH 校正后仍有 1 种（杯状细胞）过线。</li>
</ul>
<p>如果只看这张表，结论会是「<b>应激颗粒程序在 UC 黏膜多种细胞类型中显著升高</b>」——
一篇可以投出去的阳性论文。但它的 donor 级测序深度是：</p>""")
    A(tbl(["组别", "donor 数", "深度中位数（UMI，pseudobulk 加权）", "Mann–Whitney p", "用深度区分两组的 AUC"],
          [["UC", "8", "11,231", "0.029", "0.839"],
           ["正常", "7", "5,189", "", ""]]))
    A(f"""<p><b>深度与疾病在 donor 级高度共线</b>；在细胞类型内部，log10(nUMI) 对疾病标签的 R² 最高 0.65
（杯状细胞 0.65、动脉内皮 0.55、Th17 0.60），中位 0.29。此时残差法不再可靠——
把深度回归掉，等于把与深度共线的疾病效应一起回归掉。</p>
<p>所以改用<b>不依赖模型假设的做法</b>：把每个细胞抽稀到相同 1,000 UMI
（丢弃不足者，保留 241,220 / 320,979 = 75% 细胞），重新打分、重新检验。</p>""")
    A(img("fig_D7_ds4_depth_control", "图 D7｜ds4 三种口径的效应量对照（25 种两种口径共有的细胞类型，按收缩幅度排序）。"
          "左：原始评分——14/26 种 p&lt;0.05，几乎全部向上；中：残差法——收缩到 ≈0；"
          "右：深度匹配抽稀——同样收缩到 ≈0，且深度在两组完全相等（中位 1000 vs 1000）。"))
    rows = []
    for c in ORDER_D4:
        for key, lab in [("SGRG", "原始"), ("SGRG_adj", "残差"), ("@rar", "抽稀")]:
            if key == "@rar":
                r = next((x for x in t4r if x["cell_type"] == c and x["score"] == "SGRG"), None)
                if not r:
                    continue
                d, p, q = r["diff"], r["p_perm"], r["q_BH"]
            else:
                r = next((x for x in t4 if x["cell_type"] == c and x["score"] == key), None)
                if not r:
                    continue
                d, p, q = r["diff"], r["p_perm"], r["q_BH"]
            rows.append([SHORT_D4.get(c, c), lab, nf(d, 4), psci(p), psci(q)])
    A(tbl(["细胞类型（收缩幅度前 8）", "口径", "Δ (UC−正常)", "p", "q(BH)"], rows[:24]))
    A(f"""<div class="ok"><b>判定：验证队列同样触发失败条件第 2 条，且证据比发现队列强得多。</b>
深度匹配后 25 种细胞类型无一通过 FDR（最小 q=0.434），效应量中位数从 0.048 → 0.008（收缩 6 倍），
方向也不再一致（12/25 转为下降）。原始口径下的 14 个「显著」是<b>测序深度的伪影</b>，
不是应激颗粒程序的激活。</div>
<p class="small">两条独立证据链指向同一结论：① 残差法与抽稀法给出几乎相同的收缩幅度；
② 收缩幅度最大的细胞类型，正是深度-疾病共线性最强的那些（动脉内皮 R²=0.55、杯状细胞 R²=0.65）。
若「UC 升高」是真实的生物学，它不应在两种不同的去混杂方法下同时消失，也不应与深度的共线性程度排序一致。</p>
</div>""")
    A("<div class='card'><h3>5.10 跨数据集汇总：四种口径，同一个零</h3>")
    A(img("fig_D6_cross_dataset_forest", "图 D6｜四个数据集、32 个「细胞类型 × 数据集」组合的效应量。"
          "左：原始评分；右：深度控制后（ds1/ds2/ds3 用残差法，ds4 用抽稀法）。"
          "唯一在深度控制后仍 p&lt;0.05 的是 ds4 的浆细胞（UC 中<b>更低</b>，q=0.33）。"))
    A("""<p>把四套数据放在一起看，模式是稳定的：<b>原始口径下散布着若干 nominally significant 的点，
深度控制后全部向零收缩</b>，而且没有一个点在 BH 校正后存活。四套数据、三种去混杂方法、
两个独立的研究团队（Garrido-Trigo / Oliver / Tearle）、合计约 46 万细胞——
没有一处出现可重复的面板级升高。</p></div>""")

# 5.11 随机基因集负对照
if rnd_raw and rnd_ds:
    rr, rd_ = rnd_raw, rnd_ds
    sgrg_raw = next(x for x in rr["real"] if x["geneset"] == "SGRG")
    sgrg_ds = next(x for x in rd_["real"] if x["geneset"] == "SGRG")
    exp1, exp2 = rr["expected_false_positives"], rd_["expected_false_positives"]
    A(f"""<div class='card'><h3>5.11 随机基因集负对照：把「阳性」钉死为伪影</h3>
<p>5.9 的解释留下一个必须堵死的后门：<b>「抽稀后全阴性，会不会是抽稀本身把信号洗掉了？」</b>
抽稀必然损失信息，若只有这一条证据，审稿人有理由质疑。所以做了正式的负对照——
生成 <b>200 个与 SGRG 表达分层匹配的随机基因集</b>（同等大小 K=103、同表达分箱构成、
排除 SGRG 基因本身），让它们<b>通过与 SGRG 完全相同的处理链</b>（同一套细胞、同一归一化、
同一对照基线、同一检验）。因为评分对基因集是线性的，200 个集与 SGRG 在一次矩阵乘法中并行算完。</p>
<h4>结果一：未校正口径下，任何基因集都「显著」</h4>""")
    A(tbl(["口径（26 种细胞类型，期望假阳性 1.3）", "随机集显著数（中位）", "膨胀倍数", "SGRG", "SGRG 在随机分布中的经验 p"],
          [["Mann–Whitney U", "6（5–95%: 3–10）", "4.6×", "5", "0.856"],
           ["均差置换检验", "11（5–95%: 5–16）", "8.5×", "13", "0.488"],
           ["效应量 |Δ| 中位", "0.037", "—", "0.044", "0.19（81 百分位）"]]))
    A(f"""<p><b>SGRG 的全部指标都落在随机基因集的零分布之内，且低于中位数。</b>
最刺眼的是方向性：随机基因集在 UC 中「升高」的细胞类型中位数是 <b>23/26</b>——
未校正口径下任何基因集都会「几乎全部细胞类型升高」，这正是深度混杂的指纹，与基因集内容无关。</p>
<h4>结果二：抽稀口径下，假阳性和方向偏差同步消失</h4>""")
    A(tbl(["口径（25 种细胞类型，期望假阳性 1.2）", "随机集显著数（中位）", "膨胀倍数", "SGRG", "经验 p"],
          [["Mann–Whitney U", "2", "1.6×", "3", "0.388"],
           ["均差置换检验", "3", "2.4×", "2", "0.751"],
           ["方向性 Δ&gt;0 中位", "13/25（≈一半，无偏）", "—", "12/25", "—"]]))
    A(img("fig_D8_random_control", "图 D8｜随机基因集负对照（未校正口径）。"
          "(A)(B) 200 个随机基因集的「显著细胞类型数」分布：Mann–Whitney 口径中位 6（期望 1.3，膨胀 4.6×）、"
          "置换口径中位 11（膨胀 8.5×），SGRG（红线）都落在分布中部。"
          "(C) SGRG 的效应量位于随机分布第 81 百分位。(D) 逐细胞类型：SGRG（红钻）几乎全部落在 "
          "200 个随机集的 5–95% 区间内。"))
    A(img("fig_D9_randcalib", "图 D9｜深度匹配前后的校准对比（以随机集为真零）。"
          "(A)(B) 假阳性膨胀从未校正的 4.6×/8.5× 回到抽稀后的 1.6×/2.4×，接近名义水平；"
          "(C) 细胞类型内 UC/对照方差比从 0.53 恢复到 0.72（1.0=两组可交换）；"
          "(D) 方向偏差从 23/26「全向上」恢复到 13/25（≈一半）。"
          "抽稀同时移除了假阳性与方向偏差——所以「抽稀后全阴性」是真实的，不是检验力损失。"))
    A(f"""<h4>结果三：一个意外收获——置换检验本身在这个数据上不校准</h4>
<p>同一份 pseudobulk 数据，置换检验给出 14/26 显著而 Mann–Whitney 给出 5/26；24/26 个细胞类型中
置换检验的 p 系统性更小。机制：UC 的 pseudobulk 方差只有对照的一半（sd 比中位 0.49）——
UC 细胞更多、更深，平均后噪声更小，两组<b>不可交换</b>，而均差置换检验假设可交换。
以随机集为真零直接测出：置换检验膨胀 8.5 倍、Mann–Whitney 膨胀 4.6 倍。
也就是说 <b>5.9 的「14/26」里没有一行是 SGRG 特有的</b>——随机基因集也能拿到 11/26。</p>
<div class="ok"><b>判定：负对照把混杂解释从「合理」升级为「完备」。</b>
未校正口径下的全部「阳性」（显著数、方向一致性、最小 p、效应量）都能被 200 个随机基因集复现；
抽稀口径下这些伪影与随机集同步归零——抽稀没有制造假阴性，它移除的是假阳性。
「SGRG 在 UC 升高」至此在单细胞分辨率上被排除到与随机噪声无法区分的程度。</div>
<p class="small">注：阳性对照（Hallmark 炎症/IFN-γ 通路）在 donor 级也未超出随机分布
（IFN-γ 7/25 vs 随机中位 2，是唯一高于随机的真实通路，与其在 UC 中的先验一致但同样不达 FDR）——
这与「donor 级对比被技术变异主导」一致，不单独解读为生物学。</p>
</div>""")

else:
    A("""<div class="card"><h3>5.9 大样本验证（ds4，320,979 细胞）——进行中</h3>
<p>Tearle 2026（<i>Nat Commun</i>）的 320,979 个结肠黏膜细胞正在下载（2.7 GB），
该队列是 <b>正常 / UC / PSC-UC 三种状态、四个结肠部位</b>，是本方案里唯一有足够 UC donor 数的数据源。
分析脚本与判据已固定，数据到位即可直接跑出同一组表格。</p>
</div>""")

# ================= 5.12 来源审计 =================
BA = {r["dataset"]: r for r in baudit}
if BA:
    A("""<div class='card'><h3>5.12 来源审计：把混杂推到比深度更上游的一层</h3>
<p>5.9–5.11 把 ds4 的「阳性」解释为深度伪影，并做了完备的负对照。<b>但这条解释有一个从未检查过的前提：
「正常」对照组本身是不是一个同质的组。</b>本轮对四套数据统一做了来源审计
（<code>partD_batch_audit.py</code>），结果如下：</p>""")
    rows = []
    for ds in ["ds1_garrido", "ds2_myeloid", "ds3_largeintestine", "ds4_tearle"]:
        r = BA.get(ds)
        if not r:
            continue
        same = (r.get("frac_ctrl_cells_samebatch_assay") if ds == "ds1_garrido" else
                r.get("frac_ctrl_cells_samebatch_study") if ds != "ds4_tearle" else
                r.get("frac_ctrl_cells_samebatch_batch_condition"))
        rows.append([
            {"ds1_garrido": "ds1 Garrido-Trigo", "ds2_myeloid": "ds2 Oliver 髓系",
             "ds3_largeintestine": "ds3 Oliver 大肠", "ds4_tearle": "ds4 Tearle PSC/UC"}[ds],
            r["n_uc_donors"], r["n_ctrl_donors"],
            "%.2f×" % float(r["depth_fold_uc_over_ctrl"]),
            ("%d%%" % round(100 * float(r["frac_ctrl_cells_primary"]))) if r.get("frac_ctrl_cells_primary") not in (None, "", "nan") else "—",
            "%d%%" % round(100 * float(same)) if same not in (None, "", "nan") else "—",
            r.get("n_batch_study") or r.get("n_batch_batch_condition") or r.get("n_batch_assay"),
            r.get("n_batch_mixed_study") or r.get("n_batch_mixed_batch_condition") or r.get("n_batch_mixed_assay")])
    A(tbl(["数据集", "UC donor", "对照 donor", "深度比 (UC/对照)",
           "对照中本研究原产细胞", "对照中「与病例同批次」细胞", "批次数", "同时含病例与对照的批次"], rows))
    A(f"""<div class="warn"><b>ds4：对照组的 94.4% 不是这项研究的对照。</b>
CELLxGENE 的 <code>is_primary_data</code> 字段把 320,979 个细胞一分为二：
<b>病例侧（UC 62,968 / PSC 69,510）100% 是本研究原产</b>；
而「正常」侧只有 <b>10,519 个细胞（HC-1 / HC-2，5.6%）</b>是原产——
其余 <b>177,982 个细胞来自 5 个复用 donor（411c/414c/417c/432c/440c）</b>，
批次名 A3x-*/CRC-*，为切除标本、10x 5' v2 化学版本、细胞中位深度 289–9,769 UMI。
原文集合描述写的是「healthy control individuals (n=2)」——与 5.6% 这个数字完全对得上。</div>
<h4>重算：上一轮的「深度混杂」其实是「对照组成」</h4>""")
    A(tbl(["对照定义", "对照 donor 数", "细胞级中位深度", "UC/对照 深度比", "说明"],
          [["全部 normal（上一轮口径）", "7", "1,305", "4.35×", "94.4% 为复用低深数据"],
           ["仅复用部分（is_primary=False）", "5", "1,252", "4.53×", "全部来自 A3x/CRC 切除批次"],
           ["仅本研究对照（is_primary=True）", "2", "5,440", "<b>1.04×</b>", "HC-1（2025 批次）+ HC-2（CRC 批次）"]]))
    A("""<p><b>深度差在把对照组换成研究自己的对照后彻底消失（1.04×）。</b>
也就是说 5.9 里「UC 深度是正常 1.65 倍、donor 级 AUC=0.84」这个混杂，其上游原因是
<b>对照被填入了不同来源、不同流程、低一个数量级深度的复用数据</b>；
抽稀到 1,000 UMI 之所以能「消除」全部信号，是因为它把两组强行拉到同一个低深度上——
这是在拿病例数据与别人的数据比，而不是在拿病例与对照比。</p>
<h4>ds2 / ds3：同样的问题，只是更隐蔽</h4>
<p>Oliver 2024 是 11–25 个研究的聚合图谱。UC 只出现在其中
<b>2–3 个研究</b>（ds2：Kinchen2018 / Parikh2019 / Uzzan2022，共 9 donor；ds3：同三个研究，共 8 donor），
而对照来自全部 11–25 个研究。<b>对照细胞中只有 1.5%（ds2）/ 19.3%（ds3）来自「也含病例」的研究</b>——
我们此前的 ds2/ds3 检验实质上是「UC（3 个研究）vs 正常（20+ 个研究）」的跨研究对比。
ds1 是唯一没这个问题的（单一研究、病例对照都在 10x 3' v3 内、90% 对照细胞与病例同化学版本），
但它只有 6 donor/组。</p>""")
    A(img("fig_D10_provenance_confounding", "图 D10｜对照组的来源组成如何决定结论。"
          "(A) ds4 细胞级深度分布：UC / PSC / 本研究对照 / 复用「正常」——复用组低一个数量级。"
          "(B) donor 级中位深度：换成研究自己的对照后，UC/对照深度比从 4.53× 塌缩到 1.04×。"
          "(C) 四套数据的对照同源性：ds1 90%、ds3 19%、ds4 4%、ds2 1.5%。"
          "(D) 同一份数据、两种对照定义的 SGRG 效应量：灰（朴素口径）21/26 种细胞类型 UC 更高；"
          "红（同源口径）只有 8/23——<b>方向都翻了</b>。"))

# ================= 5.13 检验力裁决 =================
if prov and spv:
    A("""<div class='card'><h3>5.13 检验力裁决：三种口径，三种死法</h3>
<p>既然「UC vs 全部 normal」不可用，那还有什么口径能回答这个问题？把所有可行的设计在同一遍扫描里都跑了一遍
（同一套细胞、同一套 200 个表达分层匹配的随机基因集做零分布校准），结果：</p>""")
    SNAME = {"A_naive_all_normal": "A_naive", "B_primary_only": "B_primary", "C_inflamed": "C_inflamed"}

    def _row(cname, lab, ndesc):
        r = [x for x in spv if x["contrast"] == SNAME.get(cname, cname)]
        if not r:
            return None
        real = [x for x in r if x["kind"] == "real" and x["geneset"] == "SGRG"][0]
        rand = [x for x in r if x["kind"] == "random"]
        pos = [x for x in r if x["kind"] == "positive"]
        n_ct = int(float(real["n_ct"]))
        import numpy as _np
        rmed = float(_np.median([float(x["n_p05"]) for x in rand])) if rand else float("nan")
        pmed = float(_np.median([float(x["n_p05"]) for x in pos])) if pos else float("nan")
        # 方向性（C 口径用 mean_diff）
        T = tpv if cname != "C_inflamed" else tin
        TNAME = {"A_naive_all_normal": "A_naive_all_normal",
                 "B_primary_only": "B_primary_only",
                 "C_inflamed": "C_within_donor_inflamed"}
        kk = "SGRG"
        d = [x for x in T if x["contrast"] == TNAME.get(cname, cname) and x["geneset"] == kk]
        up = int(sum(1 for x in d if (float(x["diff"]) if "diff" in x else float(x.get("mean_diff", 0))) > 0))
        pp = real.get("n_p05_perm")
        return [lab, str(n_ct), "%d (置换 %s)" % (int(float(real["n_p05"])),
                ("—" if pp in (None, "", "nan") else str(int(float(pp))))),
                "%.0f（期望 %.1f）" % (rmed, n_ct * 0.05),
                "%d" % int(pmed) if pmed == pmed else "—",
                "%d/%d" % (up, len(d))]
    rr = [x for x in (_row("A_naive_all_normal", "A 朴素：UC vs 全部 normal", ""),
                      _row("B_primary_only", "B 同源：UC vs 本研究自己的对照（n=2）", ""),
                      _row("C_inflamed", "C 受控：同 donor 炎症 vs 非炎症（4 对）", "")) if x]
    A(tbl(["口径", "细胞类型", "SGRG 显著数", "随机集显著数（中位）",
           "阳性对照显著数", "方向（病例更高）"], rr))
    A(f"""<p><b>三个口径各有致命缺陷，而且都能被阳性对照诊断出来：</b></p>
<ul>
<li><b>A（朴素）</b>：SGRG 21/26 种细胞类型 UC 更高、5/26 显著——但随机基因集中位也有 6/26（膨胀 4.6×），
SGRG 经验 p=0.86，落在随机分布之内。这是 5.9–5.11 已经钉死的伪影。</li>
<li><b>B（同源对照，n=2）</b>：随机集的显著数中位 1（膨胀 <b>0.9×</b>，校准良好），
SGRG 1/23、经验 p=0.85。<b>但五个炎症通路阳性对照全部 0/23 显著</b>——
连 NF-κB / IFN-γ / 缺氧 / UPR 这些在 UC 必然升高的程序都测不出来，
说明这个口径<b>没有检验力</b>（2 例对照撑不起 23 种细胞类型的组间比较）。
方向也从 A 的 21/26 反转为 8/23。<b>「0 显著」在这里读作「测不出」，不是「不存在」。</b></li>
<li><b>C（同 donor 炎症配对，4 对）</b>：这是唯一在结构上受控的设计
（同一个 donor、同一个批次、同一次解离，且 inflamed/noninflamed 深度比中位仅 1.17）。
方向上与先验一致：五个炎症通路在 13–15/15 种细胞类型中升高、SGRG 在 11/15 中升高。
但 <b>4 例配对设计的最小可达双侧 p = 2/16 = 0.125</b>——
所有检验（包括阳性对照）都被卡在这个地板上，0/15 显著。<b>同样是零检验力。</b></li>
</ul>""")
    A(img("fig_D11_power_diagnosis", "图 D11｜检验力裁决。"
          "(A) 三种口径下 SGRG 的显著数（红）对照 200 个随机基因集的零分布（灰）：A 口径落在分布中部；"
          "B/C 口径随机集显著数≈0（校准良好），但那是因为整个检验没有检验力。"
          "(B) 阳性对照诊断：A 口径能测出 2–7 个炎症通路；<b>B 与 C 口径连一个炎症通路都测不出来</b>。"
          "(C) 配对设计的检验力天花板：4 对 donor 的最小可达 p=0.125。"
          "(D) 唯一受控设计的逐细胞类型配对差：方向与先验一致，但全部被 p 地板卡住。"))
    A("""<div class="ok"><b>判定（修订 5.9–5.11 的结论）：</b>
Part D 的四套数据<b>在结构上无法回答「SGRG 面板在 UC 黏膜是否被激活」</b>：
<ul>
<li>ds1：设计干净，但 6 donor/组，检验力不足以排除中等效应（这是 5.5 就写明的）；</li>
<li>ds2 / ds3：病例与对照来自不同研究，疾病对比跨研究混杂；</li>
<li>ds4：对照组 94.4% 是复用的外部数据，疾病与批次完全嵌套（35 个批次只有 1 个混合），
把对照换成研究自己的 2 例后深度差消失，但检验力归零；唯一受控的同 donor 炎症配对只有 4 例。</li>
</ul>
因此，此前所有「阴性 p 值」的正确读法是<b>「在这些设计下无法检出」</b>，而不是「面板未被激活」。
反过来，5.9 的「14/26 显著」也不能被复活——它在随机集零分布之内，且方向随对照定义翻转。</div>
<p class="small">方法学要点：单细胞图谱里 <code>is_primary_data</code> 标记的是「该细胞是否属于本集合的主数据」，
不是「是否属于这项研究的对照」。当一个 atlas 把外部对照并进来时，对照会被无声地替换成低深度、
不同流程的数据，而下游的所有深度校正、批次校正都建立在一个错误的分组之上。
这个检查成本是零（读一列 obs），但四套主流 UC 图谱里三套都没过。</p>
</div>""")
# ================= 5.14 技术轴清单（步 1：A2 宽表的等构对应物） =================
inv = tload("partD_technical_inventory.tsv")
if inv:
    IM = {r["dataset"]: r for r in inv}

    def _g(d, k, fmt="%.3f", dash="—"):
        try:
            v = float(IM[d][k])
        except Exception:
            return dash
        return dash if v != v else (fmt % v)

    A("""<div class='card'><h3>5.14 技术轴清单：把 A2 的队列级宽表搬到单细胞分辨率</h3>
<p>A2 在 bulk 队列上做的事——给每个队列量出「技术轴占据度」与「技术轴-分组共线度」——
现在对四套单细胞数据各做了一份等构的清单（<code>partD_technical_inventory.py</code>，
不重跑 h5ad，全部从既有结果表汇总）。这一步把 A2 与 Part D 从「比喻关系」焊成「方法学对称」：
两张宽表可以用同一套列名并排陈述。</p>""")
    A(tbl(["数据集", "细胞", "UC/对照 donor", "深度比 donor级 / 细胞级", "深度 AUC",
           "ρ(评分,深度) 中位/最大", "R²(深度~疾病) 中位/最大",
           "主轴技术占据", "主轴疾病占据", "批次嵌套"],
          [[d,
            "%s" % f"{int(float(IM[d]['n_cell_pseudobulk'])):,}",
            "%s/%s" % (IM[d]["n_uc_donor"], IM[d]["n_ctrl_donor"]),
            "%s / %s" % (_g(d, "depth_fold_uc_over_ctrl", "%.2f×"), _g(d, "depth_fold_celllevel", "%.2f×")),
            _g(d, "depth_separability"),
            "%s / %s" % (_g(d, "rho_score_nUMI_med"), _g(d, "rho_score_nUMI_max")),
            "%s / %s" % (_g(d, "r2_depth_disease_med"), _g(d, "r2_depth_disease_max")),
            _g(d, "pc_tech_share"), _g(d, "pc_disease_share"),
            "是（%s）" % IM[d]["batch_worst_var"] if IM[d]["batch_nested_flag"] == "1" else "否"]
          for d in ["ds1", "ds2", "ds3", "ds4"]]))
    A("""<p><b>清单读出三条结构性事实：</b></p>
<ul>
<li><b>「技术轴占据主轴」与「技术轴与分组共线」是两个独立的指标</b>：
ds1 的主轴技术占据 0.40（PC1 与深度的 R²=0.64），但深度对疾病的 AUC 只有 0.56（无共线）——
技术轴占了主轴却没混杂分组；ds4 两者都高（占据 0.48、共线 AUC 0.91）。
A2 里同样有这两类队列。任何「占不占 PC1」的单一指标都不能代替对共线性的直接测量。</li>
<li><b>评分-深度耦合是普适的，共线不是</b>：四套数据的细胞类型内 ρ(面板评分, nUMI) 中位
0.58–0.81、最高 0.98——评分随深度上升这件事在所有数据都成立；
但「深度在 donor 级与疾病共线」（R² 中位 0.30）只有 ds4。所以 ds4 是唯一一个
「耦合 + 共线」同时成立的验证队列，也是唯一产出大规模名义阳性的队列。</li>
<li><b>「UC 测得更深」不是稳健事实</b>：donor 级中位与细胞级中位两个口径下，
ds1 为 0.86× vs 1.52×、ds3 为 1.28× vs 0.65×（方向都反了）；ds4 为 1.77× vs 4.35×，
而换成研究自己的对照后是 1.04×（5.12）。深度比是一个<b>口径敏感的描述量</b>，
引用时必须同时给出口径，不能作为论点的承重墙。</li>
</ul>""")
    A(img("fig_AD4_technical_inventory",
          "图 AD4｜单细胞数据集级技术轴清单。(A) 深度比在 donor 级与细胞级两个口径下的对比"
          "（绿 = ds4 仅用本研究自有对照），标注为深度区分疾病的 AUC——深度比不稳健，"
          "且「UC 更深」在 ds1/ds2/ds3 都不成立或反向。(B) 主轴技术占据 × 深度-疾病可分性："
          "ds4 是唯一落在右上（高占据 + 高共线）的数据集，灰线为深度口径的敏感性范围。"
          "(C) ds4 四道干预下的名义显著细胞类型数：名义 14/26 → 抽稀 1/25 → 仅自有对照 1/23 → "
          "同 donor 炎症配对 0/15（虚线 = 5% 期望假阳性）。"))
    A(tbl(["干预", "技术轴变化", "名义显著（p&lt;0.05）", "FDR&lt;0.05", "随机集校准"],
          [["抽稀到 1,000 UMI", "两组深度强行相等（1000 vs 1000）",
            "14/26 → 1/25", "1 → 0（最小 q=0.397）", "随机集同步归零（中位 11 → 3）"],
           ["对照换成研究自有（n=2）", "深度比 4.35× → 1.04×",
            "14/26 → 1/23", "0", "随机集中位 1（膨胀 0.9×，校准良好）"],
           ["同 donor 炎症配对（4 对）", "深度比 1.17×，同一批次同一解离",
            "14/26 → 0/15", "0", "p 地板 0.125（4 对配对的最小可达双侧 p）"]]))
    A("""<div class="ok"><b>这一节的价值：把「干预-响应」从 ds4 的个案升级为一张可复用的清单。</b>
A2 的宽表说「bulk 队列的批次占据度可以逐队列量出来」；本节说「单细胞数据集也可以，且量出来的结果
恰好解释了哪个数据集会产出名义阳性（ds4）、哪个不会（ds1/ds2/ds3）」。
两份宽表 + 两个干预-响应证据，构成方法学主论点的全部数据支撑。</div>
<p class="small">口径注意：主轴占据度（pc_tech_share）在 ds1 因 donor×细胞类型仅 12×5 而不稳定
（深度口径换成细胞数加权后从 0.40 掉到 0.09），ds2/ds3/ds4 两口径差 &lt;0.03。
明细与全部列定义见 <code>results/partD_technical_inventory.tsv</code> 与
<code>partD_technical_inventory_meta.json</code>。</p>
</div>""")

# ================= 5.15 公共资源可用性全库审计 =================
cxg_summ = jload("partD_resfeas_cxg_summary.json", {})
geo_summ = jload("partD_resfeas_geo_summary.json", {})
verdict = jload("partD_resfeas_verdict.json", {})
ledger = tload("partD_resfeas_ledger.tsv")
if cxg_summ and geo_summ:
    A("""<div class='card'><h3>5.15 公共资源可用性全库审计：「你真的找遍了吗」</h3>
<p>5.12-5.14 的结论是「公共单细胞资源在结构上无法做出有检验力、批次受控的检验」。
审稿人的第一句反问必然是：<b>「你真的找遍了吗？也许有一套有检验力的数据你没看。」</b>
本节把这句话变成可审计的枚举——不是检索关键词数一篇综述，而是把两个公共库逐条拆开，
给每个候选资源过一遍四条闸门。</p>
<p><b>四条闸门（缺一不可）：</b>① assay 测 mRNA（否则 SGRG 这类编码基因面板根本打不了分）；
② 含 UC；③ 有受控对比且达到检验力（同供体炎症/非炎症配对 ≥6 对，或同研究真对照每组 ≥5 例）；
④ 在体组织（排除类器官/细胞系）。检验力门槛来自配对设计的数学事实：
双侧符号翻转检验最小可达 p = 2/2<sup>n</sup>（n=4 → 0.125，n=6 → 0.031）。</p>""")
    A(tbl(["库", "枚举范围", "含 UC 的资源", "有配对设计", "≥6 对", "通过全部闸门"],
          [["CELLxGENE", "%d 集合 / %d 肠道数据集" % (cxg_summ["n_collections_total"],
                                                     cxg_summ["n_datasets_gut"]),
            "%d 集合（%d 数据集）" % (cxg_summ["n_collections_with_UC"],
                                     cxg_summ["n_datasets_gut_has_UC"]),
            "—（CELLxGENE 无配对元数据，按同研究真对照判定）", "—",
            "<b>%d</b>" % verdict.get("cxg_usable_controlled_resource", 0)],
           ["GEO", "%d 人类系列（IBD×肠段）" % geo_summ["n_gse_human"],
            "%d 个 mRNA 可打分系列中含 UC" % geo_summ["n_series_pairing_and_UC_mrna_ge6"],
            "%d 个（mRNA 可打分）" % geo_summ["n_series_with_pairing_mrna"],
            "<b>%d</b>" % geo_summ["n_series_with_pairing_ge6_mrna"],
            "<b>1</b>"]]))
    A("""<p><b>枚举过程中的三个方法学陷阱（每个都改变了一次结论）：</b></p>
<ul>
<li><b>陷阱 1：系列摘要 ≠ 实际检测分子。</b>GSE272890 有 43 对同供体炎症/非炎症配对——
全库最大的配对设计，一度让「找不到达标数据」的判断动摇。但它这个 series 实际是
<i>Non-coding RNA profiling</i>（测 miRNA），摘要里写的「spatial transcriptomics」是论文里的
另一部分。按摘要关键词归类会把它虚记成可用空间数据。<b>可审计字段是 <code>!Series_type</code>
，不是摘要文本</b>。加「是否测 mRNA」闸门后，≥15 对的设计全库归零。</li>
<li><b>陷阱 2：配对 ≠ 同部位配对。</b>GSE312415/312420（空间转录组，CosMx/Xenium）里
同一供体确实同时有炎症与非炎症样本，但配对是<b>跨肠段</b>的——回肠非炎症 vs 直肠炎症，
状态与部位完全混杂，无法归因。按「同供体 + 同部位」的严格口径，空间数据同部位配对 = 0。</li>
<li><b>陷阱 3：供体标识在样本标题里，且不能激进归一化。</b>GEO 的配对信息只存在于样本标题
（characteristics 常常只有 tissue）。而把 <code>353_6_3338</code> 削成 <code>353</code>
这类归一化会凭空造出配对。本审计按原样 token 比对，再把全库仅有的配对系列人工核查一遍。</li>
</ul>""")
    A(img("fig_D12_resource_feasibility",
          "图 D12｜公共资源可用性全库审计。(A) GEO 枚举漏斗：521 个人类 IBD×肠段系列 → 376 个测 mRNA → "
          "130 个单细胞/空间 → 13 个有同供体配对 → 只有 4 个 ≥6 对。(B) 配对设计按分辨率分布："
          "达标设计全部落在 bulk（实心 = 测 mRNA）；最大配对数 GSE272890 的 43 对测的是 miRNA（空心，不可用）；"
          "单细胞最多 5 对、空间同部位 0 对，全部低于 6 对检验力门槛。(C) CELLxGENE 全库只有 3 个集合含 UC："
          "只有 ds1（6 UC + 6 对照、100% 原产）落在「对照 ≥5 例且原产率 ≥0.9」的达标区。"
          "(D) 四条闸门判据矩阵：候选资源里通过全部闸门的只有 2 个。"))
    if ledger:
        rows = []
        for r in ledger:
            mark = lambda b: "✓" if str(b) in ("True", "1", "true") else "✗"
            npair = r.get("n_pairs") or ""
            ncase = r.get("n_case") or ""
            nctrl = r.get("n_ctrl") or ""
            rows.append([r["source"], str(r["resource"])[:44], r["resolution"], r["disease"],
                         npair if npair else ("%s/%s" % (ncase, nctrl) if ncase else "—"),
                         mark(r["G1_mRNA"]), mark(r["G2_UC"]), mark(r["G3_powered_control"]),
                         mark(r["G4_in_vivo"]),
                         "<b style='color:#2e6f9e'>USABLE</b>" if r["USABLE"] in ("True", "1", "true")
                         else "fails", str(r["reason"])[:90]])
        A(tbl(["库", "资源", "分辨率", "疾病", "配对/例数", "①mRNA", "②UC", "③检验力", "④在体",
               "结论", "未过原因"], rows))
    A("""<p><b>裁决：全库只有 2 个资源能支撑有检验力且受控的 SGRG 检验。</b></p>
<ul>
<li><b>ds1（CELLxGENE，单细胞，6 UC vs 6 对照，100% 原产）</b>——已在 Part D 主线检验：阴性，
且经随机集校准（5.11）确认非抽稀假阴性。</li>
<li><b>GSE133060（GEO，Agilent GPL4133 全基因组芯片，8 对同供体炎症/非炎症 UC 黏膜）</b>——
<b>此前从未被任何一方检验过</b>。8 对的双侧符号秩最小 p = 2/256 = 0.0078，有名义检验力。
这是本审计顺手发现的一个「存量资产」：一个现成的、有检验力的、批次天然受控的 UC 配对设计，
静静躺在 bulk 层面。</li>
</ul>
<p>由此原判定需要一句精确的修订：<b>「公共单细胞资源在结构上无法做出有检验力、批次受控的检验」成立；
但公共资源整体（含 bulk）存在一个未被检验的达标设计</b>——
而检验它——结果见 <b>5.16</b>。</p>
<p class="small">脚本：<code>partD_resfeas_cxg.py</code>（CELLxGENE 全库 397 集合 → 40 个候选 →
159 个肠道数据集）、<code>partD_resfeas_geo.py</code>（GEO 四条互补查询 → 521 人类系列 →
423 个解析 SOFT，被跳过的 98 个全为 bulk 且无单细胞/空间，覆盖完整）、
<code>partD_resfeas_ledger.py</code>（裁决表 + 图）。明细：
<code>partD_resfeas_cxg_datasets.tsv</code> / <code>partD_resfeas_geo_series.tsv</code> /
<code>partD_resfeas_geo_samples.tsv</code> / <code>partD_resfeas_ledger.tsv</code>。</p>
</div>""")


# ================= 5.16 GSE133060 实际检验 =================
v133 = jload("partD_gse133060_verdict.json", {})
if v133:
    A("""<div class='card'><h3>5.16 把全库唯一达标的 UC 配对设计测掉：GSE133060</h3>
<p>5.15 的全库枚举留下一个悬而未决的资产：<b>GSE133060</b>——Agilent GPL4133 全基因组芯片，
<b>8 例 UC 患者，每人「炎症黏膜 vs 非炎症黏膜」各一针</b>，同一组织类型、同一批次，
是全库唯一「UC + ≥6 对同供体配对 + 测 mRNA」的公共资源，且从未被检验。
本节把它测掉——这是把 Part D 的故事从「我们找不到数据」推进到「我们把唯一能测的测了」的最后一步。
8 对的双侧符号秩最小 p = 2/256 = 0.0078，名义上<b>有检验力</b>。</p>
<p><b>第一步：设计有效吗？——这是全项目第一个有效力的对比。</b>
在 Part D 此前检验过的所有单细胞口径里，连 NF-κB / IFN-γ 这类已知必然变化的通路都测不出来
（零检验力，见 5.13）。而在这里，<b>5 条 Hallmark 炎症通路全部 8/8 患者同向上调
（Δz = 0.27–0.82，p = 0.0078，即检验力下限）</b>。炎症 vs 非炎症是真实的生物学对比，
该设计能检测它。这是本项目第一次具备「阳性对照通得过」的前提。</p>""")
    A(tbl(["集合","K（可测基因）","配对 Δz（炎症−非炎症）","上调患者数","名义 Wilcoxon p"],
          [["SGRG 面板","%d / %d" % (v133["panel_genes_mapped"], v133["panel_genes_in_panel"]),
            "%+.3f" % v133["panel_mean_paired_z"], "%d / 8" % v133["panel_n_up"],
            "%.4f" % v133["wilcoxon_p_nominal"]],
           ["H_ INFLAMMATORY_RESPONSE（阳性对照）","197","+0.661","8 / 8","0.0078"],
           ["H_ TNFA_SIGNALING_VIA_NFKB（阳性对照）","194","+0.781","8 / 8","0.0078"],
           ["H_ INTERFERON_GAMMA_RESPONSE（阳性对照）","191","+0.815","8 / 8","0.0078"],
           ["H_ HYPOXIA（阳性对照）","186","+0.267","8 / 8","0.0078"],
           ["H_ UNFOLDED_PROTEIN_RESPONSE（阳性对照）","105","+0.655","8 / 8","0.0078"]]))
    A("""<p><b>第二步：名义 p 依然不可信。</b>SGRG 名义 p = 0.0156，看似「显著」。
但 500 个表达分层匹配随机基因集里，<b>250 个（50%）也达到名义 p &lt; 0.05</b>——
因为在 8 对配对设计里，一个强大的全局轴（炎症浸润稀释非免疫转录组，全基因组平均 Δz = −0.101）
使得「≥7/8 患者同向」成为常态而非例外。<b>名义 p 在这里膨胀约 10 倍；
即使设计有效力，校准仍然不是可选项。</b></p>
<p><b>第三步：校准后的裁决。</b>SGRG 的 Δz = +0.208 在 500 个匹配随机集中<b>排名第一（rank 1/501）</b>：
<ul>
<li><b>单侧（方向预注册：炎症中上调）经验 p = 0.002</b>——没有任何随机集达到这个水平；</li>
<li>双侧经验 p = 0.15——零分布不对称（74 个随机集反向达到同等幅度），该口径下不可分辨；</li>
<li>扣除全局轴后（Δz = +0.309），双侧经验 p = 0.006。</li>
</ul>
<b>效应量本身很小：SGRG 的 Δz = 0.21，只有最弱阳性对照（Hypoxia，0.27）的 3/4、
最强（IFN-γ，0.82）的 1/4。</b></p>""")
    A(img("fig_D13_gse133060",
          "图 D13｜GSE133060（8 对 UC 配对）的 SGRG 检验。(A) 500 个表达分层匹配随机集的经验零分布："
          "零分布因全局轴而整体左偏（中位 −0.10），SGRG（红线）位于分布最右端、排名 1/501。"
          "(B) 效应量对比：SGRG 的 Δz = 0.21 通过随机集 q95（虚线），但仅为最弱炎症通路（Hypoxia）的 3/4。"
          "(C) 逐患者配对轨迹：7/8 患者炎症黏膜评分更高。"))
    A("""<p><b>裁决：这是全项目第一个「SGRG 不明显失败」的数据集——但只能算弱阳性。</b></p>
<ul>
<li><b>方向一致且通过校准</b>：7/8 患者上调，在表达匹配的随机集零分布中排名第 1（单侧经验 p = 0.002，
去全局轴后双侧 p = 0.006）。</li>
<li><b>但效应量小</b>（Δz = 0.21，约为炎症通路的 1/4–3/4），且双侧口径（p = 0.15）下无法排除
「随机集也能反向达到同等幅度」。</li>
<li><b>它不推翻 Part D 的主线判定</b>：单细胞层面「结构性无法检验」的结论不变；
它只是说明——当唯一一个有检验力的设计被找到并检验时，信号是<b>弱的、方向性的、依赖校准口径的</b>，
而非教科书式的强阳性。</li>
</ul>
<p><b>对论文的意义</b>：这一节同时解决了两个审稿风险——
①「框架是否永远说 no」：不是，它在有效力的设计上给出了一个通过校准的弱阳性；
②「你找到了为什么不测」：测了，结果如实报告。</p>
<p class="small">脚本：<code>partD_gse133060_test.py</code>（解析 16 个 Agilent FE 文件 → 探针聚合 →
log2 → quantile → z-score → 配对 Wilcoxon + 500 随机集校准）、<code>partD_gse133060_diag.py</code>
（零分布诊断）。明细：<code>partD_gse133060_gene_scores.tsv</code> /
<code>partD_gse133060_tests.tsv</code> / <code>partD_gse133060_verdict.json</code>。
数据：GSE133060_RAW.tar（16 个 FE 文件 + GPL4133 注释）。</p>
</div>""")

# ================= 5.17 全库来源审计（Rung 1 全库尺度） =================
lw = jload("partD_resfeas_libwide_summary.json", {})
if lw:
    A("""<div class='card'><h3>5.17 把来源审计推到全库尺度：可审计性本身与数据质量混杂</h3>
<p>5.12 的发现（ds4 对照臂 94.4% 复用外部供体）靠读<b>单个数据集的一列元数据</b>。它引出一个全库问题：
这种检查<b>有多少机会可行</b>？我们把 <code>partD_resfeas_libwide.py</code> 从 40 个肠道偏置集合
（冒烟产物，已废弃）扩到<b>全库 397 个集合 / 2,238 个数据集 / 292.8M 细胞</b>
（无需下载任何 h5ad），并把每个数据集按疾病标注分为三类：
<code>ctl_only</code>（仅对照标签）/ <code>dis_only</code>（仅疾病标签）/
<b><code>mixed</code>（同一数据集内同时含对照与疾病标签——对比存在于数据集内部，臂无法从元数据拆分）</b>。</p>
<p><b>三个数字刻画全库：</b></p>
<ul>
<li>细胞级复用是常态：<b>41.3%</b> 的细胞（人类 37.5%）非所属数据集原产；
<b>36.3%</b> 的数据集（812/2238）完全由复用细胞构成。</li>
<li><b>★ 携带病例-对照对比的数据集，恰是最无法归属的一类</b>：
346 个 mixed 数据集（占全库 15.5%）的 frac_primary 中位仅 <b>0.06</b>；
人类 mixed 数据集的 <b>49.8%</b>（159/319）完全复用（复用细胞占比 34.5%）。</li>
<li><b>★ 臂级审计几乎不可能</b>：140 个含对比集合中，只有 <b>15 个（10.7%）</b>把两臂存成独立数据集；
其余 <b>125 个（89.3%）</b>的对比藏在 mixed 数据集内部——<b>连「哪些细胞是对照」都无法从元数据知道</b>。</li>
</ul>
<p><b>★ 可审计性本身与质量混杂（本节最重要的方法学发现）。</b>这 15 个可审计集合的对照臂复用率中位
0.00——单独 deposit 对照臂的研究通常用自有对照——但范围从 0 到 <b>1.0</b>：
其中<b>一个集合的对照臂 100% 复用而病例臂 0%</b>（不对称性 = 1.0，ds4 病理的极端形态），
4/15 存在「对照更差」的不对称。因此我们能查的子集<b>恰好富集干净对照</b>，
而对 89.3% 我们只能说「无法审计」，<b>绝不能说「干净」</b>。
来源检查由此不是一道常规闸门，而是一种<b>特权</b>——最需要它的数据集恰恰是它够不到的。
这也解释了为什么 ds4 的失效必须靠逐数据集读元数据发现，而不可能由一次全库筛查拦下。</p>""")
    A(img("fig_D14_libwide_provenance",
          "图 D14｜CELLxGENE 全库来源流行度（397 集合 / 2,238 数据集 / 292.8M 细胞）。"
          "(A) 每数据集 frac_primary 分布按三类堆叠：mixed（携带对比）数据集大量堆积在 frac_primary≈0。"
          "(B) 15 个臂可分集合的病例臂 vs 对照臂复用率：一个集合对照臂 100% 复用而病例臂 0%；"
          "另有 125 个集合的对比只存在于混合数据集内部（臂不可分）。"
          "(C) 臂不可分比例按组织系统（≥5 个对比集合）：肠 100%、脑 97%、血 91%、肺 91%、肝 80%、肾 67%——"
          "该模式不是肠道特有。(D) 臂可分子集的对照/病例深度比中位 0.86（n=14），远温和于 ds4 的 4.35×。"))
    sys_rows = lw.get("by_tissue_system") or []
    if sys_rows:
        rows = [[b["tissue_system"], b["n_contrast_collections"], b["n_not_arm_separable"],
                 "%.0f%%" % (100 * b["share_not_arm_separable"])] for b in sys_rows]
        A(tbl(["组织系统", "含对比集合数", "臂不可分", "占比"], rows))
    A("""<p class="small">脚本：<code>partD_resfeas_libwide.py</code>（★ 口径修订：mixed 三分，
    修复原「混合数据集整体计入病例臂」导致的对比集合漏计；同时修复 <code>g_ctrl_pure</code>
    在 ctl_only 分支内恒为 False 的潜伏 bug）。明细：<code>partD_resfeas_libwide_datasets.tsv</code> /
    <code>partD_resfeas_libwide_collections.tsv</code> / <code>partD_resfeas_libwide_summary.json</code>。
    全部来自集合详情端点缓存（data/raw/cxg_coll_detail/，397 个 JSON，断点续抓）。</p>
<p><b>对主论文的意义</b>：这一节把 Rung 1 从「4 套数据的个案」升级为「全库流行度 + 一个结构性悖论」——
直接回应「你的来源审计只看了 4 套数据」与「框架是不是 UC 专属」两个审稿质疑；
并成为主论文新 §2.7（<i>Rung 1 at library scale: auditability is confounded with quality</i>）。</p>
</div>""")

A("""
<h2>六、与前三个 Part 的关系：这是第三次同类现象</h2>
<div class="card">
<table><thead><tr><th>部分</th><th>名义上的信号</th><th>去掉混杂后</th><th>混杂类型</th></tr></thead><tbody>
<tr><td>Part A</td><td>44/97 基因 FDR 显著；面板 AUC 0.936</td><td>面板不优于随机基因集；外部验证失败</td><td>批次-疾病混杂（Cramér's V=0.71）</td></tr>
<tr><td>Part C</td><td>9 个基因×组织 FDR 显著</td><td>有效独立工具变量全部 = 1</td><td>距离剪枝假冒 LD clumping</td></tr>
<tr><td>Part B</td><td>（无阳性）</td><td>结论未变：干净阴性</td><td>—</td></tr>
<tr><td>Part D（发现队列 ds1）</td><td>T 细胞 UC 中面板下调 p=0.020</td><td>深度校正后 p=0.189；无任何细胞类型升高</td><td>测序深度-基因集评分混杂</td></tr>
<tr><td>Part D（验证队列 ds4）</td><td>14/25 种细胞类型 UC 中升高，杯状细胞 q=0.016</td><td>抽稀到相同深度后全部消失（最小 q=0.434），效应量收缩 6 倍</td><td>深度与疾病在 donor 级共线（AUC=0.84）</td></tr>
<tr><td>Part D（负对照）</td><td>随机基因集同样「显著」6–11/26、同样 23/26 全向上</td><td>SGRG 全部指标在 200 个随机集分布之内；抽稀后伪影同步归零</td><td>深度混杂 + 置换检验在方差不等时的 I 类错误膨胀（8.5×）</td></tr>
<tr><td><b>Part D（来源审计，本轮）</b></td><td><b>UC/对照深度比 4.3×，被归因为疾病相关</b></td><td><b>把对照组换成研究自己的 2 例后深度比塌缩到 1.04×；ds2/ds3 对照池 98.5%/80.7% 来自无病例的研究</b></td><td><b>对照来源组成 + 疾病与批次完全嵌套</b></td></tr>
<tr><td><b>Part D（检验力裁决，本轮）</b></td><td><b>三种口径三种答案：A 21/26 向上、B 8/23、C 11/15</b></td><td><b>B 口径连炎症通路都测不出（0/23）、C 口径被 p 地板卡住（min p=0.125）→ 均为零检验力</b></td><td><b>设计层面的不可识别 / 无检验力</b></td></tr>
</tbody></table>
<p>四个部分、<b>六种</b>不同的技术混杂形态，指向同一个可写的论点：<b>在没有把技术协变量与对照来源剥离之前，
「基因集/通路被激活」这类结论在小样本公共数据中极不可靠</b>。Part D 补上的正是这个论点在单细胞分辨率的版本，
而且本轮把它的层级又推深了一层：<b>深度校正、随机集负对照这些「下游武器」都无力处理
「对照组本身被外部数据替换」这个上游错误</b>——因为一旦分组错了，所有的校正都是在错误的分组里做得很对。</p>
</div>
""")

# ================= 7 下一步 =================
A("""
<h2>七、下一步（按优先级）</h2>
<div class="card">
<ol>
<li><b>把「对照来源审计」写成方法学贡献的第一主角</b>：本轮的发现按层次排成一个完整的方法学阶梯——
①对照来源（<code>is_primary_data</code> / 研究归属）→ ②批次-疾病嵌套 → ③深度-评分混杂 →
④置换检验 I 类错误膨胀 → ⑤检验力地板。每一层都能用一个零成本检查 + 一个随机集校准来诊断。
四套主流 UC 图谱三套倒在第①②层，这本身就是一个可发表的系统性审计
（<i>Briefings in Bioinformatics</i> / <i>Genome Biology</i> 方法注释级别）。</li>
<li><b>补有检验力的受控数据</b>：结论已经把「需要什么数据」写得非常清楚——
同一研究内 <b>≥15 对</b>「炎症/非炎症配对活检」（同 donor 配对设计的最小可达 p 才能过 0.05），
或者 <b>≥10 例真实健康对照</b>且与病例同批次同流程。候选：GEO 检索已命中的
GSE282124/282123、GSE312415/312420（空间转录组）以及配套的 scRNA 数据；本沙箱 GEO/NCBI 通道实测可达。</li>
<li><b>把 5.13 的检验力曲线与 5.14 的技术轴清单做成可复用的诊断工具</b>：
<code>partD_technical_inventory.py</code> 已经是前一半——对任意一套新数据，
从既有结果表零成本量出「占据度 / 共线度 / 评分-深度耦合」并给出是否会产生名义阳性的预警；
后一半（输入 n_pair 或 n_case/n_ctrl 与效应量，输出最小可达 p 与所需样本量）待做——
审稿人问「为什么不显著」时可以一图作答。</li>
<li><b>论文定位随之更新</b>：Part D 不再是「单细胞分辨率的阴性验证」，而是
「<b>单细胞公共资源的可用性审计</b>」——用 SGRG 面板做一个完整的案例研究，
展示从名义阳性到检验力地板的完整崩溃链条。这与出路 A 的方法学主线是同一篇论文的两个层面。</li>
</ol>
</div>
""")

# ================= 8 局限 =================
A("""
<h2>八、局限（当前结论不能外推的部分）</h2>
<div class="card">
<ul>
<li><b>「无检验力」与「无效应」仍然不能互换</b>：本轮最重要的修订是把此前的阴性结论改写为
「在这些设计下无法检出」。但这不等于说面板真的被激活了——5.11 的随机集证据依然成立：
在能给出显著结果的口径里，SGRG 的表现与随机基因集无法区分。</li>
<li><b><code>is_primary_data</code> 的语义边界</b>：它标记的是「是否属于本集合的主数据」，
不是「是否属于这项研究的对照」。5 个复用 donor（411c–440c）有可能确实是正常肠组织的合法参考
（比如同一联盟的其他研究），把它们完全剔除是一种保守选择；但它们与病例组在批次、化学版本、
深度上系统性不同，保留它们只会引入更严重的混杂。两条路都走了一遍（5.13 的 A 与 B），结论一致。</li>
<li><b>同 donor 炎症标注的来源未核验</b>：<code>inflammation</code> 字段在 ds4 里是细胞级的，
且同一样本内可混有两种取值；原文方法部分未能取得。若它是由某种转录组评分派生的，
则 C 口径存在<b>循环论证风险</b>——缓解证据是五个 Hallmark 炎症通路在 13–15/15 种细胞类型方向一致升高
（若标签只是随机分层，不会出现这种方向性），但这不排除它与 SGRG 面板存在共同基因。</li>
<li><b>抽稀的代价</b>：深度匹配丢弃了 25% 的低深度细胞。若这 25% 有选择性，抽稀会引入新偏倚。</li>
<li><b>打分口径</b>：面板评分衡量「转录本丰度」，不是应激颗粒的实际装配状态（翻译后与相分离事件）。</li>
<li><b>ds1 的 6 donor/组</b>：pseudobulk 置换检验的最小可达 p≈0.002，对 d≈0.3 的效应检出力有限。</li>
</ul>
</div>
""")

# ================= 9 文件 =================
files = [
    ["scripts/partD_common.py", "公共模块：h5ad 读取（自动判别 counts/log 标准化）、归一化、score_genes、置换检验、BH"],
    ["scripts/partD_analyse_ds1.py", "ds1 主分析（per-cell 打分 + pseudobulk + 检验）"],
    ["scripts/partD_ds1_genes_diagnostics.py", "ds1 逐基因 pseudobulk + 深度混杂诊断 + 深度校正检验"],
    ["scripts/partD_ds1_state_analysis.py", "ds1 炎症状态五分位耦合分析"],
    ["scripts/partD_analyse_generic.py", "通用分析器（ds2/ds3 共用，带组织/疾病过滤）"],
    ["scripts/partD_analyse_stream.py", "分块流式分析器（ds4 32 万细胞：h5py 按行分块 + 两遍扫描 + 可选深度匹配抽稀）"],
    ["scripts/partD_fig_crossdataset.py", "图 D6 跨数据集森林图"],
    ["scripts/partD_fig_ds4_depthcontrol.py", "图 D7 ds4 三口径对照图"],
    ["scripts/partD_random_control.py", "随机基因集负对照（表达分层匹配 × 200，一次矩阵乘并行打分，双检验口径）"],
    ["scripts/partD_fig_randomcontrol.py", "图 D8 负对照零分布图"],
    ["scripts/partD_fig_randcalib.py", "图 D9 深度匹配前后校准对比图"],
    ["scripts/partD_batch_audit.py", "<b>来源审计（本轮）</b>：四套数据的批次-疾病嵌套、对照同源性、细胞级深度"],
    ["scripts/partD_provenance_analysis.py", "<b>按来源重做（本轮）</b>：A 朴素 / B 同源对照 / C 同 donor 炎症配对，同一遍扫描 + 200 随机集校准"],
    ["scripts/partD_fig_provenance.py", "图 D10 对照来源与深度结构"],
    ["scripts/partD_fig_power.py", "图 D11 三口径检验力裁决"],
    ["scripts/partD_technical_inventory.py", "<b>步 1 技术轴清单（本轮）</b>：A2 队列级宽表的等构对应物——四套数据的深度/批次/主轴占据度/共线度/干预后显著数，只读既有结果表，不重跑 h5ad"],
    ["scripts/partD_fig_inventory.py", "图 AD4 技术轴清单三面板（深度口径敏感性 / 占据×共线象限 / ds4 四道干预）"],
    ["scripts/partD_get_genesets.py", "MSigDB Hallmark + QuickGO 拆分 SGRG 亚组"],
    ["scripts/probe_sc_sources.py / probe_sc_channels2.py / probe_cxg_uc_detail.py", "数据通道与数据集核验"],
    ["scripts/dl_h5ad.py / dl_parallel.py", "h5ad 下载（断点续传 / 8 线程分片）"],
    ["results/fig_D1..D11 (.png/.pdf)", "图件（600 dpi，双格式）"],
    ["results/partD_ds1_* ~ partD_ds4_tearle*", "四套数据集的全部结果表（per-cell 打分、pseudobulk、检验、逐基因、诊断）"],
    ["results/partD_ds4_tearle_rand* / partD_ds4_tearle_ds1k_rand*", "负对照：200 个随机基因集的检验表、汇总、pseudobulk"],
    ["results/partD_batch_audit_*.tsv", "<b>来源审计结果（本轮）</b>：donor 级深度/批次表 + 汇总指标"],
    ["results/partD_ds4_prov_*.tsv", "<b>按来源重做的结果（本轮）</b>：A/B 口径检验表、C 口径配对检验、per-cell 打分（pickle）"],
    ["results/partD_technical_inventory.tsv / _meta.json", "<b>步 1 技术轴清单（本轮）</b>：4 数据集 × 51 列，含口径定义与稳健性警示"],
    ["results/fig_AD4_technical_inventory.png/pdf", "图 AD4（方法学合并稿主图之一）"],
]
A("<h2>九、文件清单</h2><div class='card'>" +
  tbl(["文件", "说明"], [[f"<code>{a}</code>", b] for a, b in files]) +
  "</div>")

A("</div></body></html>")

out = os.path.join(OUT, "PartD_单细胞方案启动报告.html")
open(out, "w", encoding="utf-8").write("".join(H))
print("[saved]", out, "%.2f MB" % (os.path.getsize(out) / 1e6))
