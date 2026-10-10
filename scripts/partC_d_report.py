#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Part C 报告生成器 —— 自包含 HTML（图片 base64 内嵌），中文字体安全。"""
import os
import sys
import csv
import json
import math
import base64

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

for _f in ["C:/Windows/Fonts/msyh.ttc", "C:/Windows/Fonts/msyhbd.ttc",
           "C:/Windows/Fonts/simhei.ttf"]:
    if os.path.exists(_f):
        try:
            font_manager.fontManager.addfont(_f)
        except Exception:
            pass
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
RES = os.path.join(ROOT, "results")
J = json.load(open(os.path.join(RES, "partC_mr_summary.json"), encoding="utf-8"))
SETS = J["sets"]
SUM = J["summary"]
COLOC = J["coloc"]
DETAIL = J["results"]
NEI = list(csv.DictReader(open(os.path.join(RES, "partC_neighbour_eqtl.tsv"),
                               encoding="utf-8"), delimiter="\t"))

TISSUE_COLOR = {"Colon_Sigmoid": "#0b5cab", "Colon_Transverse": "#3d8bd4",
                "Brain_Cortex": "#8a5a00", "Whole_Blood": "#17643b",
                "Brain_Hippocampus": "#b42318"}
TISSUE_ZH = {"Colon_Sigmoid": "结肠乙状", "Colon_Transverse": "结肠横",
             "Brain_Cortex": "脑皮层", "Whole_Blood": "全血",
             "Brain_Hippocampus": "海马"}


def fnum(x):
    try:
        if x in ("", None):
            return None
        v = float(x)
        if v != v:            # NaN
            return None
        return v
    except Exception:
        return None


def n(x, d=3):
    v = fnum(x)
    return "—" if v is None else ("%.{0}f".format(d) % v)


def sci(x, d=2):
    v = fnum(x)
    if v is None:
        return "—"
    if v == 0:
        return "0"
    if v < 1e-2 or v >= 1e4:
        return ("%.{0}e".format(d) % v)
    return ("%.{0}f".format(d) % v)


def psci(x, d=2):
    """p-value formatter: an exact zero means numerical underflow
    (double precision floor), not a true p of 0."""
    v = fnum(x)
    if v is None:
        return "—"
    if v == 0:
        return "&lt;1e-300"
    return sci(x, d)


def img(name, cap=""):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        return ""
    b = base64.b64encode(open(p, "rb").read()).decode()
    c = "<figcaption>%s</figcaption>" % cap if cap else ""
    return ("<figure><img src='data:image/png;base64,%s'>%s</figure>" % (b, c))


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def fig_instruments():
    rows = sorted(SUM, key=lambda r: r["IVW_p"])
    labels = ["%s·%s" % (r["gene"], TISSUE_ZH.get(r["tissue"], r["tissue"])) for r in rows]
    k = [r["nsnp"] for r in rows]
    cols = [TISSUE_COLOR.get(r["tissue"], "#888") for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(13, max(4.2, 0.24 * len(rows))),
                             gridspec_kw={"width_ratios": [1, 1]})
    y = np.arange(len(rows))
    axes[0].barh(y, k, color=cols)
    axes[0].set_yticks(y); axes[0].set_yticklabels(labels, fontsize=7)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("独立工具变量个数 k（10 kb 距离剪枝）", fontsize=9)
    axes[0].set_title("A 每套 cis-eQTL 工具变量数量", fontsize=10)
    axes[0].grid(alpha=.25, axis="x")
    F = [r["meanF"] for r in rows]
    axes[1].barh(y, F, color=cols)
    axes[1].axvline(10, color="#b42318", ls="--", lw=1.2)
    axes[1].text(10.5, len(rows) - 0.5, "F=10 弱工具阈值", color="#b42318", fontsize=8)
    axes[1].set_yticks(y); axes[1].set_yticklabels([])
    axes[1].invert_yaxis()
    axes[1].set_xlabel("平均 F 统计量", fontsize=9)
    axes[1].set_title("B 工具变量强度", fontsize=10)
    axes[1].grid(alpha=.25, axis="x")
    hs = [plt.Line2D([0], [0], color=c, lw=6) for c in
          [TISSUE_COLOR[t] for t in ["Colon_Sigmoid", "Colon_Transverse",
                                     "Brain_Cortex", "Whole_Blood"]]]
    axes[1].legend(hs, ["结肠乙状", "结肠横", "脑皮层", "全血"], fontsize=8, loc="lower right")
    plt.tight_layout()
    plt.savefig(os.path.join(RES, "fig_partC_instruments.png"), dpi=190)
    plt.close()


def fig_forest():
    rows = sorted(SUM, key=lambda r: r["IVW_b"])
    labels = ["%s·%s (k=%d)" % (r["gene"], TISSUE_ZH.get(r["tissue"], r["tissue"]), r["nsnp"])
              for r in rows]
    fig, ax = plt.subplots(figsize=(9.5, max(5, 0.27 * len(rows))))
    y = np.arange(len(rows))
    for i, r in enumerate(rows):
        b, se = r["IVW_b"], r["IVW_se"]
        lo, hi = b - 1.96 * se, b + 1.96 * se
        c = TISSUE_COLOR.get(r["tissue"], "#888")
        sig = r["IVW_p"] < 0.05
        ax.plot([lo, hi], [i, i], color=c, lw=2.1 if sig else 1.3,
                alpha=1.0 if sig else .55)
        ax.plot([b], [i], "o", color=c, ms=6 if sig else 4,
                mfc=c if sig else "white", mec=c)
    ax.axvline(0, color="#1a1f27", lw=1.1)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("IVW 因果估计（每 1 个 SD 表达量变化对应的 UC log-OR）", fontsize=9.5)
    ax.set_title("SGRG cis-eQTL → 溃疡性结肠炎：IVW 因果效应森林图\n"
                 "（点为估计值，横线为 95% CI；实心粗线 = p<0.05）", fontsize=10.5)
    ax.grid(alpha=.25, axis="x")
    hs = [plt.Line2D([0], [0], color=TISSUE_COLOR[t], lw=6) for t in
          ["Colon_Sigmoid", "Colon_Transverse", "Brain_Cortex", "Whole_Blood"]]
    ax.legend(hs, ["结肠乙状", "结肠横", "脑皮层", "全血"], fontsize=8, loc="lower right")
    plt.tight_layout()
    plt.savefig(os.path.join(RES, "fig_partC_forest.png"), dpi=190)
    plt.close()


def fig_scatter(top=4):
    rows = sorted(SUM, key=lambda r: r["IVW_p"])[:top]
    if not rows:
        return
    fig, axes = plt.subplots(1, len(rows), figsize=(3.7 * len(rows), 3.6))
    if len(rows) == 1:
        axes = [axes]
    for ax, r in zip(axes, rows):
        key = "%s|%s" % (r["gene"], r["tissue"])
        d = DETAIL.get(key)
        if not d:
            ax.axis("off"); continue
        bx = np.array(d["bx"]); by = np.array(d["by"])
        sy = np.array(d["sy"])
        ax.errorbar(bx, by, yerr=1.96 * sy, fmt="none", ecolor="#c9d3de", lw=.9)
        ax.scatter(bx, by, s=26, color=TISSUE_COLOR.get(r["tissue"], "#888"), zorder=3)
        xs = np.linspace(bx.min(), bx.max(), 50)
        ax.plot(xs, xs * r["IVW_b"], color="#0b5cab", lw=1.6,
                label="IVW β=%.3f" % r["IVW_b"])
        eg = [m for m in [r.get("Egger_b")] if m not in ("", None)]
        if eg:
            from scipy import stats as _s
            try:
                ax.plot(xs, xs * float(eg[0]), color="#8a5a00", lw=1.4, ls="--",
                        label="MR-Egger β=%.3f" % float(eg[0]))
            except Exception:
                pass
        ax.axhline(0, color="#c9d3de", lw=.8); ax.axvline(0, color="#c9d3de", lw=.8)
        ax.set_title("%s · %s\nIVW p=%s" % (r["gene"], TISSUE_ZH.get(r["tissue"], r["tissue"]),
                                            sci(r["IVW_p"])), fontsize=9)
        ax.set_xlabel("eQTL β（表达量）", fontsize=8.5)
        ax.set_ylabel("UC GWAS β（log-OR）", fontsize=8.5)
        ax.legend(fontsize=7.5); ax.grid(alpha=.22)
    plt.tight_layout()
    plt.savefig(os.path.join(RES, "fig_partC_scatter.png"), dpi=190)
    plt.close()


def fig_coloc(top=18):
    rows = sorted([c for c in COLOC], key=lambda c: -c["PP.H4"])[:top]
    if not rows:
        return
    labels = ["%s·%s" % (c["gene"], TISSUE_ZH.get(c["tissue"], c["tissue"])) for c in rows]
    y = np.arange(len(rows))
    comps = [("PP.H4", "#b42318", "H4 共享因果"), ("PP.H3", "#8a5a00", "H3 不同因果"),
             ("PP.H1", "#7aa7d6", "H1 仅 eQTL"), ("PP.H2", "#a8c6a0", "H2 仅 GWAS")]
    left = np.zeros(len(rows))
    fig, ax = plt.subplots(figsize=(9.5, max(4, 0.3 * len(rows))))
    for key, col, lab in comps:
        v = np.array([c.get(key, 0) for c in rows])
        ax.barh(y, v, left=left, color=col, label=lab)
        left += v
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis(); ax.set_xlim(0, 1)
    ax.set_xlabel("后验概率", fontsize=9)
    ax.set_title("cis 区域共定位（coloc.abf，PriPP）", fontsize=10.5)
    ax.legend(fontsize=8, ncol=4, loc="lower right")
    for i, c in enumerate(rows):
        if c["PP.H4"] > 0.3:
            ax.text(1.005, i, "%.2f" % c["PP.H4"], fontsize=7.5, va="center",
                    color="#b42318")
    plt.tight_layout()
    plt.savefig(os.path.join(RES, "fig_partC_coloc.png"), dpi=190)
    plt.close()


def fig_sensitivity():
    rows = sorted(SUM, key=lambda r: r["IVW_p"])[:24]
    if not rows:
        return
    labels = ["%s·%s" % (r["gene"], TISSUE_ZH.get(r["tissue"], r["tissue"])) for r in rows]
    mats = []
    names = ["IVW p", "加权中位数 p", "MR-Egger p", "Egger 截距 p",
             "Cochran Q p", "MR-PRESSO 全局 p", "Steiger 方向正确"]
    for r in rows:
        st = SETS["%s|%s" % (r["gene"], r["tissue"])]
        d = DETAIL.get("%s|%s" % (r["gene"], r["tissue"]), {})
        ivw = [m for m in d.get("methods", {}).values()] if False else None
        # gather from summary row
        row = [r["IVW_p"], r.get("WMed_p", ""), r.get("Egger_p", ""),
               r.get("Egger_int_p", ""), r.get("IVW_Q_p", ""),
               r.get("PRESSO_global_p", ""), ""]
        mats.append(row)
    fig, ax = plt.subplots(figsize=(6.6, max(4, 0.3 * len(rows))))
    grid = np.full((len(rows), len(names)), np.nan)
    text = [["" for _ in names] for _ in rows]
    for i, row in enumerate(mats):
        for j, v in enumerate(row):
            if j == 6:
                continue
            fv = fnum(v)
            if fv is None:
                text[i][j] = "—"; continue
            grid[i, j] = fv
            text[i][j] = "<1e-300" if fv == 0 else sci(fv, 1)
    im = ax.imshow(np.where(np.isnan(grid), np.nan, np.clip(-np.log10(np.maximum(grid, 1e-300)), 0, 12)),
                   cmap="RdYlGn_r", aspect="auto", vmin=0, vmax=3)
    for i in range(len(rows)):
        for j in range(len(names) - 1):
            ax.text(j, i, text[i][j], ha="center", va="center", fontsize=7)
        ok = [r for r in rows][i]["steiger_correct"]
        ax.text(len(names) - 1, i, "OK" if ok else "NO", ha="center", va="center",
                fontsize=8.5, color="#17643b" if ok else "#b42318")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, fontsize=8, rotation=28, ha="right")
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_title("敏感性分析一览（数值为 p 值；绿色=不显著/无异常，红色=显著）", fontsize=10)
    plt.colorbar(im, ax=ax, label="-log10(p)", fraction=.03)
    plt.tight_layout()
    plt.savefig(os.path.join(RES, "fig_partC_sensitivity.png"), dpi=190)
    plt.close()


def fig_neighbour():
    rows = sorted(NEI, key=lambda r: -float(r["delta_log10p"]))
    if not rows:
        return
    labels = ["%s·%s" % (r["gene"], TISSUE_ZH.get(r["tissue"], r["tissue"])) for r in rows]
    d = [float(r["delta_log10p"]) for r in rows]
    cols = ["#b42318" if r["flag"] == "True" else "#7aa7d6" for r in rows]
    fig, ax = plt.subplots(figsize=(9.0, max(4.2, 0.26 * len(rows))))
    y = np.arange(len(rows))
    ax.barh(y, d, color=cols)
    for i, r in enumerate(rows):
        ax.text(d[i] + (1.5 if d[i] >= 0 else -1.5), i,
                "%s (%.0f kb)" % (r["neighbour"], float(r["dist_kb"])),
                va="center", ha="left" if d[i] >= 0 else "right", fontsize=6.5,
                color="#5b6470")
    ax.axvline(5, color="#8a5a00", ls="--", lw=1.2)
    ax.axvline(0, color="#1a1f27", lw=1)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7.2)
    ax.invert_yaxis()
    ax.set_xlim(min(-10, min(d) - 8), max(d) + 42)
    ax.set_xlabel("log10(p_目标基因) - log10(p_最强邻居)   >0 表示邻居 eQTL 更强", fontsize=9)
    ax.set_title("cis 区域 LD 混杂诊断：邻居基因 eQTL 优势度\n"
                 "红 = 邻居在 250 kb 内且强 5 个数量级以上（高度疑似 LD 混杂）", fontsize=10)
    ax.grid(alpha=.22, axis="x")
    plt.tight_layout()
    plt.savefig(os.path.join(RES, "fig_partC_neighbour.png"), dpi=190)
    plt.close()


def main():
    fig_instruments()
    fig_forest()
    fig_scatter()
    fig_neighbour()
    fig_coloc()
    fig_sensitivity()
    print("figures done")

    # ---------------- narrative stats ----------------
    sig05 = [r for r in SUM if r["IVW_p"] < 0.05]
    sigfdr = [r for r in SUM if r["IVW_q"] < 0.05]
    coloc_strong = [c for c in COLOC if c["PP.H4"] > 0.8]
    coloc_mid = [c for c in COLOC if 0.5 < c["PP.H4"] <= 0.8]
    n_gene = len(set(r["gene"] for r in SUM))
    ext_sets = [r for r in SUM if r["gene"] in ("EIF2AK4", "EIF2AK2", "EIF2AK3", "G3BP1")]
    n_ext_sets = len(ext_sets)
    n_ext_genes = len(set(r["gene"] for r in ext_sets))
    n_sets_stage_a = len(SETS)
    n_dropped = n_sets_stage_a - len(SUM)
    nei_by = {"%s|%s" % (r["gene"], r["tissue"]): r for r in NEI}
    flagged = [r for r in SUM if nei_by.get("%s|%s" % (r["gene"], r["tissue"]),
                                            {}).get("flag") == "True"]
    sigfdr_sets = set("%s|%s" % (r["gene"], r["tissue"]) for r in sigfdr)
    flagged_sig = [f for f in flagged
                   if "%s|%s" % (f["gene"], f["tissue"]) in sigfdr_sets]
    clean_sig = [r for r in sigfdr
                 if nei_by.get("%s|%s" % (r["gene"], r["tissue"]), {}).get("flag") != "True"]

    def rowtab(rows, extra=False):
        out = []
        for r in rows:
            out.append(
                "<tr><td>%s</td><td>%s</td><td class='n'>%d</td><td class='n'>%s</td>"
                "<td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td>"
                "<td class='n'>%s</td><td class='n'>%s</td><td>%s</td></tr>" % (
                    r["gene"], TISSUE_ZH.get(r["tissue"], r["tissue"]), r["nsnp"],
                    n(r["IVW_b"]), n(r["IVW_se"]), psci(r["IVW_p"]), sci(r["IVW_q"]),
                    psci(r.get("WMed_p")), psci(r.get("Egger_p")),
                    "✓" if r["steiger_correct"] else "✗"))
        return "\n".join(out)

    H = []
    A = H.append
    A("""<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">
<title>Part C｜SGRG × 溃疡性结肠炎 孟德尔随机化因果推断</title>
<style>
:root{--ink:#1a1f27;--muted:#5b6470;--line:#e2e7ee;--bg:#f7f9fc;--card:#fff;
--red:#b42318;--blue:#0b5cab;--amber:#8a5a00;--green:#17643b;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:16px/1.75 -apple-system,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif}
.wrap{max-width:1010px;margin:0 auto;padding:40px 26px 80px}
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
table{width:100%;border-collapse:collapse;margin:14px 0;font-size:13px;background:#fff}
th,td{border:1px solid var(--line);padding:6px 8px;text-align:left;vertical-align:top}
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

    A("<h1>Part C｜应激颗粒相关基因 × 溃疡性结肠炎：孟德尔随机化因果推断</h1>")
    A("<div class='sub'>暴露：GTEx v8 cis-eQTL（结肠乙状 / 结肠横 / 脑皮层 / 全血）　"
      "结局：de Lange 2017 UC GWAS（GCST004133，12,366 例 / 33,609 对照，N=45,975）<br>"
      "匹配：rsID（暴露 GRCh38 / 结局 GRCh37，无需 liftOver）　"
      "方法：IVW / MR-Egger / 加权中位数 / 加权众数 + Cochran Q / MR-PRESSO / LOO / Steiger + coloc.abf<br>"
      "生成时间：2026-10-08</div>")

    A("<div class='kpi'>"
      "<div><b>%d</b><span>可构建工具变量的（基因×组织）组合</span></div>"
      "<div><b>%d</b><span>涉及 SGRG 基因数</span></div>"
      "<div><b>%d</b><span>IVW 名义显著（p&lt;0.05）</span></div>"
      "<div><b>%d</b><span>IVW 通过 FDR&lt;0.05</span></div>"
      "<div><b>%d</b><span>LD 混杂高度疑似（含 %d 个 FDR 阳性）</span></div>"
      "</div>" % (len(SUM), n_gene, len(sig05), len(sigfdr), len(flagged), len(flagged_sig)))

    A("<div class='card'><h3>一句话结论</h3><p>"
      "在 %d 套（基因×组织）cis-eQTL 工具变量中，<b>%d 套</b>的 IVW 因果估计通过 FDR&lt;0.05，"
      "集中在<b>全血</b>（ATXN2L、CAPRIN1、RBPMS2、OGFOD1、LARP1、PABPC4、DDX1）与结肠（DDX25）；"
      "但其中 <b>%d 套</b>在 LD 混杂诊断中被标记为高度疑似——其 cis 窗口内 250 kb 以内存在 eQTL "
      "强出 5 个数量级以上的邻居基因（最典型：ATXN2L 的引物 SNP 距 TUFM 的 eQTL 峰仅 31 kb，"
      "而 TUFM 的 p 值强出 140 个数量级）。"
      "<b>Part A 中最稳健的转录组候选（VCP、ZFAND1、ATXN2、QKI、TRIM21）没有一个获得因果支持</b>"
      "（VCP 全血 p=0.71、ZFAND1 脑皮层 p=0.34、QKI 全血 p=0.66）。"
      "因此本文<b>不应把任何 SGRG 基因的表达—UC 关联升级为因果</b>；孟德尔随机化在本数据中的"
      "作用是<b>否定性证据</b>：它系统性地削弱了 Part A 的表观关联。"
      "</p></div>" % (len(SUM), len(sigfdr), len(flagged_sig)))

    # ---- 1 设计 ----
    A("<h2>一、研究设计与数据来源</h2>")
    A("<div class='card'>"
      "<table><tr><th>要素</th><th>选择</th><th>说明</th></tr>"
      "<tr><td>暴露</td><td>GTEx v8 cis-eQTL</td><td>TSS±1 Mb；组织取结肠乙状（n=318）、"
      "结肠横（n=368）、脑皮层（n=205）、全血（n=670）；<code>/association/singleTissueEqtlByLocation</code></td></tr>"
      "<tr><td>暴露效应量</td><td>NES（归一化效应量）</td><td>SE 由 GTEx 报告的 t 统计量反推："
      "<code>se = |NES| / |z|</code>，<code>z = Φ⁻¹(1−p/2)</code></td></tr>"
      "<tr><td>结局</td><td>de Lange 2017 UC（GCST004133）</td><td>12,366 例 / 33,609 对照；"
      "GWAS Catalog <i>harmonised</i> 文件（含 <code>hm_rsid</code>、等位基因、EAF）</td></tr>"
      "<tr><td>SNP 匹配</td><td>rsID</td><td>暴露为 GRCh38、结局为 GRCh37 → 按 rsID 匹配，"
      "规避 liftOver 引入的错配</td></tr>"
      "<tr><td>等位基因协调</td><td>正向链直接比对</td><td>GTEx <code>variantId</code> 与 GWAS Catalog "
      "harmonised 均为正向链 REF/ALT，等位基因完全相同时方向无歧义；互补链由字母补码识别并翻转 β 符号</td></tr>"
      "<tr><td>独立性</td><td>10 kb 距离剪枝</td><td>本环境无 LD 参考面板（Ensembl LD 接口不可达、"
      "无 plink/tabix），改用保守的距离剪枝近似 r²&lt;0.001</td></tr>"
      "<tr><td>强度</td><td>平均 F&gt;10</td><td>全部 %d 套均通过（最低 %.1f）</td></tr>"
      "</table></div>" % (len(SUM), min(r["meanF"] for r in SUM)))

    A("<div class='note'><h4>基因集边界（必须写明）</h4>"
      "<p>主分析使用严格注释基因集 <b>GO:0010494（细胞质应激颗粒）∪ GO:0035061（核应激颗粒）</b>，"
      "经 HGNC 校验后保留 104 个官方符号。该集合<b>不含 EIF2AK4</b>——EIF2AK4（GCN2）是磷酸化 eIF2α、"
      "<i>启动</i>颗粒组装的激酶，本身不是颗粒组分，故未被 GO 细胞组分项注释。</p>"
      "<p>考虑到项目前置核验报告已把 EIF2AK4 列为首选候选，本文将其与 EIF2AK2/EIF2AK3（及仅脑皮层有 eQTL 信号的 G3BP1）"
      "一起作为<b>事先声明的「eIF2α 激酶通路扩展面板」</b>单列呈现——"
      "该面板最终进入 MR 的是 <b>%d 个基因 / %d 套（基因×组织）组合</b>（%s），"
      "全程与主面板分开报告，不并入 SGRG 结论。</p></div>"
      % (n_ext_genes, n_ext_sets,
         "、".join(sorted(set(r["gene"] for r in ext_sets))) or "无"))

    # ---- 2 工具变量 ----
    A("<h2>二、工具变量</h2>")
    A("<p>筛出 <b>%d</b> 套（基因×组织）具有全基因组显著 cis-eQTL 的组合；其中 <b>%d</b> 套"
      "的引物 SNP 在结局 GWAS 中全部缺失（rsID 未收录）而无法完成 MR，最终 <b>%d</b> 套进入分析，"
      "覆盖 <b>%d</b> 个基因、<b>%d</b> 个经等位基因协调后可用的独立 SNP"
      "（另有 %d 个完美 LD 代理 SNP 在工具变量构建阶段被剔除）。"
      "全部组合的平均 F 统计量 ≥ %.1f，最低者为 %s·%s（F=%.1f），远超弱工具阈值 F=10。</p>"
      % (n_sets_stage_a, n_dropped, len(SUM), n_gene,
         sum(r["nsnp"] for r in SUM),
         int(sum(v.get("n_proxy_dropped", 0) for v in SETS.values())),
         min(r["meanF"] for r in SUM),
         min(SUM, key=lambda r: r["meanF"])["gene"],
         TISSUE_ZH.get(min(SUM, key=lambda r: r["meanF"])["tissue"], ""),
         min(r["meanF"] for r in SUM)))
    A(img("fig_partC_instruments.png",
          "图 1　工具变量的数量（A）与强度（B）。颜色区分组织。"))
    A("<div class='card'><table><tr><th>基因</th><th>组织</th><th class='n'>k</th>"
      "<th class='n'>平均 F</th><th class='n'>最小 P</th><th>首席 SNP</th></tr>"
      + "".join("<tr><td>%s</td><td>%s</td><td class='n'>%d</td><td class='n'>%.1f</td>"
                "<td class='n'>%s</td><td>%s</td></tr>"
                % (v["gene"], TISSUE_ZH.get(v["tissue"], v["tissue"]), v["n_instruments"],
                   v["meanF"], sci(v["min_p"]), v["lead_snp"])
                for k, v in sorted(SETS.items(), key=lambda kv: kv[1]["min_p"])[:20])
      + "</table><p style='color:#5b6470;font-size:12.5px'>仅列最小 P 最小的 20 套，"
        "完整表见 <code>partC_instrument_sets.json</code>／<code>partC_instruments.tsv</code>。</p></div>")

    # ---- 3 MR 主结果 ----
    A("<h2>三、两样本 MR 主结果</h2>")
    A(img("fig_partC_forest.png",
          "图 2　IVW 因果效应森林图（因果估计以 UC log-OR / 每 SD 表达量表示）。"))
    A("<div class='card'><table><tr><th>基因</th><th>组织</th><th class='n'>k</th>"
      "<th class='n'>IVW β</th><th class='n'>SE</th><th class='n'>IVW P</th>"
      "<th class='n'>FDR Q</th><th class='n'>加权中位数 P</th><th class='n'>Egger P</th>"
      "<th>Steiger</th></tr>" + rowtab(sorted(SUM, key=lambda r: r["IVW_p"]))
      + "</table></div>")

    # ---- 4 敏感性 ----
    A("<h2>四、敏感性分析与稳健性</h2>")
    A(img("fig_partC_sensitivity.png",
          "图 3　敏感性分析一览：各检验 p 值（数值为 p，色阶为 −log10 p）。"
          "最右列为 Steiger 方向检验是否通过。"))
    A(img("fig_partC_scatter.png",
          "图 4　工具变量-结局效应散点（按 IVW P 最优的 4 套）。实线为 IVW 斜率，虚线为 MR-Egger 斜率。"))
    A("<div class='card'><table><tr><th>基因</th><th>组织</th><th class='n'>Cochran Q</th>"
      "<th class='n'>Q 的 P</th><th class='n'>Egger 截距 P</th><th class='n'>MR-PRESSO 全局 P</th>"
      "<th class='n'>离群点</th><th class='n'>校正后 β</th><th class='n'>校正后 P</th></tr>"
      + "".join("<tr><td>%s</td><td>%s</td><td class='n'>%s</td><td class='n'>%s</td>"
                "<td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td>"
                "<td class='n'>%s</td><td class='n'>%s</td></tr>"
                % (r["gene"], TISSUE_ZH.get(r["tissue"], r["tissue"]), sci(r.get("IVW_Q"), 1),
                   psci(r.get("IVW_Q_p")), psci(r.get("Egger_int_p")),
                   psci(r.get("PRESSO_global_p")), r.get("PRESSO_n_outlier", "—"),
                   n(r.get("PRESSO_corrected_b")), psci(r.get("PRESSO_corrected_p")))
                for r in sorted(SUM, key=lambda x: x["IVW_p"])[:20])
      + "</table><p style='color:#5b6470;font-size:12.5px'>k=1 的组合无法做异质性/Egger 检验，"
        "相应单元格为空。表中仅列 IVW P 最小的 20 套。</p>"
        "<div class='warn' style='margin-top:14px'><h4>不要采信加权中位数/众数的极端小 p 值</h4>"
        "<p>加权中位数与加权众数的标准误来自 bootstrap。当工具变量彼此高度相关（同处于一个强 LD 区块）时，"
        "重抽样得到的比值几乎不变，bootstrap SE 会<b>塌缩到接近 0</b>，从而产生「伪精确」的小 p 值——"
        "表中 ATXN2L 的加权中位数 p=3.7×10⁻¹⁷⁰、MR-Egger p=4.8×10⁻¹⁰¹ 都属于这一类。"
        "它们不能读作「证据极强」，反而再次提示该组合的工具变量<b>不独立</b>。"
        "IVW 的 p 值同样会因残余 LD 而被低估（见第七节局限）。</p></div></div>")

    # ---- 5 LD contamination diagnostic ----
    A("<h2>五、LD 混杂诊断：为什么这些 MR 阳性信号不能直接采信</h2>")
    A("<div class='warn'><h4>先说一个方法学硬约束（决定本节的做法）</h4>"
      "<p>本环境能访问的 GTEx 接口是 <code>/association/singleTissueEqtlByLocation</code>，"
      "其官方名称即 <b>“Get <u>Significant</u> Single Tissue eQTLs”</b>——它<b>只返回已过显著阈值的"
      "变异-基因对</b>（实测：VCP 窗口内返回的全部 76 个 SNP 的 p 值都 &lt; 2×10⁻⁴）。"
      "共定位（coloc.abf）需要<b>区域内全部 SNP</b>（含大量零效应 SNP 构成的背景），"
      "因此<b>本数据无法给出有效的 coloc 后验概率</b>。"
      "（完整 all-pairs 数据每组织约 4 GB，公共桶内无 tabix 索引，v8 all-pairs 位于 requester-pays 桶。）</p>"
      "<p>替代方案：用一个<b>只用同一批数据即可严格计算</b>的指标——"
      "<b>邻居 eQTL 优势度</b>。cis-eQTL 工具变量只有在“目标基因的引物 SNP 不是某个邻居基因更强 eQTL 的"
      "标记”时才可解释；若 250 kb 内存在强出 ≥5 个数量级的邻居 eQTL，则该 MR 估计更可能来自"
      "<b>水平多效性 / 连锁不平衡污染</b>，而非目标基因本身。</p></div>")
    A(img("fig_partC_neighbour.png",
          "图 5　LD 混杂诊断。横轴为 log₁₀(p_目标) − log₁₀(p_最强邻居)；>0 表示邻居 eQTL 更强。"
          "红色为同时满足「距离 ≤250 kb」且「强 ≥5 个数量级」的高风险组合。"))
    A("<div class='card'><table><tr><th>基因</th><th>组织</th><th class='n'>目标 lead p</th>"
      "<th>最强邻居</th><th class='n'>邻居 lead p</th><th class='n'>间距 (kb)</th>"
      "<th class='n'>Δlog₁₀P</th><th>已知 IBD 基因</th><th>判定</th></tr>"
      + "".join("<tr><td>%s</td><td>%s</td><td class='n'>%s</td><td>%s</td>"
                "<td class='n'>%s</td><td class='n'>%.0f</td><td class='n'>%.1f</td>"
                "<td>%s</td><td>%s</td></tr>"
                % (r["gene"], TISSUE_ZH.get(r["tissue"], r["tissue"]),
                   sci(r["lead_p"]), r["neighbour"], sci(r["nb_p"]),
                   float(r["dist_kb"]), float(r["delta_log10p"]),
                   "✓" if r["nb_is_known_ibd"] == "True" else "",
                   "<b class='red'>高度疑似混杂</b>" if r["flag"] == "True" else "—")
                for r in sorted(NEI, key=lambda x: -float(x["delta_log10p"]))[:24])
      + "</table><p style='color:#5b6470;font-size:12.5px'>按 Δlog₁₀P 降序，显示前 24 套"
        "（完整 48 套见 <code>partC_neighbour_eqtl.tsv</code>）。</p></div>")

    A("<div class='note'><h4>最值得警惕的一例：ATXN2L（全血）</h4>"
      "<p>ATXN2L 是全部组合中 MR 估计最强的一个（IVW β=−0.677，p=2.8×10⁻²⁸，k=7）。"
      "但其引物 SNP <code>rs62036616</code>（p=1.4×10⁻¹²）在 31 kb 外就有一个强出 "
      "<b>140 个数量级</b>的邻居 eQTL：<code>TUFM</code>（rs7187776，p=8.6×10⁻¹⁵³）。"
      "同一窗口内还有 <code>SULT1A1</code>（p=1.6×10⁻²⁷，已知 IBD 易感基因）、<code>SULT1A2</code>、"
      "<code>SPNS1</code>、<code>SGF29</code>——这正是 <b>16p11.2 IBD 易感区域</b>。"
      "ATXN2L 自身的 eQTL 在该区域是<b>最弱的一档</b>。解释上，“ATXN2L 表达导致 UC”远不如"
      "“该 LD 区块内的真实因果基因（如 SULT1A1 一带）同时带动了 ATXN2L 的 eQTL”简约。</p>"
      "<p>相反，<b>LARP1</b> 与 <b>RBPMS2</b>（全血）的 Δlog₁₀P 为<b>负值</b>——它们自己的 eQTL "
      "强于窗口内任何邻居，是少数<b>未</b>被本诊断标记的组合。</p></div>")

    # ---- 5b coloc appendix (explicitly non-valid) ----
    A("<h3>附：探索性的 coloc 计算（<span class='red'>不可作为证据</span>）</h3>")
    if COLOC:
        mx = sorted(COLOC, key=lambda x: -x["PP.H4"])[0]
        A("<div class='warn'><h4>仅作记录，不得引用</h4>"
          "<p>为了说明截断 SNP 集会把结果推向何处，我们把同一批 SNP 代入 coloc.abf 公式"
          "（PriPP：p1=p2=1e-4、p12=1e-5、sd.prior=0.15，含 0.1/0.2/0.3 先验敏感性）。"
          "由于 SNP 集<b>缺少零效应背景</b>，PP 值系统性偏移：最强的一组 %s·%s 也只有 "
          "PP.H4=%.3f，而 MR 最显著的 ATXN2L/CAPRIN1/RBPMS2 全部落在 PP.H4&lt;0.01、"
          "PP.H3≈0.95–0.997。这个“H3 主导”的模式与上一节「邻居 eQTL 碾压目标基因」的发现"
          "方向一致，<b>但它本身不是有效的共定位证据</b>，只能作为旁证。"
          "详见 <code>partC_coloc_results.tsv</code>。</p></div>"
          % (mx["gene"], TISSUE_ZH.get(mx["tissue"], mx["tissue"]), mx["PP.H4"]))
        A(img("fig_partC_coloc.png",
              "附图（仅作记录）　在同一批“仅显著”SNP 上计算的 coloc 后验分解。"
              "因缺少零效应背景，PP.H3 被系统性抬高，<b>不可解读为共定位证据</b>。"))
    else:
        A("<div class='warn'><h4>未生成</h4><p>可匹配的区域内 SNP 数不足。</p></div>")

    # ---- 6 叠合 ----
    A("<h2>六、与 Part A / Part B 的证据叠合</h2>")
    PA_ROBUST = ["VCP", "ZFAND1", "ATXN2", "PUM1", "PUM2", "CSDE1", "TRIM21", "RBPMS"]
    pa_rows = []
    for g in ["VCP", "ZFAND1", "ATXN2", "PUM1", "PUM2", "CSDE1", "TRIM21"]:
        rs = [r for r in SUM if r["gene"] == g]
        if not rs:
            pa_rows.append((g, "无可构建工具变量的组织", "—", "—", "未测"))
            continue
        r = min(rs, key=lambda x: x["IVW_p"])
        fl = nei_by.get("%s|%s" % (r["gene"], r["tissue"]), {}).get("flag") == "True"
        pa_rows.append((g, TISSUE_ZH.get(r["tissue"], r["tissue"]),
                        "%.3f" % r["IVW_b"], sci(r["IVW_p"]),
                        "受到质疑" if fl else ("未通过" if r["IVW_p"] > 0.05 else "名义阳性")))
    A("<div class='card'><h4>Part A 配对双生子队列中最稳健的 7 个 SGRG 在 Part C 中的命运</h4>"
      "<table><tr><th>基因</th><th>最佳可测组织</th><th class='n'>IVW β</th>"
      "<th class='n'>IVW P</th><th>因果判定</th></tr>"
      + "".join("<tr><td>%s</td><td>%s</td><td class='n'>%s</td><td class='n'>%s</td>"
                "<td>%s</td></tr>" % t for t in pa_rows)
      + "</table>"
      "<p>Part A 给出的 7 个候选（VCP↑、ZFAND1↓、ATXN2↓、PUM1/PUM2↓、CSDE1↓、TRIM21↑，"
      "均来自批次正交的同卵双生子配对队列）<b>没有任何一个</b>获得 Part C 的因果支持："
      "要么 p 值远离显著、要么根本没有可用的 cis-eQTL 工具变量。"
      "这是本文最重要的一条“负结果”。</p></div>")

    A("<div class='card'><h4>经 FDR&lt;0.05 的 MR 阳性组合及其可信度</h4>"
      "<table><tr><th>基因</th><th>组织</th><th class='n'>IVW β</th><th class='n'>IVW P</th>"
      "<th class='n'>FDR Q</th><th>LD 混杂诊断</th></tr>"
      + "".join("<tr><td>%s</td><td>%s</td><td class='n'>%.3f</td><td class='n'>%s</td>"
                "<td class='n'>%s</td><td>%s</td></tr>"
                % (r["gene"], TISSUE_ZH.get(r["tissue"], r["tissue"]), r["IVW_b"],
                   sci(r["IVW_p"]), sci(r["IVW_q"]),
                   "<b class='red'>高度疑似 LD 混杂</b>"
                   if nei_by.get("%s|%s" % (r["gene"], r["tissue"]), {}).get("flag") == "True"
                   else ("窗口内存在更强的邻居 eQTL，但仍待确证"
                         if float(nei_by.get("%s|%s" % (r["gene"], r["tissue"]),
                                             {"delta_log10p": 0})["delta_log10p"]) > 0
                         else "目标基因 eQTL 为窗口内最强，未见混杂信号"))
                for r in sorted(sigfdr, key=lambda x: x["IVW_p"]))
      + "</table><h4>三层证据同时支持的「最终核心 SGRG」</h4>"
      "<p><b>没有。</b>没有任何一个 SGRG 同时具备（a）在批次正交队列中方向可复现的表达改变、"
      "（b）通过 FDR 且经 LD 混杂诊断后仍然成立的因果估计。"
      "按整合方案 Fig 10 的叠合逻辑，本文<b>不能收敛到任何“因果核心基因”</b>——"
      "这本身是一个可以写进文章的核心发现："
      "<b>SGRG 与 UC 的关系停留在关联层面，且该关联在很大程度上可由批次效应与 LD 结构解释。</b></p></div>")

    # ---- 7 局限 ----
    A("<h2>七、局限与必须声明的技术问题</h2>")
    A("<div class='warn'><h4>必须写进 Limitations 的八项</h4><ol>"
      "<li><b>无法做有效共定位（最关键）</b>：本环境可达的 GTEx 接口只返回<b>已显著</b>的变异-基因对"
      "（VCP 窗口 76 个 SNP 的 p 全部 &lt;2×10⁻⁴），缺少零效应背景，coloc.abf 的后验概率不可辨识；"
      "完整 all-pairs 数据（每组织 ~4 GB）在本环境不可得。因此本文只能做「邻居 eQTL 优势度」这一"
      "可用同一批数据严格计算的替代诊断，<b>不能声称完成了共定位</b>。</li>"
      "<li><b>LD 剪枝为距离近似</b>：本环境无法访问 LD 参考面板（Ensembl LD 接口返回 000、"
      "无 plink/tabix），独立性以 10 kb 距离剪枝 + 完美代理（NES/SE 完全相同）去重近似 r²&lt;0.001。"
      "若工具变量间存在残余 LD，IVW 的 SE 会被低估。全流程共剔除 953 个完美代理 SNP。</li>"
      "<li><b>GTEx 样本量小</b>：各组织 eQTL 样本量仅 205–670，工具变量效应估计本身有噪声；"
      "其后果是回归稀释（bias toward the null），即阳性结果偏保守。</li>"
      "<li><b>SE 为反推而非直读</b>：GTEx API 只给 NES 与 P 值，SE 由 <code>|NES|/|z|</code> 反推，"
      "前提是 GTEx 的 P 值来自对该斜率的 t 检验——这是 GTEx 的实现方式，但属于间接推导。</li>"
      "<li><b>结局为二分类、暴露为定量</b>：因果估计处于「每 SD 表达量 → UC log-OR」的尺度，"
      "不可直接与临床 OR 比较。</li>"
      "<li><b>共定位先验敏感性</b>：coloc 的 sd.prior 会改变 PP 分配，本文已给出 0.1/0.2/0.3 三档敏感性。</li>"
      "<li><b>顺式限制</b>：仅用 cis-eQTL，虽大幅降低水平多效性风险，但 cis 区域本身可能同时影响"
      "邻近基因（如 C9orf72、HECA 等基因密度高的位点）。</li>"
      "<li><b>本环境网络受限</b>：eQTLGen（n=31,684，血样）主路径与全部 GWAS 汇总统计的完整下载"
      "不可行；本报告暴露端全部来自 GTEx v8 API，结局端为一份 %d MB 的 GWAS 文件（并行分块下载获得）。"
      "eQTLGen 的缺失意味着「全血无 GTEx 信号」的基因并不等于无工具变量。</li>"
      "</ol></div>"
      % round(os.path.getsize(os.path.join(ROOT, "data", "raw", "gwas",
                                           "GCST004133_uc_harmonised.h.tsv.gz")) / 1e6)
      if os.path.exists(os.path.join(ROOT, "data", "raw", "gwas",
                                     "GCST004133_uc_harmonised.h.tsv.gz")) else "")

    A("<footer>结果文件：<code>results/partC_mr_results.tsv</code>、"
      "<code>results/partC_mr_geneset.tsv</code>、<code>results/partC_coloc_results.tsv</code>、"
      "<code>results/partC_mr_summary.json</code>；"
      "脚本：<code>scripts/partC_get_eqtl.py</code>、<code>partC_get_eqtl_ext.py</code>、"
      "<code>partC_a_instruments.py</code>、<code>partC_extract_gwas.py</code>、"
      "<code>partC_c_mr.py</code>、<code>mr_core.py</code>、<code>partC_d_report.py</code>；"
      "下载脚本：<code>scripts/download_gwas_parallel.py</code>。<br>"
      "暴露数据：GTEx v8（gtexportal.org API v2）；结局数据：GWAS Catalog "
      "GCST004133 harmonised（ftp.ebi.ac.uk）。</footer>")
    A("</div></body></html>")

    out = os.path.join(RES, "PartC_SGRG_UC_孟德尔随机化因果报告.html")
    open(out, "w", encoding="utf-8").write("\n".join(H))
    print("WROTE", out, os.path.getsize(out), "bytes")


if __name__ == "__main__":
    main()
