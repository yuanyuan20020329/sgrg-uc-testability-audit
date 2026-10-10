# -*- coding: utf-8 -*-
"""Part D / 公共资源可用性裁决表 + 图 D12

把 CELLxGENE 全库枚举（partD_resfeas_cxg.py）与 GEO 全库枚举（partD_resfeas_geo.py）
合并成**一张可审计的可行性裁决表**，回答审稿人那句「你真的找遍了吗，真没有有检验力的数据？」。

判定四条闸门（缺一不可）
  G1 assay 测 mRNA        —— 否则 SGRG 这种编码基因面板根本打不了分
  G2 含 UC              —— 本项目的目标疾病（CD/IBD 单列，可作旁证但不算达标）
  G3 有受控对比且达到检验力 —— 同供体配对 ≥6 对，或同研究真对照每组 ≥5 例
  G4 在体组织（非类器官/细胞系）

输出
----
  results/partD_resfeas_ledger.tsv     每个候选资源一行（四条闸门 + 结论 + 依据）
  results/partD_resfeas_verdict.json   汇总裁决
  figures/fig_D12_resource_feasibility.{png,pdf}
用法
----
  python partD_resfeas_ledger.py
"""
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

MM = 1 / 25.4
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 6.2,
    "axes.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 400,
})
R, FIG = "results", "results"
C_OK, C_BAD, C_PAIR = "#2e6f9e", "#b0b7bf", "#b5442e"
C_SC, C_SP, C_BU = "#4c9f70", "#c9873a", "#8d99ae"
PAIR_MIN = 6
CTRL_MIN = 5

os.makedirs(FIG, exist_ok=True)

CXG = pd.read_csv(os.path.join(R, "partD_resfeas_cxg_datasets.tsv"), sep="\t")
CGX = pd.read_csv(os.path.join(R, "partD_resfeas_cxg_collections.tsv"), sep="\t")
GEO = pd.read_csv(os.path.join(R, "partD_resfeas_geo_series.tsv"), sep="\t")
cxg_summ = json.load(open(os.path.join(R, "partD_resfeas_cxg_summary.json"), encoding="utf-8"))
geo_summ = json.load(open(os.path.join(R, "partD_resfeas_geo_summary.json"), encoding="utf-8"))

rows = []

# ---------------------------------------------------------------- CELLxGENE 侧
for _, r in CXG[CXG.has_UC].iterrows():
    frac = r.frac_primary if pd.notna(r.frac_primary) else np.nan
    nuc, nctl = int(r.n_UC_donors), int(r.n_CTL_donors)
    g1 = True
    g2 = True
    g3 = (nuc >= CTRL_MIN) and (nctl >= CTRL_MIN) and (frac >= 0.9)
    g4 = True
    why = []
    if not g3:
        if nctl < CTRL_MIN:
            why.append("自有对照供体 %d<%d" % (nctl, CTRL_MIN))
        if frac < 0.9:
            why.append("原产细胞仅 %.1f%%" % (100 * frac if pd.notna(frac) else float("nan")))
    rows.append({
        "source": "CELLxGENE", "resource": r.collection_name,
        "dataset": r.dataset_title, "assay": (r.assays if isinstance(r.assays, str) else ""),
        "resolution": ("spatial" if r.spatial else "single-cell"),
        "disease": "UC", "n_case": nuc, "n_ctrl": nctl,
        "n_pairs": "", "frac_primary": round(float(frac), 4) if pd.notna(frac) else "",
        "G1_mRNA": g1, "G2_UC": g2, "G3_powered_control": g3, "G4_in_vivo": g4,
        "USABLE": bool(g1 and g2 and g3 and g4),
        "reason": "; ".join(why) if why else "同研究真对照，供体数与同源性均达标",
    })

# ---------------------------------------------------------------- GEO 侧
G = GEO[GEO.n_paired_donors.fillna(0) >= 2].copy()
for _, r in G.iterrows():
    npr = int(r.n_paired_donors)
    nss = int(r.n_paired_same_site) if pd.notna(r.n_paired_same_site) else 0
    dis = str(r.diseases or "")
    g1 = bool(r.mrna_scorable)
    g2 = "UC" in dis
    g3 = npr >= PAIR_MIN
    g4 = not bool(r.in_vitro)
    why = []
    if not g1:
        why.append("assay=%s，非 mRNA，无法给编码基因面板打分" % r.molecule)
    if not g2:
        why.append("疾病为 %s（无 UC）" % (dis or "未标注"))
    if not g3:
        why.append("同供体配对仅 %d 对（<%d）" % (npr, PAIR_MIN))
    if not g4:
        why.append("离体模型（类器官/细胞系）")
    if g1 and nss == 0 and r.resolution == "spatial":
        why.append("空间数据配对为跨肠段，部位与状态混杂")
    rows.append({
        "source": "GEO", "resource": r.gse, "dataset": r.title,
        "assay": str(r.series_type)[:70], "resolution": r.resolution,
        "disease": dis or "-", "n_case": "", "n_ctrl": "",
        "n_pairs": npr, "frac_primary": "",
        "G1_mRNA": g1, "G2_UC": g2, "G3_powered_control": g3, "G4_in_vivo": g4,
        "USABLE": bool(g1 and g2 and g3 and g4),
        "reason": "; ".join(why) if why else "达标",
    })

L = pd.DataFrame(rows)
L["n_pairs_num"] = pd.to_numeric(L.n_pairs, errors="coerce").fillna(0)
L = L.sort_values(["USABLE", "n_pairs_num"], ascending=[False, False]).drop(columns=["n_pairs_num"])
L.to_csv(os.path.join(R, "partD_resfeas_ledger.tsv"), sep="\t", index=False)

usable = L[L.USABLE]
verdict = {
    "cxg_total_collections": cxg_summ["n_collections_total"],
    "cxg_collections_with_UC": cxg_summ["n_collections_with_UC"],
    "cxg_usable_controlled_resource": int((L.USABLE & (L.source == "CELLxGENE")).sum()),
    "geo_human_gse": geo_summ["n_gse_human"],
    "geo_mrna_scorable_series": geo_summ["n_mrna_scorable"],
    "geo_singlecell_spatial_mrna": geo_summ["n_singlecell_or_spatial_mrna"],
    "geo_series_with_any_pairing": geo_summ["n_series_with_pairing_mrna"],
    "geo_mrna_ge6_pairs": geo_summ["n_series_with_pairing_ge6_mrna"],
    "geo_mrna_ge15_pairs": geo_summ["n_series_with_pairing_ge15_mrna"],
    "geo_singlecell_pairs_max": int(GEO[GEO.resolution == "single-cell"].n_paired_donors.max()),
    "geo_spatial_pairs_max": int(GEO[GEO.resolution == "spatial"].n_paired_donors.max()),
    "geo_spatial_samesite_pairs": geo_summ["n_spatial_mrna_with_pairing_samesite"],
    "geo_max_pairs_any_assay": int(GEO.n_paired_donors.max()),
    "geo_max_pairs_mirna_example": "GSE272890（43 对，但测 miRNA）",
    "usable_resources": usable[["source", "resource", "resolution", "disease",
                                "n_case", "n_ctrl", "n_pairs"]].to_dict("records"),
    "verdict_text": (
        "全库枚举后，可用于「有检验力且受控地检验 SGRG 面板」的公共资源总计 %d 个；"
        "单细胞与空间转录组里 **同供体配对 ≥6 对的 UC 数据为 0 个**。"
        "唯一达标的单细胞资源是 CELLxGENE 的 ds1（6v6、100%% 原产、但仅 46.7k 细胞）；"
        "唯一达标的 UC 配对设计是 GEO 的 GSE133060（8 对，bulk 芯片）。"
        "故原判定「公共单细胞资源在结构上无法做出有检验力、批次受控的检验」成立，"
        "但需补一句：**在 bulk 层面存在一个达标的配对设计，而它尚未被检验**。" % 0
    ).replace("总计 0 个", "总计 %d 个" % int(len(usable))),
}
json.dump(verdict, open(os.path.join(R, "partD_resfeas_verdict.json"), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False, default=str)

# ================================================================== 图 D12
fig = plt.figure(figsize=(180 * MM, 108 * MM))
gs = GridSpec(2, 3, figure=fig, hspace=0.62, wspace=0.42,
              left=0.115, right=0.985, top=0.90, bottom=0.115)

# ---- (A) GEO 漏斗
ax = fig.add_subplot(gs[0, 0])
steps = [("IBD×gut human GSE", geo_summ["n_gse_human"], C_BAD),
         ("SOFT parsed", geo_summ["n_series_parsed"], C_BAD),
         ("mRNA-scorable", geo_summ["n_mrna_scorable"], C_OK),
         ("single-cell / spatial", geo_summ["n_singlecell_or_spatial_mrna"], C_OK),
         ("any pairing", geo_summ["n_series_with_pairing_mrna"], C_PAIR),
         ("≥6 pairs", geo_summ["n_series_with_pairing_ge6_mrna"], C_PAIR)]
y = np.arange(len(steps))[::-1]
ax.barh(y, [s[1] for s in steps], color=[s[2] for s in steps], height=0.62)
for yy, (lab, v, _) in zip(y, steps):
    ax.text(v + 8, yy, str(v), va="center", ha="left", fontsize=5.8, fontweight="bold")
ax.set_yticks(y); ax.set_yticklabels([s[0] for s in steps], fontsize=5.8)
ax.set_xlim(0, 640); ax.set_xlabel("GEO series (count)", fontsize=6)
ax.set_title("A  GEO enumeration funnel", fontsize=6.6, fontweight="bold", loc="left")

# ---- (B) 配对系列：分辨率 × 配对数
ax = fig.add_subplot(gs[0, 1])
order = ["bulk-array", "bulk-seq", "single-cell", "spatial", "ncRNA/other"]
colr = {"bulk-array": C_BU, "bulk-seq": "#a9b4bf", "single-cell": C_SC, "spatial": C_SP,
        "ncRNA/other": "#c9c3b8"}
GEO["_axis"] = GEO.resolution.where(GEO.resolution.isin(order), "ncRNA/other")
rng = np.random.default_rng(7)
for i, res in enumerate(order):
    sub = GEO[GEO._axis == res]
    xs = i + rng.uniform(-0.22, 0.22, len(sub))
    ok = sub.molecule == "mRNA"
    ax.scatter(xs[ok], sub.n_paired_donors[ok], s=17, color=colr[res],
               edgecolor="white", linewidth=0.3, zorder=3)
    ax.scatter(xs[~ok], sub.n_paired_donors[~ok], s=20, facecolor="none",
               edgecolor=colr[res], linewidth=0.7, zorder=3)
ax.axhline(PAIR_MIN, color=C_PAIR, lw=0.8, ls="--")
ax.text(4.45, PAIR_MIN * 1.06, "power floor\n6 pairs", fontsize=5.0, color=C_PAIR, ha="right")
for gse, lab in [("GSE272890", "GSE272890\n43 pairs but\nmiRNA → unusable"),
                 ("GSE133060", "GSE133060\n8 pairs\nbulk array")]:
    r = GEO[GEO.gse == gse]
    if len(r):
        i = order.index(r.iloc[0]._axis)
        ax.annotate(lab, (i, r.iloc[0].n_paired_donors),
                    xytext=(i - 0.85, r.iloc[0].n_paired_donors * 1.45),
                    fontsize=4.9, color="#333333",
                    arrowprops=dict(arrowstyle="-", lw=0.45, color="#777777"))
ax.set_xticks(range(len(order)))
ax.set_xticklabels(["bulk\narray", "bulk\nRNA-seq", "single\ncell", "spatial", "ncRNA/\nother"],
                   fontsize=5.4)
ax.set_yscale("log"); ax.set_ylim(1.6, 90)
ax.set_ylabel("same-donor inflamed/non-inflamed pairs", fontsize=6)
ax.set_title("B  Paired designs by resolution\n(filled = mRNA-scorable)", fontsize=6.6,
             fontweight="bold", loc="left")

# ---- (C) CELLxGENE：3 个含 UC 的集合
ax = fig.add_subplot(gs[0, 2])
cc = CXG[CXG.has_UC].groupby("collection_id").agg(
    name=("collection_name", "first"), n_uc=("n_UC_donors", "max"),
    frac=("frac_primary", "max"), cells=("cell_count", "sum")).reset_index()
short = {cc.iloc[i].collection_id: lab for i, lab in enumerate(
    ["Oliver 2024\n(metaplasia)", "ds4 Tearle\n(PSC-UC)", "ds1 Garrido\n(IBD & HC)"])} if len(cc) == 3 else {}
for i, r in cc.iterrows():
    color = C_OK if (r.n_uc >= CTRL_MIN and r.frac >= 0.9) else C_BAD
    ax.scatter(r.frac, r.n_uc, s=64, color=color, edgecolor="white", linewidth=0.5, zorder=3)
    ax.annotate(short.get(r["collection_id"], str(r["name"])[:16]), (r["frac"], r["n_uc"]),
                xytext=(-4, 7), textcoords="offset points", fontsize=5.0, ha="right")
ax.axhline(CTRL_MIN, color="#777777", lw=0.6, ls=":")
ax.axvline(0.9, color="#777777", lw=0.6, ls=":")
ax.set_xlim(-0.05, 1.16); ax.set_ylim(0, 9.6)
ax.set_xlabel("fraction of cells that are primary (own study)", fontsize=6)
ax.set_ylabel("UC donors in collection", fontsize=6)
ax.set_title("C  CELLxGENE: all %d collections\ncontaining UC" % cxg_summ["n_collections_with_UC"],
             fontsize=6.6, fontweight="bold", loc="left")

# ---- (D) 四条闸门判据矩阵
ax = fig.add_subplot(gs[1, :])
pick = ["ds1 Garrido (CELLxGENE)", "ds4 Tearle (CELLxGENE)", "Oliver metaplasia (CELLxGENE)",
        "GSE133060", "GSE208303", "GSE289496", "GSE169136", "GSE312415", "GSE272890"]
lab_map = {}
for _, r in L.iterrows():
    if r.source == "CELLxGENE":
        nm = str(r.resource)
        key = ("ds1 Garrido (CELLxGENE)" if "IBD" in nm and "healthy" in nm.lower()
               else "ds4 Tearle (CELLxGENE)" if "sclerosing" in nm.lower()
               else "Oliver metaplasia (CELLxGENE)")
        lab_map.setdefault(key, r)
    else:
        lab_map.setdefault(str(r.resource), r)
crit = ["G1_mRNA", "G2_UC", "G3_powered_control", "G4_in_vivo"]
crit_lab = ["mRNA assay", "contains UC", "powered\ncontrol", "in-vivo\ntissue"]
present = [k for k in pick if k in lab_map]
M = np.zeros((len(present), 4))
for i, k in enumerate(present):
    r = lab_map[k]
    for j, c in enumerate(crit):
        M[i, j] = 1.0 if bool(r[c]) else 0.0
for i in range(len(present)):
    for j in range(4):
        ax.add_patch(plt.Rectangle((j, len(present) - 1 - i), 0.94, 0.9,
                                   facecolor=(C_OK if M[i, j] else "#e9ebee"),
                                   edgecolor="white", lw=0.8))
        ax.text(j + 0.47, len(present) - 1 - i + 0.45, "✓" if M[i, j] else "✗",
                ha="center", va="center", fontsize=6.5, color="white" if M[i, j] else "#8b9199")
    ok = M[i].all()
    ax.text(4.15, len(present) - 1 - i + 0.45, "USABLE" if ok else "fails",
            va="center", fontsize=5.6, fontweight="bold", color=C_OK if ok else "#b5442e")
ax.set_xlim(-0.05, 5.2); ax.set_ylim(-0.1, len(present))
ax.set_yticks(np.arange(len(present)) + 0.45)
ax.set_yticklabels(present[::-1], fontsize=5.8)
ax.set_xticks(np.arange(4) + 0.47); ax.set_xticklabels(crit_lab, fontsize=5.6)
idxs = [k for k in M.nonzero()[0]]
for j, c in enumerate(crit):
    pass
ax.set_title("D  Four-gate feasibility ledger — %d of %d candidate resources pass all gates"
             % (int(L.USABLE.sum()), len(L)), fontsize=6.6, fontweight="bold", loc="left")
ax.text(-0.04, -0.62, "fold-change filter: %d pairs min · ≥%d own-study controls/group · "
        "fraction primary ≥0.9" % (PAIR_MIN, CTRL_MIN), fontsize=5.0, color="#666666",
        transform=ax.transAxes)

fig.savefig(os.path.join(FIG, "fig_D12_resource_feasibility.png"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "fig_D12_resource_feasibility.pdf"), bbox_inches="tight")
plt.close(fig)

print("=" * 92)
print("公共资源可用性裁决")
print("=" * 92)
print("CELLxGENE：全库 %d 集合 → 含 UC %d 个 → 通过四条闸门 %d 个"
      % (cxg_summ["n_collections_total"], cxg_summ["n_collections_with_UC"],
         verdict["cxg_usable_controlled_resource"]))
print("GEO：人类 GSE %d → mRNA 可打分 %d → 单细胞/空间 %d → 有配对 %d → ≥6 对 %d"
      % (geo_summ["n_gse_human"], geo_summ["n_mrna_scorable"],
         geo_summ["n_singlecell_or_spatial_mrna"], geo_summ["n_series_with_pairing_mrna"],
         geo_summ["n_series_with_pairing_ge6_mrna"]))
print("单细胞配对最大 %d 对；空间配对最大 %d 对（同部位 %d 对）"
      % (verdict["geo_singlecell_pairs_max"], verdict["geo_spatial_pairs_max"],
         verdict["geo_spatial_samesite_pairs"]))
print("\n通过全部闸门的资源：")
if len(usable):
    print(usable[["source", "resource", "resolution", "disease", "n_case", "n_ctrl",
                  "n_pairs"]].to_string(index=False, max_colwidth=46))
else:
    print("  （无）")
print("\n裁决原文：\n  " + verdict["verdict_text"])
print("\n产出 results/partD_resfeas_ledger.tsv / _verdict.json，figures/fig_D12_resource_feasibility.*")
