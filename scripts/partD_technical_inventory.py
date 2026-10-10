# -*- coding: utf-8 -*-
"""步 1｜单细胞数据集级「技术轴清单」——A2 队列级宽表（partA2_cohorts.tsv）的等构对应物

目的
----
把 A2 在 bulk 队列上做的「技术轴占据度」度量搬到单细胞分辨率，
用**同一张表的同一批列**描述 4 套 UC 黏膜 scRNA-seq：

  技术轴 A（深度）：depth_fold / auc_depth_to_disease / r2_depth_disease / rho_score_nUMI
  技术轴 B（批次）：从 partD_batch_audit_summary.tsv 取「对照 donor 是否与病例落在同一批次」
  主成分占据度    ：pc_tech_share = PC1–3 中可被深度解释的方差占比（A2 的 auc_pc1_batch 对应物）
  下游后果        ：名义显著细胞类型数 / 干预后显著数

★ 诚实纪律
- 只读已有结果表，不重跑 2.7 GB h5ad；
- 名义显著数一律标清检验口径（置换 / MWU）与 FDR 是否校正；
- ds2/ds3 没有做抽稀，`n_ct_after_control` 留空而不是填 0。

用法
----
  python partD_technical_inventory.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from partD_common import bh  # noqa: E402

OUT = "results"
MIN_DONORS_RHO = 5      # 细胞类型内至少要几个 donor 才算 ρ / R²
COV_MIN = 0.60          # 细胞类型要覆盖多少比例的 donor 才进 PCA
N_PC = 3

DATASETS = [
    dict(tag="ds1", label="ds1 Garrido-Trigo 2024（IBD 全黏膜）",
         pb="partD_ds1_pseudobulk_score.tsv",
         tests="partD_ds1_pseudobulk_tests.tsv"),
    dict(tag="ds2", label="ds2 Oliver2024 髓系（结肠）",
         pb="partD_ds2_myeloid_pseudobulk_score.tsv",
         tests="partD_ds2_myeloid_uc_vs_hc_tests.tsv"),
    dict(tag="ds3", label="ds3 Oliver2024 大肠（上皮为主）",
         pb="partD_ds3_largeintestine_pseudobulk_score.tsv",
         tests="partD_ds3_largeintestine_uc_vs_hc_tests.tsv"),
    dict(tag="ds4", label="ds4 Tearle2026（UC/PSC 黏膜）",
         pb="partD_ds4_tearle_pseudobulk_score.tsv",
         tests="partD_ds4_tearle_uc_vs_hc_tests.tsv",
         tests_alt="partD_ds4_tearle_ds1k_uc_vs_hc_tests.tsv"),
]

log = lambda *x: print(*x, flush=True)


# ------------------------------------------------------------------ 工具
def auc_rank(pos, neg):
    """AUC = P(pos > neg)，用 Mann-Whitney U 精确算。>0.5 表示 UC 侧数值更大。"""
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    pos = pos[np.isfinite(pos)]
    neg = neg[np.isfinite(neg)]
    if len(pos) < 2 or len(neg) < 2:
        return np.nan
    U = mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return float(U / (len(pos) * len(neg)))


def r2_cont(x, y):
    """一元回归 R²（连续解释变量）。"""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan
    return float(np.corrcoef(x[m], y[m])[0, 1] ** 2)


def r2_binary(y, grp):
    """二分类变量的点二列 R²。"""
    y = np.asarray(y, float)
    g = np.asarray(grp, float)
    return r2_cont(g, y)


def donor_table(df):
    """donor 级表。

    深度定义 = 该 donor 各细胞类型「中位 nUMI」的中位数（每个细胞类型等权，先到 donor 再到组），
    即**分析单元 = donor** 的口径。若改用细胞级中位（大 donor 主导）会得到显著不同的比值，
    两种口径都写进表里（见 depth_fold_celllevel 与 meta 说明）。
    """
    rows = []
    for d, sub in df.groupby("donor"):
        depth = float(np.median(sub["nUMI"].astype(float).values))
        w = np.maximum(sub["n_cells"].astype(float).values, 1)
        depth_cw = float(np.average(sub["nUMI"].astype(float).values, weights=w))
        dis = sub["disease"].iloc[0]
        rows.append(dict(donor=d, disease=dis, depth=depth, depth_cw=depth_cw,
                         log_depth=np.log10(max(depth, 1)),
                         log_depth_cw=np.log10(max(depth_cw, 1)),
                         n_cells=int(sub["n_cells"].astype(float).sum()),
                         n_ct=int((sub["n_cells"] >= 20).sum())))
    return pd.DataFrame(rows)


def grp_of(disease):
    d = str(disease).lower()
    if "ulcer" in d:
        return "UC"
    if d in ("normal", "healthy", "control", "healthy control"):
        return "CTRL"
    return "OTHER"


# ------------------------------------------------------------------ 主计算
def inventory_one(ds):
    pb = pd.read_csv(os.path.join(OUT, ds["pb"]), sep="\t")
    dt = donor_table(pb)
    dt["grp"] = dt["disease"].map(grp_of)
    key = dt[dt.grp.isin(["UC", "CTRL"])].copy()

    uc = key[key.grp == "UC"]
    ct = key[key.grp == "CTRL"]
    rec = dict(dataset=ds["tag"], label=ds["label"],
               n_donor_total=int(dt.shape[0]), n_uc_donor=int(uc.shape[0]),
               n_ctrl_donor=int(ct.shape[0]),
               n_donor_excluded_other=int((dt.grp == "OTHER").sum()),
               n_cell_pseudobulk=int(dt.n_cells.sum()),
               n_celltype_pb=int(pb.cell_type.nunique()))

    # ---- 深度轴
    rec["med_nUMI_uc"] = float(uc.depth.median()) if len(uc) else np.nan
    rec["med_nUMI_ctrl"] = float(ct.depth.median()) if len(ct) else np.nan
    rec["depth_fold_uc_over_ctrl"] = (rec["med_nUMI_uc"] / rec["med_nUMI_ctrl"]
                                      if len(uc) and len(ct) and rec["med_nUMI_ctrl"] > 0 else np.nan)
    a = auc_rank(uc.log_depth, ct.log_depth)
    rec["auc_depth_to_disease"] = a
    rec["depth_separability"] = max(a, 1 - a) if np.isfinite(a) else np.nan
    # ---- 细胞类型内：ρ(评分, 深度) 与 R²(深度 ~ 疾病)
    rho, r2d, ct_used = [], [], 0
    for c, sub in pb.groupby("cell_type"):
        sub = sub[sub.disease.isin(dt[dt.grp.isin(["UC", "CTRL"])].disease.unique())]
        if sub.donor.nunique() < MIN_DONORS_RHO:
            continue
        ct_used += 1
        r = spearmanr(sub["nUMI"].astype(float), sub["SGRG"].astype(float)).statistic
        if np.isfinite(r):
            rho.append(abs(float(r)))
        s = sub[sub.disease.map(grp_of).isin(["UC", "CTRL"])].copy()
        if s.disease.map(grp_of).nunique() == 2 and s.disease.map(grp_of).value_counts().min() >= 2:
            rr = r2_binary(np.log10(s["nUMI"].astype(float)), (s.disease.map(grp_of) == "UC").astype(float))
            if np.isfinite(rr):
                r2d.append(float(rr))
    rec["n_ct_for_rho"] = ct_used
    rec["rho_score_nUMI_med"] = float(np.median(rho)) if rho else np.nan
    rec["rho_score_nUMI_max"] = float(np.max(rho)) if rho else np.nan
    rec["n_ct_for_r2"] = len(r2d)
    rec["r2_depth_disease_med"] = float(np.median(r2d)) if r2d else np.nan
    rec["r2_depth_disease_max"] = float(np.max(r2d)) if r2d else np.nan

    # ---- 主轴占据度（A2 的 auc_pc1_batch 对应物）
    M = pb.pivot_table(index="donor", columns="cell_type", values="SGRG")
    keep = M.columns[M.notna().mean() >= COV_MIN]
    M = M[keep]
    M = M.loc[[d for d in M.index if d in set(key.donor)]]
    if M.shape[1] >= 2 and M.shape[0] >= 5:
        Mf = M.apply(lambda s: s.fillna(s.mean()), axis=0)
        Z = (Mf - Mf.mean()) / Mf.std(ddof=0).replace(0, np.nan)
        Z = Z.loc[:, Z.notna().all(axis=0)].fillna(0.0)
        # 去掉零方差列
        Z = Z.loc[:, Z.std(ddof=0) > 1e-9]
        k = min(N_PC, min(Z.shape) - 1)
        U, S, Vt = np.linalg.svd(Z.values - Z.values.mean(0), full_matrices=False)
        ve = (S ** 2) / np.sum(S ** 2)
        pc = U[:, :k] * S[:k]
        dmap = dt.set_index("donor")
        dep = dmap.loc[Z.index, "log_depth"].values
        dep_cw = dmap.loc[Z.index, "log_depth_cw"].values
        g = (dmap.loc[Z.index, "grp"] == "UC").astype(float).values
        r2_dep = [r2_cont(pc[:, i], dep) for i in range(k)]
        r2_dep_cw = [r2_cont(pc[:, i], dep_cw) for i in range(k)]
        r2_dis = [r2_binary(pc[:, i], g) for i in range(k)]
        rec["pc_var_expl_1"] = float(ve[0])
        rec["pc_tech_share"] = float(np.nansum([ve[i] * (r2_dep[i] if np.isfinite(r2_dep[i]) else 0)
                                                for i in range(k)]))
        rec["pc_tech_share_cw"] = float(np.nansum([ve[i] * (r2_dep_cw[i] if np.isfinite(r2_dep_cw[i]) else 0)
                                                   for i in range(k)]))
        rec["pc_disease_share"] = float(np.nansum([ve[i] * (r2_dis[i] if np.isfinite(r2_dis[i]) else 0)
                                                   for i in range(k)]))
        rec["pc1_r2_depth"] = float(r2_dep[0]) if np.isfinite(r2_dep[0]) else np.nan
        rec["pc1_r2_depth_cw"] = float(r2_dep_cw[0]) if np.isfinite(r2_dep_cw[0]) else np.nan
        rec["pc1_r2_disease"] = float(r2_dis[0]) if np.isfinite(r2_dis[0]) else np.nan
        rec["main_axis_gap"] = (rec["pc1_r2_depth"] - rec["pc1_r2_disease"]
                                if np.isfinite(rec["pc1_r2_depth"]) and np.isfinite(rec["pc1_r2_disease"])
                                else np.nan)
        rec["n_pc_used"] = k
        rec["n_ct_in_pca"] = int(Z.shape[1])
        rec["n_donor_in_pca"] = int(Z.shape[0])
    else:
        for kk in ("pc_var_expl_1", "pc_tech_share", "pc_tech_share_cw", "pc_disease_share",
                   "pc1_r2_depth", "pc1_r2_depth_cw", "pc1_r2_disease", "main_axis_gap"):
            rec[kk] = np.nan
        rec["n_pc_used"] = 0
        rec["n_ct_in_pca"] = int(M.shape[1])
        rec["n_donor_in_pca"] = int(M.shape[0])

    # ---- 下游后果
    tp = os.path.join(OUT, ds["tests"])
    n_p05 = n_q05 = minp = np.nan
    n_tested_denom = np.nan
    if os.path.exists(tp):
        t = pd.read_csv(tp, sep="\t")
        if "score" in t.columns:
            t = t[t.score == "SGRG"]
        if "contrast" in t.columns:            # ds1：多对比，取 UC_vs_normal
            t = t[t.contrast == "UC_vs_normal"]
        if len(t):
            n_tested_denom = int(t["cell_type"].nunique())   # ★ 分母必须与检验表一致
            pv = t["p_perm"].astype(float)
            n_p05 = int((pv < 0.05).sum())
            minp = float(pv.min())
            qcol = "q_BH_within_score" if "q_BH_within_score" in t.columns else ("q_BH" if "q_BH" in t.columns else None)
            if qcol:
                n_q05 = int((t[qcol].astype(float) < 0.05).sum())
    rec["n_ct_p05_perm"] = n_p05
    rec["n_ct_q05"] = n_q05
    rec["min_p_perm"] = minp
    rec["n_ct_in_test_table"] = n_tested_denom
    rec["n_ct_tested"] = n_tested_denom

    # 抽稀（ds4 only）
    alt = ds.get("tests_alt")
    rec["n_ct_p05_after_rarefy"] = np.nan
    rec["n_ct_q05_after_rarefy"] = np.nan
    rec["min_p_after_rarefy"] = np.nan
    if alt and os.path.exists(os.path.join(OUT, alt)):
        t = pd.read_csv(os.path.join(OUT, alt), sep="\t")
        t = t[t.score == "SGRG"] if "score" in t.columns else t
        pv = t["p_perm"].astype(float)
        rec["n_ct_p05_after_rarefy"] = int((pv < 0.05).sum())
        qcol = "q_BH_within_score" if "q_BH_within_score" in t.columns else ("q_BH" if "q_BH" in t.columns else None)
        rec["n_ct_q05_after_rarefy"] = int((t[qcol].astype(float) < 0.05).sum()) if qcol else np.nan
        rec["min_p_after_rarefy"] = float(pv.min())
    return rec


# ------------------------------------------------------------------ 批次轴（读审计表）
AUDIT_KEY = {"ds1": "ds1_garrido", "ds2": "ds2_myeloid",
             "ds3": "ds3_largeintestine", "ds4": "ds4_tearle"}


def batch_axis():
    p = os.path.join(OUT, "partD_batch_audit_summary.tsv")
    if not os.path.exists(p):
        return {}
    b = pd.read_csv(p, sep="\t")
    out = {}
    for _, r in b.iterrows():
        vars_, share = [], []
        for c in b.columns:
            if c.startswith("n_batch_mixed_"):
                v = c[len("n_batch_mixed_"):]
                if pd.notna(r.get(c)):
                    vars_.append(v)
                    share.append(float(r.get("frac_ctrl_donors_samebatch_" + v, np.nan)))
        if not vars_:
            out[r["dataset"]] = dict(vars="", best_var="", best_share=np.nan,
                                     worst_var="", worst_share=np.nan, nested_flag=np.nan)
            continue
        s = pd.Series(share, index=vars_)
        out[r["dataset"]] = dict(
            vars=",".join(vars_),
            best_var=s.idxmax(), best_share=float(s.max()),
            worst_var=s.idxmin(), worst_share=float(s.min()),
            nested_flag=int(float(s.min()) < 0.5))
    return out


# ---- ds4 来源审计（partD_ds4_prov_summary.json + 逐细胞表实测）中的三个替代对照口径
PROV = dict(
    # 细胞级中位深度：UC 5674 / 全部对照 1305（4.35×）/ 仅本研究自有对照 5440（1.04×）
    frac_ctrl_cells_primary=10519 / 188501,
    depth_fold_celllevel_naive=5674.0 / 1305.0,
    depth_fold_celllevel_primaryctrl=5674.0 / 5440.0,
    n_ct_p05_after_primary_ctrl=1,      # 口径 B_primary：23 个细胞类型，MWU 与置换各 1 个
    n_ct_in_primary_ctrl=23,
    n_ct_p05_after_within_donor=0,      # 口径 C_inflamed：15 个细胞类型，同 donor 炎症配对
    n_ct_in_within_donor=15,
)


recs = []
for ds in DATASETS:
    log("=" * 70)
    log("计算", ds["tag"])
    r = inventory_one(ds)
    for k in ("med_nUMI_uc", "med_nUMI_ctrl", "depth_fold_uc_over_ctrl", "auc_depth_to_disease"):
        log("   %-26s %s" % (k, r[k]))
    log("   pc_tech_share=%.3f pc_disease_share=%.3f gap=%.3f | 名义 %s/%s p<0.05"
        % (r["pc_tech_share"], r["pc_disease_share"], r["main_axis_gap"],
           r["n_ct_p05_perm"], r["n_ct_tested"]))
    recs.append(r)

T = pd.DataFrame(recs)

bx = batch_axis()
T["cell_source"] = T.dataset.map({"ds1": "CELLxGENE CxG_IBD_HC_10x (Garrido-Trigo 2024)",
                                 "ds2": "CELLxGENE CxG_Oliver2024_Myeloid",
                                 "ds3": "CELLxGENE CxG_Oliver2024_LargeIntestine",
                                 "ds4": "CELLxGENE 5bd4d642 (Tearle2026)"})
T["_akey"] = T.dataset.map(AUDIT_KEY)
T["batch_vars_available"] = T._akey.map(lambda d: bx.get(d, {}).get("vars", ""))
T["batch_best_var"] = T._akey.map(lambda d: bx.get(d, {}).get("best_var", ""))
T["batch_best_ctrl_share"] = T._akey.map(lambda d: bx.get(d, {}).get("best_share", np.nan))
T["batch_worst_var"] = T._akey.map(lambda d: bx.get(d, {}).get("worst_var", ""))
T["batch_worst_ctrl_share"] = T._akey.map(lambda d: bx.get(d, {}).get("worst_share", np.nan))
T["batch_nested_flag"] = T._akey.map(lambda d: bx.get(d, {}).get("nested_flag", np.nan))

# 细胞级深度口径（来自批次审计表）与 ds4 的来源审计列
_aud = {}
p = os.path.join(OUT, "partD_batch_audit_summary.tsv")
if os.path.exists(p):
    _b = pd.read_csv(p, sep="\t")
    _aud = {r["dataset"]: r for _, r in _b.iterrows()}
T["depth_fold_celllevel"] = T._akey.map(
    lambda d: float(_aud[d]["depth_fold_uc_over_ctrl"]) if d in _aud else np.nan)
T["frac_ctrl_cells_primary"] = T.dataset.map(lambda d: PROV["frac_ctrl_cells_primary"] if d == "ds4" else np.nan)
T["depth_fold_celllevel_primaryctrl"] = T.dataset.map(
    lambda d: PROV["depth_fold_celllevel_primaryctrl"] if d == "ds4" else np.nan)
T["n_ct_p05_after_primary_ctrl"] = T.dataset.map(
    lambda d: PROV["n_ct_p05_after_primary_ctrl"] if d == "ds4" else np.nan)
T["n_ct_in_primary_ctrl"] = T.dataset.map(
    lambda d: PROV["n_ct_in_primary_ctrl"] if d == "ds4" else np.nan)
T["n_ct_p05_after_within_donor"] = T.dataset.map(
    lambda d: PROV["n_ct_p05_after_within_donor"] if d == "ds4" else np.nan)
T["n_ct_in_within_donor"] = T.dataset.map(
    lambda d: PROV["n_ct_in_within_donor"] if d == "ds4" else np.nan)

order = ["dataset", "label", "cell_source", "n_cell_pseudobulk",
         "n_donor_total", "n_uc_donor", "n_ctrl_donor", "n_donor_excluded_other",
         "n_celltype_pb", "frac_ctrl_cells_primary",
         "med_nUMI_uc", "med_nUMI_ctrl", "depth_fold_uc_over_ctrl", "depth_fold_celllevel",
         "depth_fold_celllevel_primaryctrl", "auc_depth_to_disease", "depth_separability",
         "n_ct_for_rho", "rho_score_nUMI_med", "rho_score_nUMI_max",
         "n_ct_for_r2", "r2_depth_disease_med", "r2_depth_disease_max",
         "n_donor_in_pca", "n_ct_in_pca", "n_pc_used", "pc_var_expl_1",
         "pc1_r2_depth", "pc1_r2_depth_cw", "pc1_r2_disease",
         "pc_tech_share", "pc_tech_share_cw", "pc_disease_share", "main_axis_gap",
         "batch_vars_available", "batch_best_var", "batch_best_ctrl_share",
         "batch_worst_var", "batch_worst_ctrl_share", "batch_nested_flag",
         "n_ct_tested", "n_ct_p05_perm", "n_ct_q05", "min_p_perm",
         "n_ct_p05_after_rarefy", "n_ct_q05_after_rarefy", "min_p_after_rarefy",
         "n_ct_in_primary_ctrl", "n_ct_p05_after_primary_ctrl",
         "n_ct_in_within_donor", "n_ct_p05_after_within_donor"]
T = T[[c for c in order if c in T.columns]]

# 同步写一份 provenance 单列标注（ds4 的对照组成问题）
meta = dict(
    generated_by="scripts/partD_technical_inventory.py",
    definition=dict(
        depth_donorlevel="该 donor 各细胞类型中位 nUMI 的中位数（细胞类型等权）→ 组内再取中位；分析单元 = donor",
        depth_celllevel="batch_audit 的口径：全部细胞的 nUMI 中位（大 donor 主导），仅作对照",
        auc_depth_to_disease="donor 级 log10 深度区分 UC vs 对照的 AUC（>0.5 = UC 更深）",
        rho_score_nUMI="细胞类型内 Spearman(SGRG pseudobulk 评分, nUMI)，取 |ρ| 的中位/最大；需该类型 ≥5 donor",
        r2_depth_disease="细胞类型内 log10(nUMI) ~ 疾病 的一元回归 R²；需两组各 ≥2 donor",
        pc_tech_share="donor×细胞类型 SGRG 评分矩阵 PCA 后，PC1–3 中可被 donor 深度解释的方差比例 Σ ve_k·R²_k",
        pc_tech_share_cw="同上但深度改用细胞数加权的 donor 深度（稳健性对照；ds1 因 donor×类型仅 12×5，两口径可差 0.4）",
        pc_disease_share="同上，解释变量换成疾病分组",
        main_axis_gap="pc1_r2_depth − pc1_r2_disease（正 = 主轴更偏技术轴）",
        batch_best_ctrl_share="对照 donor 中与任一 UC donor 落在同一批次的最高比例"
                               "（跨所有可用批次变量取最大 = 对『批次混杂』这一指控最有利的口径）",
        batch_nested_flag="1 = 至少有一个批次变量下，多数对照 donor 完全落在不含病例的批次里（= 批次-疾病嵌套）",
        n_ct_p05_perm="SGRG 名义显著（置换检验 p<0.05）的细胞类型数；分母 n_ct_tested 取自同一张检验表",
        n_ct_p05_after_rarefy="二项抽稀到 1,000 UMI 后的名义显著数（仅 ds4 做了）",
        n_ct_p05_after_primary_ctrl="ds4 口径 B：对照只保留本研究自有（is_primary_data=True）的 2 个 donor",
        n_ct_p05_after_within_donor="ds4 口径 C：同一 donor 内炎症黏膜 vs 非炎症黏膜配对"),
    caveats=[
        "★ 深度比随口径变化很大，必须同时给出两个口径：ds4 细胞级 4.35× vs donor 级 1.77× vs 仅自有对照 1.04×；"
        "ds1 细胞级 1.52× vs donor 级 0.86×（方向反转）→ 『UC 测得更深』在多数数据集里是组成本身的效应。",
        "★ naive 口径的名义显著数受检验方法影响：置换检验在两组 pseudobulk 方差不等时 I 类错误膨胀"
        "（ds4 实测 8.5×）。与随机基因集比较必须用同一口径（见 partD_ds4_tearle_rand_summary.json）。",
        "ds4 的「对照」在 naive 口径下 94.4% 是外部复用数据（is_primary_data=False）。"
        "见 partD_ds4_prov_summary.json 与 logs/partD_ds4_prov.log。",
        "ds2/ds3 未做抽稀、未做来源审计，相关列留空而不是填 0。",
        "PCA 稳定性：ds1 只有 5 个细胞类型 × 12 donor，主轴指标（pc_*）稳定性弱于 ds4（33 类型 × 15 donor），引用时须并列 donor/类型数。",
    ],
)
T.to_csv(os.path.join(OUT, "partD_technical_inventory.tsv"), sep="\t", index=False)
json.dump(meta, open(os.path.join(OUT, "partD_technical_inventory_meta.json"), "w",
                     encoding="utf-8"), ensure_ascii=False, indent=2)
log("=" * 70)
log("[saved] results/partD_technical_inventory.tsv  %s" % (T.shape,))
print(T.drop(columns=["label", "cell_source"]).T.to_string())
