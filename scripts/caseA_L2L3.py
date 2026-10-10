# -*- coding: utf-8 -*-
"""
案例 A · L2（批次-疾病嵌套）+ L3（技术负载）+ 机制补充。

L2：用 Gemma 的 GEEQ / batch 信息刻画本队列的批次-疾病共线程度。
    若批次信息不足，退化为「技术轴占据度」：用全基因 SVD 的 PC1..k 与疾病标签的 AUC，
    以及铁死亡面板评分与 PC/RIN 等的关系。
L3：面板评分（score_genes 均值）与「样本技术深度代理」的关系。芯片无 nUMI，改用
    (a) 全基因平均表达（检测率代理）、(b) 每样本的探针检出数、(c) PC 与 panel score 的相关。
补充：单基因 AUC 分布 —— 若大量单个基因 AUC 就 >0.95，说明「面板」不比单基因强，
    这是「基因集无特异性」的直接证据。
"""
import os, json, sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.decomposition import PCA

ROOT = r"C:\Users\H-P\WorkBuddy\应激颗粒"
EXPR = os.path.join(ROOT, "results", "caseA_expr_7673.tsv")
RES = os.path.join(ROOT, "results")
TUMOR = "liver tumor"


def main():
    df = pd.read_csv(EXPR, sep="\t")
    genes = [c for c in df.columns if c not in ("sample", "label")]
    y = (df["label"] == TUMOR).astype(int).values
    X = df[genes].values.astype(float)
    X = np.nan_to_num(X, nan=np.nanmean(X, axis=0))
    out = {}

    # ---------- L3a 单基因 AUC 分布 ----------
    aucs = np.array([roc_auc_score(y, X[:, j]) for j in range(len(genes))])
    aucs_abs = np.maximum(aucs, 1 - aucs)  # 方向无关
    out["single_gene_auc"] = {
        "n": len(genes),
        "median": float(np.median(aucs_abs)),
        "max": float(aucs_abs.max()),
        "min": float(aucs_abs.min()),
        "n_ge_0.90": int((aucs_abs >= 0.90).sum()),
        "n_ge_0.95": int((aucs_abs >= 0.95).sum()),
        "n_ge_0.99": int((aucs_abs >= 0.99).sum()),
        "top5": sorted(zip(genes, aucs_abs.tolist()), key=lambda t: -t[1])[:5],
    }
    print("L3a single-gene AUC: median=%.3f max=%.3f  #>=0.95: %d/%d  #>=0.99: %d" % (
        np.median(aucs_abs), aucs_abs.max(), out["single_gene_auc"]["n_ge_0.95"],
        len(genes), out["single_gene_auc"]["n_ge_0.99"]))

    # ---------- L2 技术轴占据度：PC 与疾病的 AUC ----------
    Xz = (X - X.mean(0)) / (X.std(0) + 1e-9)
    n_comp = min(15, Xz.shape[1] - 1)
    pca = PCA(n_components=n_comp, random_state=0)
    P = pca.fit_transform(Xz)
    pc_auc = [roc_auc_score(y, P[:, k]) for k in range(n_comp)]
    pc_auc_abs = np.maximum(pc_auc, 1 - np.array(pc_auc))
    var = pca.explained_variance_ratio_
    out["L2_technical_axis"] = {
        "pc_var_ratio": var.tolist(),
        "pc_auc_disease_abs": pc_auc_abs.tolist(),
        "pc1_auc": float(pc_auc_abs[0]),
        "max_pc_auc": float(pc_auc_abs.max()),
        "max_pc_index": int(np.argmax(pc_auc_abs)),
        "sum_var_top3": float(var[:3].sum()),
    }
    print("L2 PC-disease AUC: PC1=%.3f  max(PC%d)=%.3f  var(top3)=%.1f%%" % (
        pc_auc_abs[0], int(np.argmax(pc_auc_abs)) + 1, pc_auc_abs.max(), 100 * var[:3].sum()))

    # ---------- L3b 面板评分 vs 技术代理 ----------
    panel_score = Xz.mean(axis=1)  # 标准化后的面板平均表达
    depth_proxy = X.mean(axis=1)   # 全基因平均表达（芯片上≈整体杂交强度）
    det_proxy = (X > X.mean(0)).sum(axis=1)  # 相对检出数
    out["L3_score_vs_technical"] = {
        "rho_panelscore_vs_meanexpr": float(stats.spearmanr(panel_score, depth_proxy).statistic),
        "rho_panelscore_vs_pc1": float(stats.spearmanr(panel_score, P[:, 0]).statistic),
        "auc_meanexpr_disease": float(max(roc_auc_score(y, depth_proxy),
                                          1 - roc_auc_score(y, depth_proxy))),
        "auc_panelscore_disease": float(max(roc_auc_score(y, panel_score),
                                            1 - roc_auc_score(y, panel_score))),
    }
    print("L3 rho(panel_score, mean_expr)=%.3f  AUC(mean_expr~disease)=%.3f" % (
        out["L3_score_vs_technical"]["rho_panelscore_vs_meanexpr"],
        out["L3_score_vs_technical"]["auc_meanexpr_disease"]))

    # ---------- 机制：去掉全局轴后 AUC 还剩多少 ----------
    # 用 PC1 残差化（把疾病主轴投影掉）再重跑面板 logistic
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    Xres = Xz - np.outer(P[:, 0], pca.components_[0])
    pv = stats.ttest_ind(Xres[y == 1], Xres[y == 0], axis=0, equal_var=True)[1]
    sel = np.where(pv < 0.001)[0]
    if len(sel) >= 1:
        lr = LogisticRegression(max_iter=5000)
        lr.fit(Xres[:, sel], y)
        a_res = roc_auc_score(y, lr.predict_proba(Xres[:, sel])[:, 1])
    else:
        a_res = np.nan
    out["after_removing_global_axis"] = {
        "n_de": int(len(sel)), "auc_residual": float(a_res) if a_res == a_res else None}
    print("after removing PC1: n_de=%d  AUC=%.4f" % (len(sel), a_res))

    json.dump(out, open(os.path.join(RES, "caseA_L2L3.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("wrote results/caseA_L2L3.json")


if __name__ == "__main__":
    main()
