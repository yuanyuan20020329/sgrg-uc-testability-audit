# -*- coding: utf-8 -*-
"""
Part A · Step 2c — 全局结构与批次混杂诊断
1) 数据里是否存在与疾病状态对齐的全局主轴（PC1）？解释多少方差？AUC 多少？
2) GSE38713 的 batch 因子与疾病状态混杂到什么程度？（列联表 / Cramér's V / 批次单变量预测疾病）
3) 用 batch 标签本身作为特征，能多好地预测疾病状态？（若很高 → 诊断信号可疑）
产出 fig_global_structure_*.png / fig_batch_confound_GSE38713.png 及 JSON
"""
import os, sys, json, warnings
import numpy as np
from collections import Counter
import warnings as _w
_w.filterwarnings("ignore")
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, cross_val_predict
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OneHotEncoder
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sgrg_common as C

OUT = C.OUT
log = C.log

def cramers_v(tab):
    """手写 Cramér's V，避免 chi2_contingency 对零期望频数报错"""
    tab = np.asarray(tab, float)
    n = tab.sum()
    exp = np.outer(tab.sum(1), tab.sum(0)) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = float(np.nansum(np.where(exp > 0, (tab - exp) ** 2 / exp, 0.0)))
    r, k = tab.shape
    dof = (r - 1) * (k - 1)
    p = float(stats.chi2.sf(chi2, dof)) if dof > 0 else float("nan")
    denom = n * (min(r, k) - 1)
    return (float(np.sqrt(chi2 / denom)) if denom > 0 else float("nan")), p, chi2, dof

def build(name, matrix, lab2g):
    header, rows = C.load_matrix(matrix)
    sym_col = C.find_sym_col(header)
    sample_cols = list(range(6, len(header)))
    grp = np.array([C.resolve_label(C.col_keys(header[j]), lab2g) for j in sample_cols])
    M = []
    for r in rows:
        try: v = np.array([float(r[c]) for c in sample_cols], float)
        except (ValueError, IndexError): continue
        if np.isfinite(v).mean() >= .9: M.append(v)
    M = np.vstack(M)
    return M, grp

def run(name, matrix, lab2g, batch2=None):
    M, grp = build(name, matrix, lab2g)
    if M.shape[0] < 1000:
        log(f"[{name}] 矩阵不完整（{M.shape[0]} 行），跳过"); return None
    log(f"[{name}] genes={M.shape[0]} samples={M.shape[1]} | {dict(Counter(grp))}")
    is_uc = np.array([str(g).startswith("UC") for g in grp])
    is_ct = np.array([str(g) == "Control" for g in grp])
    keep = is_uc | is_ct
    X = M[:, keep].T
    y = np.array([1 if u else 0 for u in is_uc[keep]])
    disp = np.array(["UC" if u else "Control" for u in is_uc[keep]])
    Xz = (X - X.mean(0)) / (X.std(0) + 1e-9)
    U, S, _ = np.linalg.svd(Xz, full_matrices=False)
    pc = U * S
    var = S**2 / np.sum(S**2)
    aucs = [round(max(roc_auc_score(y, pc[:, k]), 1 - roc_auc_score(y, pc[:, k])), 4)
            for k in range(min(5, pc.shape[1]))]
    rho, p = stats.spearmanr(pc[:, 0], np.arange(len(y)))
    log(f"[{name}] PC1 var={var[0]:.3f} | AUC PC1..5={aucs} | PC1~order rho={rho:.3f}")

    res = dict(dataset=name, n_genes=int(M.shape[0]), n_samples=int(keep.sum()),
               n_uc=int(y.sum()), n_ctrl=int((y == 0).sum()),
               pc1_var=round(float(var[0]), 4), pc_auc_vs_disease=aucs,
               pc1_order_spearman=round(float(rho), 4))

    # ---------------- 批次混杂（GSE38713） ----------------
    batch = None
    if batch2 is not None:
        hdr = C.load_matrix(matrix)[0]
        bl = [C.resolve_label(C.col_keys(hdr[6 + i]), batch2) for i in range(len(hdr) - 6)]
        batch = np.array(bl)[keep]
        log(f"[{name}] batch: {dict(Counter(batch))}")
        cats = sorted(set(batch))
        tab = np.array([[int(((batch == b) & (y == 1)).sum()), int(((batch == b) & (y == 0)).sum())]
                        for b in cats])
        v, pv, chi2, dof = cramers_v(tab)
        log(f"[{name}] 批次×疾病 Cramér's V={v:.3f} chi2={chi2:.1f} dof={dof} p={pv:.3g}")
        for b, row in zip(cats, tab):
            log(f"    {b}:  UC={row[0]}  Ctrl={row[1]}")
        # 用 batch one-hot 预测疾病（交叉验证 AUC）
        oh = OneHotEncoder(sparse_output=False).fit_transform(batch.reshape(-1, 1))
        n = len(y)
        cv = LeaveOneOut() if n <= 60 else StratifiedKFold(5, shuffle=True, random_state=0)
        clf = LogisticRegression(max_iter=2000)
        try:
            prob = cross_val_predict(clf, oh, y, cv=cv, method="predict_proba")[:, 1]
            auc_batch = float(roc_auc_score(y, prob))
        except Exception:
            auc_batch = float("nan")
        log(f"[{name}] 仅用 batch 预测疾病 CV AUC = {auc_batch:.3f}")
        res.update(batch_cramers_v=round(v, 4), batch_chi2=round(float(chi2), 2),
                   batch_dof=int(dof), batch_chi2_p=float(pv),
                   batch_cats={b: dict(UC=int(r0), Control=int(r1)) for b, (r0, r1) in zip(cats, tab)},
                   batch_only_auc=round(auc_batch, 4))

        # 批次图
        fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), dpi=300)
        ax = axes[0]
        idx = np.arange(len(cats)); w = .38
        ax.bar(idx - w / 2, tab[:, 0], w, color="#b42318", label="UC")
        ax.bar(idx + w / 2, tab[:, 1], w, color="#7f9db9", label="Control")
        ax.set_xticks(idx); ax.set_xticklabels([c.replace("Batch_", "B") for c in cats],
                                               rotation=45, ha="right", fontsize=7.5)
        ax.set_ylabel("# samples", fontsize=9)
        ax.set_title(f"{name}: disease × batch\nCramér's V={v:.2f}", fontsize=9.5)
        ax.legend(fontsize=8, frameon=False); ax.tick_params(labelsize=8)
        ax = axes[1]
        for b in cats:
            m = batch == b
            ax.scatter(np.full(m.sum(), np.where(np.array(cats) == b)[0][0]) + np.random.default_rng(2).normal(0, .06, m.sum()),
                       pc[keep][m, 0], s=20, alpha=.8, lw=0)
        ax.set_xticks(range(len(cats)))
        ax.set_xticklabels([c.replace("Batch_", "B") for c in cats], rotation=45, ha="right", fontsize=7.5)
        ax.set_ylabel("PC1 score", fontsize=9)
        ax.set_title("PC1 anchored to batch", fontsize=9.5)
        ax.tick_params(labelsize=8)
        for a in axes:
            for s in ["top", "right"]: a.spines[s].set_visible(False)
        fig.tight_layout()
        pb = os.path.join(OUT, f"fig_batch_confound_{name}.png")
        fig.savefig(pb, dpi=300); plt.close(fig); log("[write]", pb)
        res["batch_fig"] = os.path.basename(pb)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9), dpi=300)
    ax = axes[0]
    v10 = var[:10] * 100
    ax.bar(range(1, len(v10) + 1), v10, color="#0b5cab", alpha=.85)
    ax.set_xlabel("principal component", fontsize=9)
    ax.set_ylabel("% variance explained", fontsize=9)
    ax.set_title(f"{name} · scree (PC1 = {var[0]*100:.1f}%)", fontsize=9.5)
    ax.tick_params(labelsize=8)
    ax = axes[1]
    for g, c in [("Control", "#7f9db9"), ("UC", "#b42318")]:
        v = pc[keep][disp == g, 0]
        ax.scatter(np.full(len(v), 0 if g == "Control" else 1) + np.random.default_rng(1).normal(0, .05, len(v)),
                   v, s=22, color=c, alpha=.8, lw=0, label=f"{g} (n={len(v)})")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Control", "UC"], fontsize=9)
    ax.set_ylabel("PC1 score", fontsize=9)
    ax.set_title(f"PC1 vs disease · AUC={aucs[0]:.3f}", fontsize=9.5)
    ax.legend(fontsize=8, frameon=False); ax.tick_params(labelsize=8)
    for a in axes:
        for s in ["top", "right"]: a.spines[s].set_visible(False)
    fig.tight_layout()
    pth = os.path.join(OUT, f"fig_global_structure_{name}.png")
    fig.savefig(pth, dpi=300); plt.close(fig); log("[write]", pth)
    res["global_fig"] = os.path.basename(pth)
    return res

if __name__ == "__main__":
    out = {}
    if os.path.exists(C.MAT_38713):
        r = run("GSE38713", C.MAT_38713, C.groups_gse38713(), batch2=C.batch_gse38713())
        if r: out["GSE38713"] = r
    if os.path.exists(C.MAT_11223):
        r = run("GSE11223", C.MAT_11223, C.groups_gse11223())
        if r: out["GSE11223"] = r
    p = os.path.join(OUT, "partA_step2c_global_structure.json")
    json.dump(out, open(p, "w", encoding="utf-8"), ensure_ascii=False, indent=2); log("[write]", p)
    log("[done]")
