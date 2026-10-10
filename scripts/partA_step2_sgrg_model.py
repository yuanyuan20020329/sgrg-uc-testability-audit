# -*- coding: utf-8 -*-
"""
Part A · Step 2  —  UC 应激颗粒相关基因(SGRG) 诊断模型
========================================================================
设计：双向交叉队列 (bidirectional cross-cohort) 验证
  方向 1 : 训练 GSE11223 (UC vs 严格健康对照)  → 外部验证 GSE38713
  方向 2 : 训练 GSE38713                        → 外部验证 GSE11223

特征  : 两队列共同可定量的 SGRG（HGNC 批准符号）
模型  : L1 / L2 惩罚逻辑回归；正则强度 C 由训练集内部交叉验证择优（测试集不参与）
评价  : 训练集交叉验证 AUC（n<=60 用 LOO，否则 5 折分层）＋ 外部队列独立 AUC
        标签置换检验给出经验 p
标准化 : 队列内逐基因 z-score（跨平台迁移的常规妥协，见报告局限）
"""
import os, re, json, sys, warnings
from collections import Counter
import numpy as np
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sgrg_common as C

OUT = C.OUT
log = C.log
RNG = np.random.default_rng(20261008)

# ================================================================ 队列装配
class Cohort:
    def __init__(self, name, expr, groups):
        self.name = name; self.expr = expr; self.groups = groups

def build_cohort(name, matrix, lab2g):
    header, rows = C.load_matrix(matrix)
    sym_col = C.find_sym_col(header)
    sample_cols = list(range(6, len(header)))          # Probe/Sequence/GeneSymbol/GeneName/GemmaId/NCBIid
    keys_per_col = [C.col_keys(header[j]) for j in sample_cols]
    expr, _ = C.collapse_genes(header, rows, sym_col, sample_cols)
    groups = [C.resolve_label(ks, lab2g) for ks in keys_per_col]
    return Cohort(name, expr, groups)

# ================================================================ CV 工具
def make_pipe(C_):
    return Pipeline([("sc", StandardScaler()),
                     ("lr", LogisticRegression(C=C_, solver="liblinear", max_iter=5000))])

def cv_predict_proba(X, y, C_, penalty, seed=0):
    pipe = make_pipe(C_)
    pipe.named_steps["lr"].set_params(penalty=penalty)
    n = X.shape[0]
    cv = LeaveOneOut() if n <= 60 else StratifiedKFold(5, shuffle=True, random_state=seed)
    oof = np.full(n, np.nan)
    for tr, te in cv.split(X, y):
        pipe.fit(X[tr], y[tr])
        oof[te] = pipe.predict_proba(X[te])[:, 1]
    return oof

def tune_C(X, y, penalties=("l1", "l2"),
           grid=(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 3.0)):
    per = {}
    for pen in penalties:
        for C_ in grid:
            oof = cv_predict_proba(X, y, C_, pen)
            if np.isnan(oof).any():
                continue
            auc = roc_auc_score(y, oof)
            if pen not in per or auc > per[pen][0]:
                per[pen] = (auc, C_, oof)
    best_pen = max(per, key=lambda k: per[k][0])
    return per, best_pen, per[best_pen]

def metrics(y, score):
    fpr, tpr, thrs = roc_curve(y, score)
    auc = roc_auc_score(y, score)
    thr = thrs[int(np.argmax(tpr - fpr))]
    pred = (score >= thr).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    return dict(auc=float(auc), thr=float(thr),
                acc=float((tp + tn) / (tp + tn + fp + fn)),
                sens=float(tp / (tp + fn)) if (tp + fn) else float("nan"),
                spec=float(tn / (tn + fp)) if (tn + fp) else float("nan"),
                tp=int(tp), fp=int(fp), tn=int(tn), fn=int(fn)), (fpr, tpr)

def permutation_p(X, y, observed_auc, n_perm=200, C_=1.0, penalty="l2"):
    rng = np.random.default_rng(7); cnt = 0
    for _ in range(n_perm):
        yp = rng.permutation(y)
        try:
            oof = cv_predict_proba(X, yp, C_, penalty)
            if np.isnan(oof).any(): continue
            if roc_auc_score(yp, oof) >= observed_auc: cnt += 1
        except Exception:
            continue
    return (cnt + 1) / (n_perm + 1)

# ================================================================ main
def main():
    sgrg = C.load_sgrg()
    log(f"[SGRG] panel = {len(sgrg)}")

    specs = [("GSE11223", C.MAT_11223, C.groups_gse11223()),
             ("GSE38713", C.MAT_38713, C.groups_gse38713())]
    cohorts = {}
    for nm, mat, lab2g in specs:
        if not os.path.exists(mat):
            log(f"[skip] {nm}: matrix missing"); continue
        co = build_cohort(nm, mat, lab2g)
        if len(co.expr) < 1000:
            log(f"[skip] {nm}: 矩阵不完整（仅解析出 {len(co.expr)} 个基因，疑为下载截断）")
            continue
        log(f"[cohort] {nm}: groups={dict(Counter(co.groups))} | genes={len(co.expr)}")
        cohorts[nm] = co

    common = None
    for nm, co in cohorts.items():
        meas = set(co.expr) & set(sgrg)
        common = meas if common is None else (common & meas)
    feats = sorted(common)
    log(f"[features] common measurable SGRG = {len(feats)}")

    def matrix_of(co):
        i_uc = [i for i, g in enumerate(co.groups) if g.startswith("UC")]
        i_ct = [i for i, g in enumerate(co.groups) if g == "Control"]
        use = i_uc + i_ct
        X = np.vstack([co.expr[g][use] for g in feats]).T
        y = np.array([1] * len(i_uc) + [0] * len(i_ct))
        return X, y, use

    data = {nm: matrix_of(co) for nm, co in cohorts.items()}
    for nm, (X, y, _) in data.items():
        log(f"[data] {nm}: X={X.shape} | pos={int(y.sum())} neg={int((y == 0).sum())}")

    results, roc_panels = {}, []
    for tr_name, co in cohorts.items():
        Xtr, ytr, _ = data[tr_name]
        per, pen, (auc_cv, Cbest, oof) = tune_C(Xtr, ytr)
        per_str = " | ".join(f"{k.upper()} C={v[1]} AUC={v[0]:.3f}" for k, v in per.items())
        log(f"[tune] {tr_name}: best {pen.upper()} C={Cbest} | CV AUC={auc_cv:.3f}  [{per_str}]")
        m_tr, (fpr_tr, tpr_tr) = metrics(ytr, oof)
        p_perm = permutation_p(Xtr, ytr, auc_cv, n_perm=200, C_=Cbest, penalty=pen)

        e = dict(train=tr_name, n_train=int(len(ytr)), n_pos=int(ytr.sum()),
                 best_C=Cbest, penalty=pen,
                 cv_auc=m_tr["auc"], cv_acc=m_tr["acc"], cv_sens=m_tr["sens"], cv_spec=m_tr["spec"],
                 permutation_p=p_perm,
                 cv_auc_l1=round(per["l1"][0], 4) if "l1" in per else None,
                 cv_auc_l2=round(per["l2"][0], 4) if "l2" in per else None)
        roc_panels.append((f"{tr_name} train CV ({pen.upper()} C={Cbest}) AUC={m_tr['auc']:.3f}",
                           fpr_tr, tpr_tr, m_tr["auc"], "#0b5cab"))

        te_names = [n for n in cohorts if n != tr_name]
        if te_names:
            tn = te_names[0]
            Xte, yte, use_te = data[tn]
            pipe = make_pipe(Cbest); pipe.named_steps["lr"].set_params(penalty=pen)
            pipe.fit(Xtr, ytr)
            # 跨平台：测试队列内部逐基因标准化后再套用训练系数
            Xte_z = (Xte - Xte.mean(0)) / (Xte.std(0) + 1e-9)
            score_te = pipe.decision_function(Xte_z)
            m_te, (fpr_te, tpr_te) = metrics(yte, score_te)
            e.update(dict(test=tn, n_test=int(len(yte)), n_pos_test=int(yte.sum()),
                          ext_auc=m_te["auc"], ext_acc=m_te["acc"], ext_sens=m_te["sens"],
                          ext_spec=m_te["spec"], ext_tp=m_te["tp"], ext_fp=m_te["fp"],
                          ext_tn=m_te["tn"], ext_fn=m_te["fn"]))
            log(f"[ext] {tr_name} -> {tn}: AUC={m_te['auc']:.3f} acc={m_te['acc']:.3f} "
                f"sens={m_te['sens']:.3f} spec={m_te['spec']:.3f}")
            roc_panels.append((f"{tn} external (trained on {tr_name}) AUC={m_te['auc']:.3f}",
                               fpr_te, tpr_te, m_te["auc"], "#b42318"))
            e["_score"] = score_te.tolist()
            e["_grp"] = [cohorts[tn].groups[i] for i in use_te]

        # 签名：L1 在其最优 C 上全量拟合
        sig_pen = "l1" if "l1" in per else pen
        sig_C = per[sig_pen][1]
        pf = make_pipe(sig_C); pf.named_steps["lr"].set_params(penalty=sig_pen)
        pf.fit(Xtr, ytr)
        coef = pf.named_steps["lr"].coef_.ravel()
        order = np.argsort(-np.abs(coef))
        sig = [dict(gene=feats[i], coef=float(coef[i])) for i in order if abs(coef[i]) > 1e-8]
        e["signature"] = sig
        e["signature_penalty"] = sig_pen
        e["signature_C"] = sig_C
        e["signature_n"] = len(sig)
        log(f"[signature] {tr_name}: {sig_pen.upper()} C={sig_C} -> {len(sig)} genes | "
            + ", ".join(f"{s['gene']}({s['coef']:+.2f})" for s in sig[:8]))
        results[tr_name] = e

    # ---------------- ROC ----------------
    fig, ax = plt.subplots(figsize=(5.8, 5.3), dpi=300)
    for lab, fpr, tpr, auc, col in roc_panels:
        ls = "-" if "train CV" in lab else "--"
        ax.plot(fpr, tpr, lw=1.7, ls=ls, color=col, label=lab)
    ax.plot([0, 1], [0, 1], ":", lw=.9, color="#8a9199")
    ax.set_xlabel("False positive rate", fontsize=9)
    ax.set_ylabel("True positive rate", fontsize=9)
    ax.set_title("SGRG model · ROC (solid = train CV, dashed = external)", fontsize=9.5)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.tick_params(labelsize=8)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    fig.tight_layout()
    p_roc = os.path.join(OUT, "fig_ROC_SGRG_diagnostic_model.png")
    fig.savefig(p_roc, dpi=300); plt.close(fig); log("[write]", p_roc)

    # ---------------- 风险评分 ----------------
    for tr_name, e in results.items():
        if "_score" not in e: continue
        sc = np.array(e["_score"]); grp = np.array(e["_grp"])
        order = ["Control", "UC_remission", "UC_active_noninflamed", "UC_active_inflamed", "UC"]
        present = [g for g in order if (grp == g).any()]
        present += [g for g in sorted(set(grp)) if g not in present]
        cols = {"Control": "#7f9db9", "UC_remission": "#c9b458",
                "UC_active_noninflamed": "#e08a3c", "UC_active_inflamed": "#b42318", "UC": "#b42318"}
        fig, ax = plt.subplots(figsize=(6.2, 4.2), dpi=300)
        vals = [sc[grp == g] for g in present]
        bp = ax.boxplot(vals, patch_artist=True, widths=.55, medianprops=dict(color="#111", lw=1.2))
        for patch, g in zip(bp["boxes"], present):
            patch.set_facecolor(cols.get(g, "#9aa0a6")); patch.set_alpha(.75)
        for i, g in enumerate(present):
            ax.scatter(np.full(len(vals[i]), i + 1) + RNG.normal(0, .06, len(vals[i])),
                       vals[i], s=15, color="#333", alpha=.6, lw=0, zorder=3)
        ax.set_xticks(range(1, len(present) + 1))
        ax.set_xticklabels([f"{g}\n(n={len(sc[grp==g])})" for g in present], fontsize=7.5)
        ax.axhline(0, ls="--", lw=.8, color="#999")
        ax.set_ylabel("SGRG model decision score", fontsize=9)
        ax.set_title(f"{e.get('test','?')} (external cohort) · risk score by mucosal group\n"
                     f"model trained on {tr_name}", fontsize=9)
        ax.tick_params(labelsize=8)
        for s in ["top", "right"]: ax.spines[s].set_visible(False)
        fig.tight_layout()
        pr = os.path.join(OUT, f"fig_riskscore_by_group_{e.get('test','x')}.png")
        fig.savefig(pr, dpi=300); plt.close(fig); log("[write]", pr)
        pt = os.path.join(OUT, f"partA_step2_riskscore_{e.get('test','x')}_extby_{tr_name}.tsv")
        with open(pt, "w", encoding="utf-8") as f:
            f.write("Group\tDecisionScore\n")
            for g, v in zip(grp, sc): f.write(f"{g}\t{v:.4f}\n")
        results[tr_name]["riskscore_tsv"] = os.path.basename(pt)

    # ---------------- 签名图 ----------------
    for tr_name, e in results.items():
        sig = e["signature"][:25]
        if not sig: continue
        fig, ax = plt.subplots(figsize=(5.8, max(3.0, .27 * len(sig))), dpi=300)
        gs = [s["gene"] for s in sig][::-1]; cv = [s["coef"] for s in sig][::-1]
        ax.barh(range(len(sig)), cv, color=["#b42318" if c > 0 else "#0b5cab" for c in cv], height=.72)
        ax.set_yticks(range(len(sig))); ax.set_yticklabels(gs, fontsize=7.5)
        ax.axvline(0, color="#444", lw=.8)
        ax.set_xlabel("standardized coefficient (log-odds per SD)", fontsize=8.5)
        ax.set_title(f"SGRG signature · trained on {tr_name}\n"
                     f"({e['signature_penalty'].upper()}, C={e['signature_C']}; "
                     f"{len(e['signature'])} of {len(feats)} genes retained)", fontsize=9)
        ax.tick_params(labelsize=8)
        for s in ["top", "right"]: ax.spines[s].set_visible(False)
        fig.tight_layout()
        ps = os.path.join(OUT, f"fig_signature_{tr_name}.png")
        fig.savefig(ps, dpi=300); plt.close(fig); log("[write]", ps)

    # ---------------- 落盘 ----------------
    clean = {k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")} for k, v in results.items()}
    pj = os.path.join(OUT, "partA_step2_model_metrics.json")
    json.dump(clean, open(pj, "w", encoding="utf-8"), ensure_ascii=False, indent=2); log("[write]", pj)
    pst = os.path.join(OUT, "partA_step2_signatures.tsv")
    with open(pst, "w", encoding="utf-8") as f:
        f.write("train_cohort\tGeneSymbol\tcoef\n")
        for k, e in results.items():
            for s in e["signature"]:
                f.write(f"{k}\t{s['gene']}\t{s['coef']:.5f}\n")
    log("[write]", pst)
    log("[done]")

if __name__ == "__main__":
    main()
