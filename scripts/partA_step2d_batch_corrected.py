# -*- coding: utf-8 -*-
"""
Part A · Step 2d — 去批次后 SGRG 面板判别力的变化
做法：逐基因减去其所在 batch 的均值（最简单可解释的批次中心化），
      再比较 SGRG 面板与随机同规模面板的交叉验证 AUC。
注意：当疾病状态与批次混杂时，去批次会连同被混杂的那部分"疾病信号"一起移除，
      因此本检验回答的是"若不依赖批次，还剩多少判别力"，而不是"真实效应有多大"。
"""
import os, sys, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sgrg_common as C

OUT = C.OUT; log = C.log
RNG = np.random.default_rng(999)
C_FIX, PEN_FIX = 0.01, "l2"

def cv_auc(X, y, C_=C_FIX, penalty=PEN_FIX, seed=0):
    pipe = Pipeline([("sc", StandardScaler()),
                     ("lr", LogisticRegression(C=C_, penalty=penalty, solver="liblinear", max_iter=5000))])
    n = X.shape[0]
    cv = LeaveOneOut() if n <= 60 else StratifiedKFold(5, shuffle=True, random_state=seed)
    oof = np.full(n, np.nan)
    for tr, te in cv.split(X, y):
        pipe.fit(X[tr], y[tr]); oof[te] = pipe.predict_proba(X[te])[:, 1]
    return roc_auc_score(y, oof)

def main():
    matrix = C.MAT_38713
    header, rows = C.load_matrix(matrix)
    sym_col = C.find_sym_col(header)
    sample_cols = list(range(6, len(header)))
    expr, _ = C.collapse_genes(header, rows, sym_col, sample_cols)
    grp = np.array([C.resolve_label(C.col_keys(header[j]), C.groups_gse38713()) for j in sample_cols])
    bmap = C.batch_gse38713()
    batch = np.array([C.resolve_label(C.col_keys(header[j]), bmap) for j in sample_cols])

    is_uc = np.array([str(g).startswith("UC") for g in grp])
    is_ct = np.array([str(g) == "Control" for g in grp])
    keep = is_uc | is_ct
    y = np.array([1 if u else 0 for u in is_uc[keep]])
    bt = batch[keep]
    log(f"[samples] n={keep.sum()} (UC {y.sum()} / Ctrl {(y==0).sum()}) | batches={sorted(set(bt))}")

    genes = [g for g, v in expr.items() if np.isfinite(v).mean() >= .9]
    G = np.vstack([expr[g][keep] for g in genes])          # genes x samples
    sgrg = set(C.load_sgrg())
    feat = [g for g in genes if g in sgrg]
    idx = {g: i for i, g in enumerate(genes)}
    log(f"[genes] pool={len(genes)} | SGRG={len(feat)}")

    def run_variant(tag, Xall):
        Xs = Xall[[idx[g] for g in feat]].T
        a_s = cv_auc(Xs, y)
        others = [g for g in genes if g not in sgrg]
        aucs = []
        for i in range(200):
            sel = RNG.choice(len(others), size=len(feat), replace=False)
            X = Xall[[idx[others[j]] for j in sel]].T
            aucs.append(cv_auc(X, y, seed=i))
        aucs = np.array(aucs)
        pct = float((aucs < a_s).mean() * 100)
        log(f"[{tag}] SGRG AUC={a_s:.3f} | random median={np.median(aucs):.3f} "
            f"p95={np.percentile(aucs,95):.3f} | percentile={pct:.1f}")
        return dict(tag=tag, auc_sgrg=round(float(a_s), 4),
                    rand_median=round(float(np.median(aucs)), 4),
                    rand_p95=round(float(np.percentile(aucs, 95)), 4),
                    percentile=round(pct, 1)), aucs, a_s

    raw = run_variant("raw", G)
    # 去批次：逐基因减所在 batch 均值
    Gc = G.copy()
    for b in sorted(set(bt)):
        m = bt == b
        if m.sum() >= 2:
            Gc[:, m] = Gc[:, m] - G[:, m].mean(axis=1, keepdims=True)
    corr = run_variant("batch_centered", Gc)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9), dpi=300, sharey=True)
    for ax, (res, aucs, a_s, title) in zip(
            axes, [(raw[0], raw[1], raw[2], "raw (as reported)"),
                   (corr[0], corr[1], corr[2], "after batch centering")]):
        ax.hist(aucs, bins=28, color="#9db4c9", edgecolor="#fff", lw=.5)
        ax.axvline(a_s, color="#b42318", lw=2.0, label=f"SGRG={a_s:.3f}")
        ax.axvline(np.median(aucs), color="#5b6470", ls="--", lw=1.0,
                   label=f"random median={np.median(aucs):.3f}")
        ax.set_title(f"GSE38713 · {title}", fontsize=9.5)
        ax.set_xlabel("cross-validated AUC", fontsize=9)
        ax.legend(fontsize=7.5, frameon=False)
        ax.tick_params(labelsize=8)
        for s in ["top", "right"]: ax.spines[s].set_visible(False)
    axes[0].set_ylabel("count", fontsize=9)
    fig.tight_layout()
    p = os.path.join(OUT, "fig_batch_corrected_specificity_GSE38713.png")
    fig.savefig(p, dpi=300); plt.close(fig); log("[write]", p)

    out = dict(dataset="GSE38713", n_samples=int(keep.sum()), n_uc=int(y.sum()),
               n_ctrl=int((y == 0).sum()), n_sgrg=len(feat), gene_pool=len(genes),
               raw=raw[0], batch_centered=corr[0],
               note="batch centering 在疾病-批次混杂时会一并移除被混杂的疾病信号，故为下界")
    pj = os.path.join(OUT, "partA_step2d_batch_corrected.json")
    json.dump(out, open(pj, "w", encoding="utf-8"), ensure_ascii=False, indent=2); log("[write]", pj)
    log("[done]")

if __name__ == "__main__":
    main()
