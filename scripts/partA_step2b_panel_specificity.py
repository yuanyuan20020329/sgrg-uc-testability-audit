# -*- coding: utf-8 -*-
"""
Part A · Step 2b — 面板特异性检验
同规模随机基因面板 vs SGRG 面板，比较交叉验证 AUC。
回答："这个 AUC 是 SGRG 特有的，还是任何同规模基因集都能做到？"
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

OUT = C.OUT
log = C.log
RNG = np.random.default_rng(12345)
C_FIX, PEN_FIX = 0.01, "l2"      # 与 Step2 在 GSE38713 上选出的最优配置一致，保证可比

def cv_auc(X, y, C_=C_FIX, penalty=PEN_FIX, seed=0):
    pipe = Pipeline([("sc", StandardScaler()),
                     ("lr", LogisticRegression(C=C_, penalty=penalty, solver="liblinear", max_iter=5000))])
    n = X.shape[0]
    cv = LeaveOneOut() if n <= 60 else StratifiedKFold(5, shuffle=True, random_state=seed)
    oof = np.full(n, np.nan)
    for tr, te in cv.split(X, y):
        pipe.fit(X[tr], y[tr]); oof[te] = pipe.predict_proba(X[te])[:, 1]
    return roc_auc_score(y, oof)

def run(name, matrix, lab2g, n_rand=200):
    header, rows = C.load_matrix(matrix)
    sym_col = C.find_sym_col(header)
    sample_cols = list(range(6, len(header)))
    expr, _ = C.collapse_genes(header, rows, sym_col, sample_cols)
    if len(expr) < 1000:
        log(f"[{name}] 矩阵不完整（{len(expr)} 基因），跳过"); return None
    grp = [C.resolve_label(C.col_keys(header[j]), lab2g) for j in sample_cols]
    sgrg = set(C.load_sgrg())

    i_uc = [i for i, g in enumerate(grp) if g.startswith("UC")]
    i_ct = [i for i, g in enumerate(grp) if g == "Control"]
    use = i_uc + i_ct
    y = np.array([1] * len(i_uc) + [0] * len(i_ct))

    pool = sorted(g for g, v in expr.items() if np.isfinite(v).mean() >= .9)
    feat = [g for g in pool if g in sgrg]
    log(f"[{name}] n={len(use)} (UC {len(i_uc)} / Ctrl {len(i_ct)}) | pool={len(pool)} | SGRG={len(feat)}")

    Xs = np.vstack([expr[g][use] for g in feat]).T
    auc_sgrg = cv_auc(Xs, y)
    log(f"[{name}] SGRG panel CV AUC = {auc_sgrg:.3f}")

    others = [g for g in pool if g not in sgrg]
    k = len(feat); aucs = []
    for i in range(n_rand):
        sel = RNG.choice(len(others), size=k, replace=False)
        X = np.vstack([expr[others[j]][use] for j in sel]).T
        aucs.append(cv_auc(X, y, seed=i))
    aucs = np.array(aucs)
    pct = float((aucs < auc_sgrg).mean() * 100)
    log(f"[{name}] random({n_rand}): median={np.median(aucs):.3f} p95={np.percentile(aucs,95):.3f} "
        f"max={aucs.max():.3f} | SGRG percentile={pct:.1f}")

    fig, ax = plt.subplots(figsize=(6.0, 4.0), dpi=300)
    ax.hist(aucs, bins=28, color="#9db4c9", edgecolor="#fff", lw=.5)
    ax.axvline(auc_sgrg, color="#b42318", lw=2.0, label=f"SGRG panel AUC={auc_sgrg:.3f}")
    ax.axvline(np.median(aucs), color="#5b6470", ls="--", lw=1.0, label=f"random median={np.median(aucs):.3f}")
    ax.set_xlabel(f"cross-validated AUC  ({n_rand} random {k}-gene panels)", fontsize=9)
    ax.set_ylabel("count", fontsize=9)
    ax.set_title(f"{name} · SGRG panel vs random gene panels", fontsize=10)
    ax.legend(fontsize=8, frameon=False); ax.tick_params(labelsize=8)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    fig.tight_layout()
    p = os.path.join(OUT, f"fig_panel_specificity_{name}.png")
    fig.savefig(p, dpi=300); plt.close(fig); log("[write]", p)

    return dict(dataset=name, n_samples=int(len(use)), n_uc=len(i_uc), n_ctrl=len(i_ct),
                n_features=k, gene_pool=len(pool), auc_sgrg=round(float(auc_sgrg), 4),
                rand_median=round(float(np.median(aucs)), 4),
                rand_p95=round(float(np.percentile(aucs, 95)), 4),
                rand_max=round(float(aucs.max()), 4), n_random=n_rand,
                sgrg_percentile=round(pct, 1), model=f"{PEN_FIX.upper()} C={C_FIX}")

if __name__ == "__main__":
    out = {}
    for nm, mat, gf in [("GSE38713", C.MAT_38713, C.groups_gse38713()),
                        ("GSE11223", C.MAT_11223, C.groups_gse11223())]:
        if not os.path.exists(mat):
            log(f"[skip] {nm}"); continue
        r = run(nm, mat, gf)
        if r: out[nm] = r
    p = os.path.join(OUT, "partA_step2b_panel_specificity.json")
    json.dump(out, open(p, "w", encoding="utf-8"), ensure_ascii=False, indent=2); log("[write]", p)
    log("[done]")
