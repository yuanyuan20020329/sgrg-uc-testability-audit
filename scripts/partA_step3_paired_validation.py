# -*- coding: utf-8 -*-
"""
Part A · Step 3 —— 换用批次平衡队列（GSE22619）重做验证

为什么换：
  GSE38713 的疾病状态与处理批次高度混杂（Cramér's V=0.712，7 批中 4 批纯 UC），
  去批次后 SGRG 面板 AUC 从 0.936 塌到 0.58 —— 单队列证据不支持诊断模型。
  GSE22619 是 UC 不一致型同卵双生子队列（10 对，每对 1 例 UC + 1 例健康同胞），
  且每个处理批次都含完整配对（Batch_01 2/2、Batch_02 4/4、Batch_03 4/4），
  批次与疾病正交。因此：
    (1) 配对差 (UC − 健康同胞) 同时消掉批次效应与遗传背景；
    (2) 该队列可作为独立验证集，检验 GSE38713 训练出的签名是否真实可迁移；
    (3) 配对居中后重新做面板特异性检验，看判别力是否还在。

面板特异性沿用 Step2b 的固定配置（L2, C=0.01, LOO），使新旧数字可直接对比。

产出：
  results/partA_step3_paired_DE.tsv
  results/partA_step3_summary.json
  results/partA_step3_panel_specificity.json
  results/partA_step3_cross_cohort.json
  results/fig_GSE22619_paired_volcano.png
  results/fig_GSE22619_batch_balance.png
  results/fig_step3_crosscohort_ROC.png
"""
import os, sys, json, math, csv, warnings
from collections import Counter
import numpy as np
warnings.filterwarnings("ignore")
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
from sklearn.metrics import roc_auc_score, roc_curve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sgrg_common as C

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = C.OUT
log = C.log
RNG = np.random.default_rng(20261008)
C_FIX, PEN_FIX = 0.01, "l2"          # 与 Step2b 一致，保证可比


# ------------------------------------------------------------------ 载入
def load_cohort(matrix_path, groups_map, batch_map=None, pair_map=None, name=""):
    """返回 SGRG 子矩阵 + 全基因池（供随机面板抽样）"""
    header, rows = C.load_matrix(matrix_path)
    sym_col = C.find_sym_col(header)
    sample_cols = list(range(6, len(header)))
    expr, _ = C.collapse_genes(header, rows, sym_col, sample_cols)
    if len(expr) < 1000:
        raise RuntimeError(f"{name} 矩阵不完整（仅 {len(expr)} 基因）")

    labs = [C.col_keys(header[j]) for j in sample_cols]
    grp = np.array([C.resolve_label(ks, groups_map) for ks in labs], dtype=object)
    bat = np.array([C.resolve_label(ks, batch_map) if batch_map else "?" for ks in labs], dtype=object)
    pr = np.array([C.resolve_label(ks, pair_map) if pair_map else "?" for ks in labs], dtype=object)
    # 防御：若映射值是元组/列表，取首元素，避免被摊成 2D 数组
    for arr in (grp, bat, pr):
        for i, v in enumerate(arr):
            if isinstance(v, (tuple, list)):
                arr[i] = v[0] if len(v) else "?"

    pool = sorted(g for g, v in expr.items() if np.isfinite(v).mean() >= .9)
    Xall = np.vstack([expr[g] for g in pool])          # all_genes x samples
    sgrg = set(C.load_sgrg())
    s_idx = [i for i, g in enumerate(pool) if g in sgrg]
    log(f"[{name}] pool={len(pool)} | SGRG measured={len(s_idx)}/{len(sgrg)} "
        f"| samples={Xall.shape[1]} | groups={dict(Counter(grp))} | batches={dict(Counter(bat))}")
    return dict(name=name, pool=pool, Xall=Xall, s_idx=np.array(s_idx),
                s_genes=[pool[i] for i in s_idx], groups=grp, batch=bat, pair=pr)


def auc_ci(y, s, B=400):
    y = np.asarray(y); s = np.asarray(s); n = len(y); v = []
    for _ in range(B):
        idx = RNG.integers(0, n, n)
        if len(set(y[idx])) < 2:
            continue
        v.append(roc_auc_score(y[idx], s[idx]))
    return (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if v else (float("nan"),) * 2


def metrics(y, s):
    y = np.asarray(y); s = np.asarray(s)
    auc = roc_auc_score(y, s)
    pred = (s >= np.median(s)).astype(int)
    tp = int(((pred == 1) & (y == 1)).sum()); fn = int(((pred == 0) & (y == 1)).sum())
    tn = int(((pred == 0) & (y == 0)).sum()); fp = int(((pred == 1) & (y == 0)).sum())
    lo, hi = auc_ci(y, s)
    return dict(auc=round(float(auc), 4), acc=round((tp + tn) / len(y), 4),
                sens=round(tp / max(tp + fn, 1), 4), spec=round(tn / max(tn + fp, 1), 4),
                tp=tp, fn=fn, tn=tn, fp=fp, auc_lo=round(lo, 4), auc_hi=round(hi, 4))


# ------------------------------------------------------------------ 1) 配对差异表达
def paired_de(co):
    g, pr, X = co["groups"], co["pair"], co["Xall"][co["s_idx"], :]
    genes = co["s_genes"]
    pairs, rows_ = sorted(set(pr) - {"?"}), []
    for p in pairs:
        iu = np.where((pr == p) & (g == "UC"))[0]
        ic = np.where((pr == p) & (g == "Control"))[0]
        if len(iu) == 1 and len(ic) == 1:
            rows_.append(X[:, iu[0]] - X[:, ic[0]])
    D = np.vstack(rows_)
    log(f"[paired] usable twin pairs = {len(pairs)}")
    res = []
    for j, gene in enumerate(genes):
        d = D[:, j]; d = d[np.isfinite(d)]
        if len(d) < 5:
            continue
        try:
            _, pw = stats.wilcoxon(d, alternative="two-sided")
        except Exception:
            pw = float("nan")
        _, pt = stats.ttest_1samp(d, 0.0)
        res.append(dict(gene=gene, n_pairs=len(d), mean_diff=float(np.mean(d)),
                        median_diff=float(np.median(d)),
                        cohen_dz=float(np.mean(d) / (np.std(d, ddof=1) + 1e-12)),
                        wilcoxon_p=float(pw), ttest_p=float(pt)))
    q = C.bh_fdr([r["wilcoxon_p"] for r in res])
    for r, qq in zip(res, q):
        r["wilcoxon_FDR"] = float(qq)
    res.sort(key=lambda r: r["wilcoxon_p"])
    return res, len(pairs)


# ------------------------------------------------------------------ 2) 批次诊断
def batch_diag(co, y):
    """批次×疾病 的混杂诊断。y 为 0/1 疾病标签（与样本列对齐），
    因此对不同队列的组名体系（UC / UC_active_inflamed）都成立。"""
    b = co["batch"]
    y = np.asarray(y).astype(int)
    cats = sorted(set(b))
    tab = np.array([[int(((b == c) & (y == 1)).sum()), int(((b == c) & (y == 0)).sum())]
                    for c in cats])
    n = int(len(y))
    try:
        chi2, p, dof, _ = stats.chi2_contingency(tab)
    except Exception:
        chi2, p, dof = float("nan"), float("nan"), 0
    k = min(tab.shape) - 1
    v = math.sqrt(chi2 / (n * k)) if (k > 0 and np.isfinite(chi2) and n > 0) else float("nan")
    Xb = np.zeros((n, len(cats)))
    for i, c in enumerate(cats):
        Xb[:, i] = (b == c).astype(float)
    oof = np.zeros(n)
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=7).split(Xb, y):
        if len(set(y[tr])) < 2:
            oof[te] = float(np.mean(y[tr])); continue
        m = Pipeline([("sc", StandardScaler()), ("lr", LogisticRegression(max_iter=2000))])
        m.fit(Xb[tr], y[tr]); oof[te] = m.predict_proba(Xb[te])[:, 1]
    auc_b = float(roc_auc_score(y, oof)) if len(set(y)) > 1 else float("nan")
    log(f"[batch:{co['name']}] Cramér's V={v:.4f} (chi2={chi2:.2f}, p={p:.3g}) | batch-only AUC={auc_b:.3f}")
    return dict(batch_cats={c: dict(Case=int(t[0]), Control=int(t[1])) for c, t in zip(cats, tab)},
                cramers_v=round(float(v), 4), chi2=round(float(chi2), 3),
                chi2_p=float(p), chi2_dof=int(dof),
                batch_only_auc=round(auc_b, 4), n=n, n_case=int(y.sum()))


# ------------------------------------------------------------------ 3) 面板特异性
def panel_specificity(Xall, pool, s_idx, y, label="", n_rand=200):
    def cv_auc(M):                                # M: samples x genes
        n = M.shape[0]
        cv = LeaveOneOut() if n <= 60 else StratifiedKFold(5, shuffle=True, random_state=3)
        oof = np.full(n, np.nan)
        for tr, te in cv.split(M, y):
            if len(set(y[tr])) < 2:
                continue
            m = Pipeline([("sc", StandardScaler()),
                          ("lr", LogisticRegression(C=C_FIX, penalty=PEN_FIX,
                                                    solver="liblinear", max_iter=5000))])
            m.fit(M[tr], y[tr]); oof[te] = m.predict_proba(M[te])[:, 1]
        if np.isnan(oof).any():
            oof[np.isnan(oof)] = 0.5
        return roc_auc_score(y, oof)

    Xs = Xall[s_idx, :].T
    obs = cv_auc(Xs)
    k = len(s_idx)
    others = np.array([i for i in range(len(pool)) if i not in set(s_idx.tolist())])
    rnd = np.array([cv_auc(Xall[RNG.choice(others, k, replace=False), :].T) for _ in range(n_rand)])
    pct = float((rnd < obs).mean() * 100)
    log(f"[panel:{label}] observed={obs:.3f} | random median={np.median(rnd):.3f} "
        f"p95={np.percentile(rnd,95):.3f} | percentile={pct:.1f} | "
        f"emp_p={(rnd>=obs).mean():.4f} (k={k}, pool={len(others)})")
    return dict(label=label, k=k, pool=int(len(others)), n_random=n_rand,
                observed=round(float(obs), 4), rand_median=round(float(np.median(rnd)), 4),
                rand_p95=round(float(np.percentile(rnd, 95)), 4),
                percentile=round(pct, 1), empirical_p=round(float((rnd >= obs).mean()), 4),
                randoms=[round(float(x), 4) for x in rnd])


# ------------------------------------------------------------------ 3b) 跨队列面板特异性
def _transfer(Xtr, ytr, Xte, yte, rows, Cv=C_FIX, pen=PEN_FIX):
    Xs = Xtr[rows, :].T; Xt = Xte[rows, :].T
    m = Pipeline([("sc", StandardScaler()),
                  ("lr", LogisticRegression(C=Cv, penalty=pen, solver="liblinear", max_iter=5000))])
    m.fit(Xs, ytr)
    s = m.predict_proba(Xt)[:, 1]
    return roc_auc_score(yte, s), s


def cross_panel_specificity(co_tr, y_tr, co_te, y_te, sgrg_set, label="", n_rand=200,
                            Cv=C_FIX, pen=PEN_FIX):
    """在训练队列上拟合「同规模随机基因面板」，**用同一模型配置**迁移到独立队列测 AUC。
    这是本场景下唯一有判别力的特异性检验：
      n=20 / p=82 时同队列内 CV 会让几乎所有基因集都逼近 AUC≈1（本次实测随机面板中位 0.955），
      因此「同队列 CV AUC」不能用来论证 SGRG 的特异性；跨队列迁移才可以。
    注意：随机面板必须与 SGRG 面板使用同一个 (penalty, C)，否则比较不公平。
    """
    tr_pool, te_pool = co_tr["pool"], co_te["pool"]
    common_pool = sorted(set(tr_pool) & set(te_pool))
    i_tr = {g: i for i, g in enumerate(tr_pool)}
    i_te = {g: i for i, g in enumerate(te_pool)}
    Xtr = co_tr["Xall"][[i_tr[g] for g in common_pool], :]
    Xte = co_te["Xall"][[i_te[g] for g in common_pool], :]
    s_rows = np.array([i for i, g in enumerate(common_pool) if g in sgrg_set])
    others = np.array([i for i, g in enumerate(common_pool) if g not in sgrg_set])

    auc_obs, s_obs = _transfer(Xtr, y_tr, Xte, y_te, s_rows, Cv=Cv, pen=pen)
    rnd = np.array([_transfer(Xtr, y_tr, Xte, y_te,
                              RNG.choice(others, len(s_rows), replace=False), Cv=Cv, pen=pen)[0]
                    for _ in range(n_rand)])
    pct = float((rnd < auc_obs).mean() * 100)
    log(f"[cross-panel:{label}] SGRG 面板迁移 AUC={auc_obs:.3f} | 随机面板中位={np.median(rnd):.3f} "
        f"p95={np.percentile(rnd,95):.3f} | {pct:.0f} 百分位 | 经验 p={(rnd>=auc_obs).mean():.4f} "
        f"(k={len(s_rows)}, pool={len(others)}, {pen.upper()} C={Cv})")
    return dict(label=label, k=int(len(s_rows)), pool=int(len(others)), n_random=n_rand,
                model=f"{pen.upper()} C={Cv}",
                sgrg_transfer_auc=round(float(auc_obs), 4),
                rand_median=round(float(np.median(rnd)), 4),
                rand_p95=round(float(np.percentile(rnd, 95)), 4),
                percentile=round(pct, 1),
                empirical_p=round(float((rnd >= auc_obs).mean()), 4),
                randoms=[round(float(x), 4) for x in rnd])


# ------------------------------------------------------------------ 4) 跨队列模型
def fit_train_eval(Xtr, ytr, Xte, yte, grid=(0.01, 0.03, 0.1, 0.3, 1.0)):
    best = None
    for pen in ("l1", "l2"):
        for Cv in grid:
            n = len(ytr); oof = np.zeros(n)
            for tr, te in StratifiedKFold(5, shuffle=True, random_state=13).split(Xtr.T, ytr):
                m = Pipeline([("sc", StandardScaler()),
                              ("lr", LogisticRegression(C=Cv, penalty=pen,
                                                        solver="liblinear", max_iter=5000))])
                m.fit(Xtr[:, tr].T, ytr[tr]); oof[te] = m.predict_proba(Xtr[:, te].T)[:, 1]
            a = roc_auc_score(ytr, oof)
            if best is None or a > best[0]:
                best = (a, Cv, pen)
    cv, Cv, pen = best
    m = Pipeline([("sc", StandardScaler()),
                  ("lr", LogisticRegression(C=Cv, penalty=pen, solver="liblinear", max_iter=5000))])
    m.fit(Xtr.T, ytr)
    s = m.predict_proba(Xte.T)[:, 1]
    coef = m.named_steps["lr"].coef_.ravel() * m.named_steps["sc"].scale_
    return dict(cv_auc=round(float(cv), 4), C=Cv, penalty=pen,
                n_nonzero=int((np.abs(coef) > 1e-8).sum()),
                ext=metrics(yte, s), scores=[round(float(x), 5) for x in s],
                coef=[round(float(x), 5) for x in coef])


# ------------------------------------------------------------------ main
def main():
    # ---- 选矩阵：优先完整文件；缺失时允许使用「按行截断的前缀」并显式标注
    matrix_status = "FULL"
    if os.path.exists(C.MAT_22619):
        mat226 = C.MAT_22619
    elif os.path.exists(C.MAT_22619_PARTIAL):
        mat226 = C.MAT_22619_PARTIAL
        matrix_status = "TRUNCATED_PREFIX"
        log("!! 使用 GSE22619 的部分矩阵（下载被截断，仅保留按行前缀）—— 结论必须带此限定")
    else:
        log("!! GSE22619 矩阵不存在，退出"); return

    co387 = load_cohort(C.MAT_38713, C.groups_gse38713(), C.batch_gse38713(), name="GSE38713")
    co226 = load_cohort(mat226, C.groups_gse22619(), C.batch_gse22619(),
                        C.pairs_gse22619(), name="GSE22619")
    if matrix_status != "FULL":
        log(f"[{matrix_status}] 该队列 SGRG 可用 {len(co226['s_genes'])}/104（前缀非随机，缺失基因集中在阵列后段探针）")

    g387 = co387["groups"]
    y226 = (co226["groups"] == "UC").astype(int)
    y387_act = np.isin(g387, ["UC_active_inflamed", "UC_active_noninflamed", "UC_remission"]).astype(int)
    y387_inf = (g387 == "UC_active_inflamed").astype(int)

    # 公共 SGRG（跨队列可比）
    common = sorted(set(co387["s_genes"]) & set(co226["s_genes"]))
    i387 = [co387["s_genes"].index(g) for g in common]
    i226 = [co226["s_genes"].index(g) for g in common]
    X387 = co387["Xall"][co387["s_idx"][i387], :]
    X226 = co226["Xall"][co226["s_idx"][i226], :]
    log(f"[common] SGRG measured in both cohorts: {len(common)}")

    # 1) 配对 DE
    paired, n_pairs = paired_de(co226)
    pmap = {r["gene"]: r for r in paired}

    de387 = {}
    p387 = os.path.join(OUT, "partA_GSE38713_SGRG_DE.tsv")
    if os.path.exists(p387):
        for r in csv.DictReader(open(p387, encoding="utf-8"), delimiter="\t"):
            de387[r["GeneSymbol"]] = dict(log2FC=float(r["log2FC_inflamed_vs_ctrl"]),
                                          FDR=float(r["FDR_inflamed"]))
    cand = [g for g, v in de387.items() if v["FDR"] < 0.05 and g in pmap]
    n_same = sum(1 for g in cand if np.sign(de387[g]["log2FC"]) == np.sign(pmap[g]["mean_diff"]))
    n_cand = len(cand)
    bt = stats.binomtest(n_same, n_cand, 0.5, alternative="greater") if n_cand else None
    log(f"[replication] GSE38713 显著基因在 GSE22619 配对中方向一致: {n_same}/{n_cand} "
        f"= {100*n_same/max(n_cand,1):.0f}% | 二项 p={(bt.pvalue if bt else float('nan')):.4g}")

    # 全部可比基因的整体方向一致性（图 5 与 summary 共用）
    _pts = [(g, de387[g]["log2FC"], pmap[g]["mean_diff"], pmap[g]["wilcoxon_p"])
            for g in pmap if g in de387]
    _xs = np.array([p[1] for p in _pts]); _ys = np.array([p[2] for p in _pts])
    conc_n_same = int((np.sign(_xs) == np.sign(_ys)).sum()); conc_n = len(_pts)
    conc_rho, conc_rp = stats.spearmanr(_xs, _ys)
    log(f"[concordance] 全部可比基因 {conc_n_same}/{conc_n} = {100*conc_n_same/conc_n:.0f}% 同向 "
        f"| Spearman ρ={conc_rho:.3f} p={conc_rp:.2g}")

    # 2) 批次诊断
    diag226 = batch_diag(co226, y226)
    diag387 = batch_diag(co387, y387_act)

    # 3) 面板特异性（未配对 / 配对居中）
    m226 = np.isin(co226["groups"], ["UC", "Control"])
    Xall226 = co226["Xall"][:, m226]
    spec_raw = panel_specificity(Xall226, co226["pool"],
                                 np.array([i for i, g in enumerate(co226["pool"])
                                           if g in co226["s_genes"]]),
                                 (co226["groups"][m226] == "UC").astype(int),
                                 label="GSE22619 原始表达")

    # 配对居中：每对减去对均值 → 消批次 + 消个体
    pr = co226["pair"][m226]
    Xc = Xall226.astype(float).copy()
    for p in sorted(set(pr) - {"?"}):
        idx = np.where(pr == p)[0]
        if len(idx) == 2:
            Xc[:, idx] -= Xc[:, idx].mean(axis=1, keepdims=True)
    # 居中后「全基因均值」也归零 → 必须改成按基因 z 标准化后同口径评估（cv_auc 内部会标准化）
    spec_cent = panel_specificity(Xc, co226["pool"],
                                  np.array([i for i, g in enumerate(co226["pool"])
                                            if g in co226["s_genes"]]),
                                  (co226["groups"][m226] == "UC").astype(int),
                                  label="GSE22619 配对居中")

    # 3b) 跨队列面板特异性（唯一有判别力的口径）
    sg_set = set(C.load_sgrg())
    m387_inf = np.isin(g387, ["UC_active_inflamed", "Control"])
    tr387 = dict(pool=co387["pool"], Xall=co387["Xall"][:, m387_inf])
    te226 = dict(pool=co226["pool"], Xall=co226["Xall"])
    tr226 = dict(pool=co226["pool"], Xall=co226["Xall"])
    te387 = dict(pool=co387["pool"], Xall=co387["Xall"][:, m387_inf])
    cross_spec = {
        "fixed_L2C0.01 · 38713(inflamed)→22619": cross_panel_specificity(
            tr387, y387_inf[m387_inf].astype(int), te226, y226.astype(int), sg_set,
            label="38713 inflamed → 22619 (固定配置)", n_rand=200, Cv=0.01, pen="l2"),
        "fixed_L2C0.01 · 22619→38713(inflamed)": cross_panel_specificity(
            tr226, y226.astype(int), te387, y387_inf[m387_inf].astype(int), sg_set,
            label="22619 → 38713 inflamed (固定配置)", n_rand=200, Cv=0.01, pen="l2"),
    }

    # 配对检验的可达下限：Wilcoxon n=10 的最小双侧 p = 2/2^10；BH 对 m 个检验的最小阈值 = 0.05/m
    m_tests = int(sum(1 for r in paired if np.isfinite(r["wilcoxon_p"])))
    min_p = 2 / 2 ** n_pairs
    bh_min = 0.05 / max(m_tests, 1)
    fdr_reachable = bool(min_p <= bh_min)
    log(f"[power] Wilcoxon 最小可达 p={min_p:.5f} ({n_pairs} 对) | BH 最小阈值=0.05/{m_tests}={bh_min:.5f} "
        f"→ FDR<0.05 {'可达' if fdr_reachable else '在数学上不可达'}")

    # 4) 跨队列模型
    cross = {}
    cross["38713_allUC→22619"] = fit_train_eval(X387, y387_act, X226, y226)
    r = cross["38713_allUC→22619"]
    log(f"[cross] 38713(allUC)→22619: CV={r['cv_auc']:.3f} ({r['penalty'].upper()},C={r['C']}) "
        f"| 外部 AUC={r['ext']['auc']:.3f} [{r['ext']['auc_lo']},{r['ext']['auc_hi']}] "
        f"| 敏={r['ext']['sens']:.3f} 特={r['ext']['spec']:.3f} | 非零={r['n_nonzero']}")

    cross["38713_inflamed→22619"] = fit_train_eval(X387[:, np.isin(g387, ["UC_active_inflamed", "Control"])],
                                                   (g387[np.isin(g387, ["UC_active_inflamed", "Control"])]
                                                    == "UC_active_inflamed").astype(int), X226, y226)
    r = cross["38713_inflamed→22619"]
    log(f"[cross] 38713(炎症区)→22619: CV={r['cv_auc']:.3f} | 外部 AUC={r['ext']['auc']:.3f}")

    cross["22619→38713_inflamed"] = fit_train_eval(X226, y226, X387, y387_inf)
    r = cross["22619→38713_inflamed"]
    log(f"[cross] 22619→38713(炎症区): CV={r['cv_auc']:.3f} | 外部 AUC={r['ext']['auc']:.3f}")

    # 签名（GSE38713 训练所得，看有多少落进配对验证的方向）
    sig = [(common[i], c) for i, c in enumerate(cross["38713_allUC→22619"]["coef"])
           if abs(c) > 1e-8]
    log("[signature 38713→22619] " + ", ".join(f"{g}({c:+.2f})" for g, c in sig[:10]))

    # 稀疏签名（L1 最优模型，仅保留非零基因）
    best = cross["38713_inflamed→22619"]
    sparse = [(common[i], c) for i, c in enumerate(best["coef"]) if abs(c) > 1e-8]
    sparse.sort(key=lambda t: -abs(t[1]))
    log(f"[sparse signature L1 C={best['C']}] {len(sparse)} genes: "
        + ", ".join(f"{g}({c:+.3f})" for g, c in sparse))

    # 3c) 公平对照：随机面板用与最优模型**相同**的配置做跨队列迁移
    pen_best = best["penalty"]; C_best = best["C"]
    cross_spec[f"tuned_{pen_best.upper()}C{C_best} · 38713(inflamed)→22619"] = \
        cross_panel_specificity(tr387, y387_inf[m387_inf].astype(int), te226, y226.astype(int),
                                sg_set, label=f"38713 inflamed → 22619 (同最优配置 {pen_best.upper()} C={C_best})",
                                n_rand=200, Cv=C_best, pen=pen_best)

    # ---- 落盘
    with open(os.path.join(OUT, "partA_step3_paired_DE.tsv"), "w", encoding="utf-8") as f:
        f.write("\t".join(["GeneSymbol", "n_pairs", "mean_diff_UC_minus_twin", "median_diff",
                           "cohen_dz", "wilcoxon_p", "wilcoxon_FDR", "ttest_p",
                           "GSE38713_log2FC", "GSE38713_FDR", "direction_match"]) + "\n")
        for r in paired:
            a = de387.get(r["gene"])
            dm = "" if a is None else ("+" if np.sign(a["log2FC"]) == np.sign(r["mean_diff"]) else "-")
            f.write("\t".join([r["gene"], str(r["n_pairs"]), f'{r["mean_diff"]:.4f}',
                               f'{r["median_diff"]:.4f}', f'{r["cohen_dz"]:.4f}',
                               f'{r["wilcoxon_p"]:.4e}', f'{r["wilcoxon_FDR"]:.4e}',
                               f'{r["ttest_p"]:.4e}',
                               "" if a is None else f'{a["log2FC"]:.4f}',
                               "" if a is None else f'{a["FDR"]:.4e}', dm]) + "\n")

    summary = dict(
        matrix_status_of_GSE22619=matrix_status,
        cohort_22619=dict(n_samples=int(Xall226.shape[1]), n_uc=int(y226.sum()),
                          n_control=int((y226 == 0).sum()), n_pairs=n_pairs,
                          sgrg_measured=len(co226["s_genes"]),
                          gene_pool=len(co226["pool"])),
        cohort_38713=dict(n_samples=int(co387["Xall"].shape[1]),
                          sgrg_measured=len(co387["s_genes"]),
                          n_uc_active=int(y387_act.sum()), n_control=int((g387 == "Control").sum())),
        common_genes=len(common),
        paired_DE=dict(n_FDR05=int(sum(1 for r in paired if r["wilcoxon_FDR"] < 0.05)),
                       n_nominal_p05=int(sum(1 for r in paired if r["wilcoxon_p"] < 0.05)),
                       n_up=int(sum(1 for r in paired if r["wilcoxon_FDR"] < 0.05 and r["mean_diff"] > 0)),
                       n_down=int(sum(1 for r in paired if r["wilcoxon_FDR"] < 0.05 and r["mean_diff"] < 0)),
                       min_possible_p=round(min_p, 6),
                       bh_min_threshold=round(bh_min, 6),
                       fdr_reachable=fdr_reachable,
                       n_tests=m_tests,
                       n_nominal_p05_and_dz_above_0_3=int(sum(
                           1 for r in paired if r["wilcoxon_p"] < 0.05 and abs(r["cohen_dz"]) >= 0.3)),
                       top=[dict(gene=r["gene"], dz=round(r["cohen_dz"], 3),
                                 mean_diff=round(r["mean_diff"], 3),
                                 wilcoxon_p=r["wilcoxon_p"], wilcoxon_FDR=r["wilcoxon_FDR"])
                            for r in paired[:20]]),
        replication=dict(n_GSE38713_sig=n_cand, n_same_direction=n_same,
                         pct=round(100 * n_same / max(n_cand, 1), 1),
                         binom_p=(float(bt.pvalue) if bt else None)),
        concordance=dict(n_same=conc_n_same, n_total=conc_n,
                         pct=round(100 * conc_n_same / max(conc_n, 1), 1),
                         rho=round(float(conc_rho), 4), spearman_p=float(conc_rp)),
        batch_diag=dict(GSE22619=diag226, GSE38713=diag387),
        panel_specificity=dict(raw=spec_raw, pair_centered=spec_cent,
                               caveat="同队列 CV 在 n=20/p=82 下所有基因集都逼近 AUC≈1，无判别力"),
        cross_cohort_panel_specificity=cross_spec,
        cross_cohort={k: dict(cv_auc=v["cv_auc"], penalty=v["penalty"], C=v["C"],
                              n_nonzero=v["n_nonzero"], ext=v["ext"]) for k, v in cross.items()},
        signature_38713_to_22619=[dict(gene=g, coef=c) for g, c in sig],
        sparse_signature=dict(model=f"{pen_best.upper()} C={C_best}", n_genes=len(sparse),
                              genes=[dict(gene=g, coef=c) for g, c in sparse]),
    )
    json.dump(summary, open(os.path.join(OUT, "partA_step3_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    json.dump(dict(raw=spec_raw, pair_centered=spec_cent, cross_cohort=cross_spec),
              open(os.path.join(OUT, "partA_step3_panel_specificity.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    json.dump({k: dict(cv_auc=v["cv_auc"], penalty=v["penalty"], C=v["C"],
                       n_nonzero=v["n_nonzero"], ext=v["ext"], scores=v["scores"],
                       genes=common, coef=v["coef"]) for k, v in cross.items()},
              open(os.path.join(OUT, "partA_step3_cross_cohort.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)

    # ---- 图 1：配对火山
    fig, ax = plt.subplots(figsize=(6.6, 5.0), dpi=170)
    md = np.array([r["mean_diff"] for r in paired])
    lp = -np.log10(np.clip([r["wilcoxon_p"] for r in paired], 1e-300, 1))
    fq = np.array([r["wilcoxon_FDR"] for r in paired])
    ax.scatter(md[fq >= 0.05], lp[fq >= 0.05], s=16, c="#9fb3c8", alpha=.85, label="FDR ≥ 0.05")
    up = (fq < 0.05) & (md > 0); dn = (fq < 0.05) & (md < 0)
    if up.any(): ax.scatter(md[up], lp[up], s=28, c="#c0392b", label=f"FDR<0.05 ↑ ({int(up.sum())})")
    if dn.any(): ax.scatter(md[dn], lp[dn], s=28, c="#1e6f9f", label=f"FDR<0.05 ↓ ({int(dn.sum())})")
    for r in paired[:8]:
        ax.annotate(r["gene"], (r["mean_diff"], -math.log10(max(r["wilcoxon_p"], 1e-300))),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.axvline(0, color="#888", lw=.8, ls="--")
    ax.set_xlabel("within-pair difference  (UC − healthy co-twin)")
    ax.set_ylabel("−log10 Wilcoxon p")
    ax.set_title(f"GSE22619 · paired SGRG differential expression ({n_pairs} twin pairs)\n"
                 "pairing removes both batch effect and genetic background", fontsize=9)
    ax.legend(fontsize=7); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_GSE22619_paired_volcano.png")); plt.close(fig)

    # ---- 图 2：批次平衡对照
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), dpi=170)
    for ax, (lbl, dg) in zip(axes, [("GSE22619  batch-balanced", diag226),
                                    ("GSE38713  batch-confounded", diag387)]):
        cats = list(dg["batch_cats"].keys())
        uc = [dg["batch_cats"][c]["Case"] for c in cats]
        ct = [dg["batch_cats"][c]["Control"] for c in cats]
        x = np.arange(len(cats)); w = .38
        ax.bar(x - w / 2, uc, w, color="#c0392b", label="UC")
        ax.bar(x + w / 2, ct, w, color="#7f9db9", label="Control")
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace("Batch_", "B").split("_")[0] for c in cats], fontsize=8)
        ax.set_ylabel("n samples")
        ax.set_title(f"{lbl}\nCramér's V={dg['cramers_v']:.3f} · batch-only AUC={dg['batch_only_auc']:.3f}",
                     fontsize=9)
        ax.legend(fontsize=7); ax.grid(alpha=.25, axis="y")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_GSE22619_batch_balance.png")); plt.close(fig)

    # ---- 图 3：跨队列 ROC
    fig, ax = plt.subplots(figsize=(5.8, 5.0), dpi=170)
    for lbl, key, truth, c in [("38713(all-UC) → 22619", "38713_allUC→22619", y226, "#c0392b"),
                               ("38713(inflamed) → 22619", "38713_inflamed→22619", y226, "#e08a5b"),
                               ("22619 → 38713(inflamed)", "22619→38713_inflamed", y387_inf, "#1e6f9f")]:
        s = np.array(cross[key]["scores"])
        fpr, tpr, _ = roc_curve(truth, s)
        ax.plot(fpr, tpr, color=c, lw=1.9, label=f"{lbl}   AUC={cross[key]['ext']['auc']:.3f}")
    ax.plot([0, 1], [0, 1], "--", color="#999", lw=.9)
    ax.set_xlabel("1 − specificity"); ax.set_ylabel("sensitivity")
    ax.set_title("Cross-cohort external validation\n(trained on one cohort, tested on the other)", fontsize=9)
    ax.legend(fontsize=7, loc="lower right"); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_step3_crosscohort_ROC.png")); plt.close(fig)

    # ---- 图 4：跨队列面板特异性
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.1), dpi=170)
    for ax, (k, v) in zip(axes, cross_spec.items()):
        r = np.array(v["randoms"])
        ax.hist(r, bins=26, color="#9db4c9", edgecolor="#fff", lw=.5)
        ax.axvline(v["sgrg_transfer_auc"], color="#b42318", lw=2.0,
                   label=f"SGRG 面板 = {v['sgrg_transfer_auc']:.3f}")
        ax.axvline(float(np.median(r)), color="#5b6470", ls="--", lw=1.0,
                   label=f"随机面板中位 = {np.median(r):.3f}")
        ax.set_title(f"{k}\n{v['percentile']:.0f} 百分位 · 经验 p={v['empirical_p']:.3f}", fontsize=9)
        ax.set_xlabel("在独立队列上的 AUC"); ax.set_ylabel("count")
        ax.legend(fontsize=7, frameon=False); ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_step3_cross_panel_specificity.png")); plt.close(fig)

    # ---- 图 5：方向一致性散点（GSE38713 log2FC vs GSE22619 配对差）
    pts = _pts
    xs = _xs; ys = _ys
    ps = np.array([p[3] for p in pts])
    n_s, n_t = conc_n_same, conc_n
    rho, rp = conc_rho, conc_rp
    fig, ax = plt.subplots(figsize=(6.4, 5.6), dpi=170)
    ax.axhline(0, color="#bbb", lw=.8); ax.axvline(0, color="#bbb", lw=.8)
    ax.plot([-2, 2], [-2, 2], ls=":", color="#999", lw=.9)
    sig = ps < 0.05
    ax.scatter(xs[~sig], ys[~sig], s=22, c="#9fb3c8",
               label=f"配对检验 p ≥ 0.05  (n={int((~sig).sum())})")
    ax.scatter(xs[sig], ys[sig], s=34, c="#c0392b",
               label=f"配对检验 p < 0.05  (n={int(sig.sum())})")
    for g, x, y, p in pts:
        if p < 0.02 or abs(x) > 0.6 or abs(y) > 0.4:
            ax.annotate(g, (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("GSE38713  log2FC  (inflamed UC vs control, unpaired)")
    ax.set_ylabel("GSE22619  within-pair difference  (UC − healthy co-twin)")
    ax.set_title(f"方向一致性：{n_s}/{n_t} = {100*n_s/n_t:.0f}%   Spearman ρ={rho:.3f} (p={rp:.2g})\n"
                 f"同一批 SGRG 基因在两个队列中的位移方向", fontsize=9)
    ax.legend(fontsize=7); ax.grid(alpha=.22)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_step3_direction_concordance.png")); plt.close(fig)
    log(f"[concordance] 全部可比基因 {n_s}/{n_t} = {100*n_s/n_t:.0f}% 同向 | Spearman ρ={rho:.3f} p={rp:.2g}")

    log("\n=== Step3 完成 ===")
    log(f"配对DE: FDR<0.05 {summary['paired_DE']['n_FDR05']} (nominal p<0.05 {summary['paired_DE']['n_nominal_p05']})"
        f" | FDR 可达={fdr_reachable} (最小 p={min_p:.5f} vs BH 阈值 {bh_min:.5f})")
    log(f"方向复现: {n_same}/{n_cand} ({summary['replication']['pct']}%) p={summary['replication']['binom_p']:.3g}")
    log(f"批次: 22619 V={diag226['cramers_v']:.3f}/AUC={diag226['batch_only_auc']:.3f}  |  "
        f"38713 V={diag387['cramers_v']:.3f}/AUC={diag387['batch_only_auc']:.3f}")
    log(f"同队列面板(无判别力): 原始 {spec_raw['observed']:.3f} (随机中位 {spec_raw['rand_median']:.3f}) "
        f"| 配对居中 {spec_cent['observed']:.3f} (随机中位 {spec_cent['rand_median']:.3f})")
    for k, v in cross_spec.items():
        log(f"跨队列面板特异性 [{k}]: SGRG={v['sgrg_transfer_auc']:.3f} "
            f"vs 随机中位 {v['rand_median']:.3f} ({v['percentile']:.0f}百分位, p={v['empirical_p']:.3f})")


if __name__ == "__main__":
    main()
