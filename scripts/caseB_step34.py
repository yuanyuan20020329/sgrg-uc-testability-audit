# -*- coding: utf-8 -*-
"""案例 B Step 3 + Step 4：复现 Wang 2024 (Sci Rep) 的 PD 焦亡诊断模型，并执行审计阶梯。

数据通道：Gemma REST（NCBI 全域封锁）。论文的 5 个 bulk 队列全部在 Gemma 中：
   训练 GSE7621(549) / GSE20164(3887) / GSE20163(3886) / GSE8397(873)
   验证 GSE49036(11617)
论文口径：diagnostic score = -0.029*POLR2K - 0.018*TIMM8B；
   训练①同集 AUC 目标 0.823；验证 GSE49036 目标 0.706。

审计层：
  L2  每队列单独 AUC vs 合并 AUC（模型是否吃「队列间差异」）
      且逐队列 batchEffect 状态来自 Gemma 策展字段
  L3  诊断评分 vs 全基因平均表达（全局技术轴）的耦合
  L4  随机 2 基因对的同集 AUC 零分布（与面板同协议）
  L5  样本量与方向一致性（设计检验力）
"""
import sys
import os
import re
import json
import time
import random
import collections
import numpy as np
import pandas as pd

sys.path.insert(0, "scripts")
from caseA_gemma import all_samples, _fetch, get_dataset

OUT = "results"
CACHE = "data/raw/caseB/gemma"
os.makedirs(CACHE, exist_ok=True)

DATASETS = {549: "GSE7621", 3887: "GSE20164", 3886: "GSE20163",
            873: "GSE8397", 11617: "GSE49036"}
TRAIN = [549, 3887, 3886, 873]
VAL = [11617]
PANEL = {"POLR2K": 5440, "TIMM8B": 26521}
# 论文的固定权重
W_POLR2K, W_TIMM8B = -0.029, -0.018


# ---------------- 标签 ----------------
def label_of(rec):
    """从 characteristics + 样本名联合判定疾病标签。"""
    ch = collections.defaultdict(list)
    for c in (rec.get("sample") or {}).get("characteristics") or []:
        ch[str(c.get("category"))].append(str(c.get("value")))
    blob = " ".join(" ".join(v) for v in ch.values()).lower()
    nm = (rec.get("name") or "").lower()
    allb = blob + " || " + nm
    if "incidental lewy" in allb or "ilbd" in allb.split():
        return "ILBD"
    if "parkinson" in allb or re.search(r"\bpd\b", allb):
        return "PD"
    if "control" in allb or "normal" in allb or "healthy" in allb:
        return "CTL"
    return None


def build_sample_table():
    recs = []
    for did, nm in DATASETS.items():
        rows = all_samples(did)
        for r in rows:
            recs.append({
                "dataset": nm, "did": did,
                "expr_key": r.get("name"),
                "gsm": r.get("shortName"),
                "label": label_of(r),
                "platform": (r.get("arrayDesign") or {}).get("shortName"),
            })
    df = pd.DataFrame(recs)
    return df


# ---------------- 表达 ----------------
def _fetch_chunk(did, part):
    """取一个基因块，返回 geneExpressionLevels 列表；失败返回 None。"""
    u = "https://gemma.msl.ubc.ca/rest/v2/datasets/%d/expressions/genes/%s?consolidate=average" % (
        did, ",".join(str(x) for x in part))
    try:
        d = json.loads(_fetch(u, timeout=900))
        return (d.get("data") or [{}])[0].get("geneExpressionLevels") or []
    except Exception:
        return None


def fetch_genes(did, entrez_ids, chunk=25):
    """返回 {entrez: {'symbol':..,'expr':{sample:value}}}。带缓存；失败时二分重试。"""
    cache = os.path.join(CACHE, "expr_%d.json" % did)
    if os.path.exists(cache):
        try:
            raw = json.load(open(cache, encoding="utf-8"))
            return {int(k): v for k, v in raw.items()}
        except Exception:
            pass

    res = {}
    ids = [int(x) for x in entrez_ids]

    def absorb(gel):
        for g in gel:
            gid = g.get("geneNcbiId")
            sym = g.get("geneOfficialSymbol")
            vecs = g.get("vectors") or []
            if not vecs:
                continue
            acc = collections.defaultdict(list)
            for v in vecs:
                for k, val in (v.get("bioAssayExpressionLevels") or {}).items():
                    if val is not None:
                        acc[k].append(float(val))
            res[gid] = {"symbol": sym,
                        "expr": {k: float(np.mean(v)) for k, v in acc.items()}}

    def work(part, depth=0):
        if not part:
            return
        gel = _fetch_chunk(did, part)
        if gel is not None:
            absorb(gel)
            time.sleep(0.2)
            return
        if len(part) == 1:
            return                      # 该基因确实取不到，放弃
        mid = len(part) // 2
        if depth < 4:
            work(part[:mid], depth + 1)
            work(part[mid:], depth + 1)
        else:
            print("      skip %d genes (too deep)" % len(part))

    for i in range(0, len(ids), chunk):
        work(ids[i:i + chunk])

    json.dump({str(k): v for k, v in res.items()},
              open(cache, "w", encoding="utf-8"))
    return res


def gene_pool(n=150, seed=20261010):
    """从 MSigDB 缓存取基因宇宙，随机抽 n 个（排除面板基因）。"""
    pool = set()
    for f in ["data/raw/caseA/genesets/h.all.v2023.2.Hs.symbols.gmt",
              "data/raw/caseA/genesets/c2.cp.v2023.1.Hs.symbols.gmt"]:
        if not os.path.exists(f):
            continue
        for line in open(f, encoding="utf-8"):
            parts = line.rstrip("\n").split("\t")
            pool.update(parts[2:])
    pool -= set(PANEL)
    pool = sorted(pool)
    random.seed(seed)
    return random.sample(pool, min(n, len(pool)))


def symbol2entrez(symbols):
    import urllib.parse
    out = {}
    for i in range(0, len(symbols), 400):
        part = symbols[i:i + 400]
        u = "https://mygene.info/v3/query?" + urllib.parse.urlencode({
            "q": "symbol:(" + " OR ".join(part) + ")",
            "species": "human", "fields": "symbol,entrezgene", "size": 1000})
        r = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
        d = json.loads(urllib.request.urlopen(r, timeout=90).read())
        for h in d.get("hits", []):
            if h.get("symbol") and h.get("entrezgene"):
                out[h["symbol"]] = int(h["entrezgene"])
        time.sleep(0.3)
    return out


# ---------------- 主流程 ----------------
def main():
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score

    print("### 1. 样本表 ###")
    df = build_sample_table()
    df.to_csv(os.path.join(OUT, "caseB_samples.tsv"), sep="\t", index=False)
    print(df.groupby(["dataset", "label"]).size().to_string())
    lab_ok = df[df["label"].isin(["PD", "CTL"])].copy()
    print("\n可用样本（PD+CTL）: %d ; PD=%d CTL=%d" % (
        len(lab_ok), (lab_ok.label == "PD").sum(), (lab_ok.label == "CTL").sum()))

    print("\n### 2. 基因启动集 ###")
    pyro_genes = []
    for line in open("data/raw/caseA/genesets/c2.cp.v2023.1.Hs.symbols.gmt", encoding="utf-8"):
        if line.startswith("REACTOME_PYROPTOSIS\t"):
            pyro_genes = line.rstrip("\n").split("\t")[2:]
            break
    print("REACTOME_PYROPTOSIS n =", len(pyro_genes))
    pool_syms = gene_pool(150)
    print("随机池 n =", len(pool_syms))
    s2e = symbol2entrez(sorted(set(pyro_genes) | set(pool_syms)))
    pool_e = [s2e[s] for s in pool_syms if s in s2e]
    pyro_e = [s2e[s] for s in pyro_genes if s in s2e]
    print("映射成功：池 %d / 焦亡 %d" % (len(pool_e), len(pyro_e)))

    all_ids = sorted(set(pool_e) | set(pyro_e) | set(PANEL.values()))

    print("\n### 3. 取表达 ###")
    expr = {}          # symbol -> {expr_key: value}
    per_ds = {}
    for did, nm in DATASETS.items():
        print("   fetching %s (%d genes)..." % (nm, len(all_ids)))
        g = fetch_genes(did, all_ids)
        per_ds[nm] = g
        print("      got %d genes" % len(g))

    # 组装：以 (dataset, expr_key) 为行
    rows = []
    meta = {}
    for nm, g in per_ds.items():
        # 该数据集的样本名集合
        keys = set()
        for gid, v in g.items():
            keys |= set(v["expr"].keys())
        meta[nm] = keys
    # 用标签表的 expr_key 对齐
    for nm in DATASETS.values():
        sub = lab_ok[lab_ok.dataset == nm]
        g = per_ds[nm]
        for _, r in sub.iterrows():
            k = r.expr_key
            rec = {"dataset": nm, "expr_key": k, "label": r.label}
            ok = False
            for gid, v in g.items():
                sym = v["symbol"]
                if sym:
                    rec[sym] = v["expr"].get(k, np.nan)
                    if k in v["expr"]:
                        ok = True
            rec["_has_data"] = ok
            rows.append(rec)
    X = pd.DataFrame(rows)
    X.to_csv(os.path.join(OUT, "caseB_expr_table.tsv"), sep="\t", index=False)
    print("   表达表: %d 行 × %d 列" % X.shape)
    print("   有数据的样本:", int(X["_has_data"].sum()))

    # ---------- Step 3 复现 ----------
    print("\n### 4. Step 3 复现（论文口径）###")
    for c in ["POLR2K", "TIMM8B"]:
        if c not in X.columns:
            print("!! 缺少", c)
    X = X[X["_has_data"]].copy()
    X["score"] = W_POLR2K * X["POLR2K"] + W_TIMM8B * X["TIMM8B"]

    # 全局表达代理（先算，后面 tr/va 才能引用）
    pool_in = [s for s in pool_syms if s in X.columns]
    X["global_mean"] = X[pool_in].mean(axis=1)

    def auc(sub, col="score"):
        s = sub.dropna(subset=[col])
        if s.label.nunique() < 2 or len(s) < 4:
            return None
        return roc_auc_score((s.label == "PD").astype(int), s[col])

    tr = X[X.dataset.isin([DATASETS[d] for d in TRAIN])]
    va = X[X.dataset.isin([DATASETS[d] for d in VAL])]
    res = {}
    res["train_pooled_AUC"] = auc(tr)
    res["val_GSE49036_AUC"] = auc(va)
    print("   训练（4 集合并）同集 AUC = %.4f   [论文 0.823]" % (res["train_pooled_AUC"] or -1))
    print("   验证 GSE49036 AUC        = %.4f   [论文 0.706]" % (res["val_GSE49036_AUC"] or -1))

    # 单基因 AUC（与论文 0.742 / 0.777 对照，判断差异来自数据还是权重）
    print("\n### 4b. 单基因 AUC（诊断差异来源定位）###")
    single = {}
    for g, paper in [("POLR2K", 0.742), ("TIMM8B", 0.777)]:
        a_tr = auc(tr, g)
        a_va = auc(va, g)
        single[g] = {"train": a_tr, "val": a_va, "paper_train": paper}
        print("   %-8s 训练 AUC=%.4f (论文 %.3f) | 验证 AUC=%.4f" % (
            g, a_tr or -1, paper, a_va or -1))
    res["single_gene_AUC"] = single

    # 重拟合权重（同样本集 logistic）——检验「固定权重是否抗尺度」
    t2 = tr.dropna(subset=["POLR2K", "TIMM8B", "label"])
    lg = LogisticRegression(max_iter=2000)
    lg.fit(t2[["POLR2K", "TIMM8B"]].values, (t2.label == "PD").astype(int).values)
    coef = dict(zip(["POLR2K", "TIMM8B"], lg.coef_[0]))
    X["score_refit"] = X["POLR2K"] * coef["POLR2K"] + X["TIMM8B"] * coef["TIMM8B"]
    # 重新取子集（X 已新增列）
    tr = X[X.dataset.isin([DATASETS[d] for d in TRAIN])]
    va = X[X.dataset.isin([DATASETS[d] for d in VAL])]
    res["refit_coef"] = coef
    res["refit_train_AUC"] = auc(tr, "score_refit")
    res["refit_val_AUC"] = auc(va, "score_refit")
    print("   重拟合权重 %s → 训练 AUC=%.4f | 验证 AUC=%.4f" % (
        {k: round(v, 4) for k, v in coef.items()},
        res["refit_train_AUC"] or -1, res["refit_val_AUC"] or -1))
    print("   （论文固定权重 %.3f / %.3f）" % (W_POLR2K, W_TIMM8B))

    # ---------- L2 每队列 vs 合并 ----------
    print("\n### 5. L2 每队列 AUC vs 合并 AUC ###")
    per = {}
    for nm in DATASETS.values():
        sub = X[X.dataset == nm]
        a = auc(sub)
        per[nm] = a
        n_pd = int((sub.label == "PD").sum())
        n_ctl = int((sub.label == "CTL").sum())
        print("   %-9s n=%2d (PD %2d / CTL %2d)  AUC = %s" % (
            nm, len(sub), n_pd, n_ctl, "%.4f" % a if a else "n/a"))
    res["per_dataset_AUC"] = per

    # ---------- L3 全局表达耦合 ----------
    print("\n### 6. L3 诊断评分 vs 全局表达（技术轴）###")
    for nm in list(DATASETS.values()):
        sub = X[X.dataset == nm]
        if len(sub) > 3:
            rho = sub["score"].corr(sub["global_mean"], method="spearman")
            print("   %-9s rho(score, global mean) = %.3f  (n=%d)" % (nm, rho, len(sub)))
    rho_all = X["score"].corr(X["global_mean"], method="spearman")
    print("   全部样本 rho = %.3f" % rho_all)
    # 组内（去队列）耦合：对每队列内中心化后再算
    Xz = X.copy()
    for c in ["score", "global_mean"]:
        Xz[c + "_z"] = Xz.groupby("dataset")[c].transform(
            lambda s: (s - s.mean()) / (s.std() + 1e-9))
    rho_within = Xz["score_z"].corr(Xz["global_mean_z"], method="spearman")
    print("   ★ 队列内中心化后 rho = %.3f  （纯「队列间」耦合 = %.3f）" % (
        rho_within, rho_all))
    res["rho_score_globalmean"] = float(rho_all)
    res["rho_score_globalmean_within_dataset"] = float(rho_within)
    res["global_mean_AUC_train"] = auc(tr, "global_mean")
    res["global_mean_AUC_val"] = auc(va, "global_mean")
    print("   全局均值单独分病：训练 AUC=%.4f | 验证 AUC=%.4f" % (
        res["global_mean_AUC_train"] or -1, res["global_mean_AUC_val"] or -1))

    # ---------- L4 随机 2 基因对 ----------
    print("\n### 7. L4 随机 2 基因对同集 AUC 零分布 ###")
    random.seed(7)
    tr2 = tr.dropna(subset=["POLR2K", "TIMM8B", "label"])
    y_tr = (tr2.label == "PD").astype(int).values
    pool_ok = [s for s in pool_in if tr2[s].notna().all() and tr2[s].std() > 0]
    print("   可用池基因 %d" % len(pool_ok))

    def pair_auc(g1, g2, sub, y):
        Xg = sub[[g1, g2]].values
        clf = LogisticRegression(max_iter=2000)
        clf.fit(Xg, y)
        p = clf.predict_proba(Xg)[:, 1]
        return roc_auc_score(y, p)

    panel_insample = pair_auc("POLR2K", "TIMM8B", tr2, y_tr)
    print("   面板（同协议 logistic 同集）AUC = %.4f" % panel_insample)
    nrand = 500
    rand_aucs = []
    for _ in range(nrand):
        g1, g2 = random.sample(pool_ok, 2)
        try:
            rand_aucs.append(pair_auc(g1, g2, tr2, y_tr))
        except Exception:
            pass
    rand_aucs = np.array(rand_aucs)
    emp_p = (1 + int((rand_aucs >= panel_insample).sum())) / (len(rand_aucs) + 1)
    print("   随机对 AUC: 中位 %.4f  q95 %.4f  max %.4f" % (
        np.median(rand_aucs), np.quantile(rand_aucs, 0.95), rand_aucs.max()))
    print("   面板经验 p（同集）= %.4f" % emp_p)
    res["L4"] = {"panel_insample_AUC": float(panel_insample),
                 "rand_median": float(np.median(rand_aucs)),
                 "rand_q95": float(np.quantile(rand_aucs, 0.95)),
                 "rand_max": float(rand_aucs.max()),
                 "emp_p": float(emp_p), "n_rand": len(rand_aucs)}

    # ---------- L5 设计 ----------
    print("\n### 8. L5 设计/检验力 ###")
    l5 = {
        "train_n": int(len(tr)), "train_pd": int((tr.label == "PD").sum()),
        "val_n": int(len(va)), "val_pd": int((va.label == "PD").sum()),
        "panel_size": 2,
        "gain_train_to_val": (res["train_pooled_AUC"] or 0) - (res["val_GSE49036_AUC"] or 0),
    }
    print("   训练 n=%d (PD %d) → 同集 AUC %.3f；验证 n=%d (PD %d) → AUC %.3f；落差 %.3f" % (
        l5["train_n"], l5["train_pd"], res["train_pooled_AUC"] or 0,
        l5["val_n"], l5["val_pd"], res["val_GSE49036_AUC"] or 0, l5["gain_train_to_val"]))
    res["L5"] = l5

    # ---------- sc 臂设计混杂（论文文本，无需数据） ----------
    print("\n### 9. sc 臂设计混杂（论文文本 + Cramér's V）###")
    cont = pd.DataFrame(
        {"GSE140231": [0, 12], "GSE157783": [5, 6], "GSE204796": [0, 9]},
        index=["PD", "healthy"])
    print(cont.to_string())
    from scipy.stats import chi2_contingency
    chi2, p, dof, exp = chi2_contingency(cont.values)
    n = cont.values.sum()
    v = np.sqrt(chi2 / (n * min(cont.shape[0] - 1, cont.shape[1] - 1)))
    ctl_from_casefree = (12 + 9) / 27
    print("   chi2=%.2f  p=%.4g  Cramér's V=%.3f" % (chi2, p, v))
    print("   对照臂来自「无病例数据集」的比例 = %d/27 = %.3f" % (12 + 9, ctl_from_casefree))
    res["sc_arm_confound"] = {"chi2": float(chi2), "p": float(p), "cramers_v": float(v),
                              "ctl_from_casefree_datasets": ctl_from_casefree}

    res["replication"] = {"train_pooled": res["train_pooled_AUC"],
                          "val": res["val_GSE49036_AUC"]}
    json.dump(res, open(os.path.join(OUT, "caseB_step34_summary.json"), "w",
                        encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
    print("\nwrote results/caseB_step34_summary.json")


if __name__ == "__main__":
    main()
