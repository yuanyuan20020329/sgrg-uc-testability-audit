# -*- coding: utf-8 -*-
"""
案例 A · 在**原文自己的队列 GSE65372** 上做三口径对比（伪重复量化）—— 第二轮。

数据：GSE65372（GPL14951 / Illumina HumanHT-12 v3），54 样本，18 患者，
      39 tumor(Y) / 15 non-tumor(N)，14 患者 Y&N 俱全。
注释：illuminaHumanv3.db 的 ExtraInfo.SymbolReannotated（49,576 探针）。

三口径（对同一批铁死亡基因）：
  P1 原文口径：54 样本，独立样本 Student's t（把患者内重复当独立）
  P2 配对口径：仅 14 个同患者 Y&N；tumor 侧先按患者取均值 -> 配对 t（n=14）
  P3 患者级聚合：每位患者 tumor 均值 vs normal 均值（同样 14 对）—— 与 P2 等价，用作核对
再加：
  P4 患者内平均后再做 AUC（同集 logistic），与 P1 的 AUC 对比

输出 results/caseA_GSE65372_paired.json / .md
"""
import gzip, os, json, sqlite3, collections, math
import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

ROOT = r"C:\Users\H-P\WorkBuddy\应激颗粒"
SERIES = os.path.join(ROOT, "data", "raw", "caseA", "GSE65372_series_matrix.txt.gz")
SQLITE = os.path.join(ROOT, "data", "raw", "caseA", "annot", "illuminaHumanv3",
                      "illuminaHumanv3.db", "inst", "extdata", "illuminaHumanv3.sqlite")
FGENES = os.path.join(ROOT, "results", "caseA_ferroptosis_genes.txt")
RES = os.path.join(ROOT, "results")


def load_series():
    txt = gzip.decompress(open(SERIES, "rb").read()).decode("utf-8", "replace")
    lines = [l.rstrip("\n") for l in txt.split("\n")]

    def getline(pfx):
        for l in lines:
            if l.startswith(pfx + "\t"):
                return [x.strip('"') for x in l.split("\t")[1:]]
        return None

    charlines = [l for l in lines if l.startswith("!Sample_characteristics_ch1")]
    gsm = getline("!Sample_geo_accession")
    pids = [x.strip('"').split(": ")[-1] for x in charlines[0].split("\t")[1:]]
    tns = [x.strip('"').split(": ")[-1] for x in charlines[1].split("\t")[1:]]
    start = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
    samples = [x.strip('"') for x in lines[start + 1].split("\t")[1:]]
    probe_rows = {}
    for l in lines[start + 2:]:
        if l.startswith("!series_matrix_table_end"):
            break
        p = l.split("\t")
        if len(p) < 2:
            continue
        probe_rows[p[0].strip('"')] = np.array(
            [float(x) if x not in ("", "NA", "null") else np.nan for x in p[1:]], float)
    return gsm, pids, tns, samples, probe_rows


def load_symbol_map():
    con = sqlite3.connect(SQLITE)
    cur = con.cursor()
    cur.execute("SELECT IlluminaID, SymbolReannotated, EntrezReannotated FROM ExtraInfo")
    m = {}
    for pid, sym, eg in cur.fetchall():
        if sym:
            m.setdefault(pid, []).append(str(sym))
    con.close()
    return m


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def main():
    gsm, pids, tns, samples, probe_rows = load_series()
    symmap = load_symbol_map()
    fgenes = [l.strip() for l in open(FGENES, encoding="utf-8") if l.strip()]
    fset = set(fgenes)
    print("samples=%d  patients=%d  probes=%d  ferroptosis genes=%d" % (
        len(samples), len(set(pids)), len(probe_rows), len(fgenes)))

    # 探针 -> 铁死亡基因（取每个基因的最好探针=方差最大）
    gene_probes = collections.defaultdict(list)
    for pid, syms in symmap.items():
        if pid not in probe_rows:
            continue
        for s in syms:
            if s in fset:
                gene_probes[s].append(pid)
    print("mapped ferroptosis genes with >=1 probe:", len(gene_probes))

    # 每个基因选方差最大的探针
    gene_vec = {}
    for g, ps in gene_probes.items():
        best, bv = None, -1
        for p in ps:
            v = np.nanvar(probe_rows[p])
            if v > bv:
                best, bv = p, v
        gene_vec[g] = probe_rows[best]
    genes = sorted(gene_vec)
    X = np.column_stack([gene_vec[g] for g in genes])  # samples x genes
    print("gene-level matrix:", X.shape)

    y = np.array([1 if t == "Y" else 0 for t in tns])
    pid = np.array(pids)
    print("Y=%d N=%d" % (y.sum(), (1 - y).sum()))

    out = {"platform": "GPL14951", "n_samples": len(samples), "n_patients": len(set(pids)),
           "n_genes_tested": len(genes), "n_tumor": int(y.sum()), "n_normal": int((1 - y).sum())}

    # ---------- P1 原文口径：独立样本 t（54 样本） ----------
    p1 = stats.ttest_ind(X[y == 1], X[y == 0], axis=0, equal_var=True)[1]
    p1f = bh_fdr(p1)
    d1 = np.abs(stats.ttest_ind(X[y == 1], X[y == 0], axis=0, equal_var=True)[0])

    def auc_insample(Xm, ym):
        Xs = StandardScaler().fit_transform(Xm)
        lr = LogisticRegression(max_iter=5000, C=1.0)
        lr.fit(Xs, ym)
        return roc_auc_score(ym, lr.predict_proba(Xs)[:, 1])

    sel1 = np.where(p1 < 0.001)[0]
    a1 = auc_insample(X[:, sel1], y) if len(sel1) >= 2 else np.nan

    # ---------- P2/P3 配对口径：患者级 ----------
    # 每位患者：tumor 均值（多样本取平均）、normal 均值
    pat_t, pat_n, keep = {}, {}, []
    for pt in sorted(set(pid)):
        mt = (pid == pt) & (y == 1)
        mn = (pid == pt) & (y == 0)
        if mt.sum() >= 1 and mn.sum() >= 1:
            pat_t[pt] = np.nanmean(X[mt], axis=0)
            pat_n[pt] = np.nanmean(X[mn], axis=0)
            keep.append(pt)
    Pt = np.array([pat_t[k] for k in keep])
    Pn = np.array([pat_n[k] for k in keep])
    n_pairs = len(keep)
    print("paired patients:", n_pairs)
    t2, p2 = stats.ttest_rel(Pt, Pn, axis=0)
    p2 = np.asarray(p2, float)
    delta = Pt - Pn

    # 配对 AUC：以每对差值方向做「面板评分」，用配对 logistic（患者级）
    Xp = np.vstack([Pt, Pn])
    yp = np.array([1] * n_pairs + [0] * n_pairs)
    sel2 = np.where(p2 < 0.001)[0]
    a2 = auc_insample(Xp[:, sel2], yp) if len(sel2) >= 2 else np.nan

    out["P1_unpaired_paper_caliber"] = {
        "n_per_group": [int(y.sum()), int((1 - y).sum())],
        "n_de_p001": int((p1 < 0.001).sum()),
        "n_de_fdr05": int((p1f < 0.05).sum()),
        "n_fdr_sig_using_fdr_as_original_claimed": None,
        "auc_insample_de001": float(a1) if a1 == a1 else None,
        "top_genes": [genes[i] for i in np.argsort(p1)[:10]],
    }
    out["P2_paired_patient_level"] = {
        "n_pairs": n_pairs,
        "n_de_p001": int((p2 < 0.001).sum()),
        "n_de_p05": int((p2 < 0.05).sum()),
        "auc_insample_de001": float(a2) if a2 == a2 else None,
        "n_de05_by_patient_power": int((p2 < 0.05).sum()),
        "min_attainable_p_two_sided_n14": 2 / (2 ** n_pairs),
        "top_genes": [genes[i] for i in np.argsort(p2)[:10]],
        "mean_abs_delta_z": float(np.nanmean(np.abs(delta / (np.nanstd(delta, axis=0) + 1e-9)))),
    }

    # 伪重复膨胀因子
    y_counts = collections.Counter(pid[y == 1])
    out["pseudoreplication"] = {
        "n_tumor_samples": int(y.sum()),
        "n_patients_with_tumor": len(y_counts),
        "inflation_factor": float(y.sum() / len(y_counts)),
        "max_tumor_per_patient": int(max(y_counts.values())),
    }

    # 对比：显著基因数的口径差异
    out["comparison"] = {
        "n_de_P1": int((p1 < 0.001).sum()),
        "n_de_P2": int((p2 < 0.001).sum()),
        "auc_P1": float(a1) if a1 == a1 else None,
        "auc_P2": float(a2) if a2 == a2 else None,
        "overlap_genes_p001": sorted(set(genes[i] for i in np.where(p1 < 0.001)[0]) &
                                     set(genes[i] for i in np.where(p2 < 0.001)[0])),
    }

    json.dump(out, open(os.path.join(RES, "caseA_GSE65372_paired.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("\n=== 结果 ===")
    print("P1 原文口径(54): DE(p<0.001)=%d  AUC=%.4f" % (out["P1_unpaired_paper_caliber"]["n_de_p001"], a1))
    print("P2 配对口径(n=%d): DE(p<0.001)=%d  DE(p<0.05)=%d  AUC=%.4f" % (
        n_pairs, out["P2_paired_patient_level"]["n_de_p001"],
        out["P2_paired_patient_level"]["n_de_p05"], a2))
    print("伪重复膨胀因子: %.2f (%d 个瘤样本来自 %d 个患者)" % (
        out["pseudoreplication"]["inflation_factor"], y.sum(), len(y_counts)))
    print("wrote results/caseA_GSE65372_paired.json")


if __name__ == "__main__":
    main()
