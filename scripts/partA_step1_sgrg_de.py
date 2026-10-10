# -*- coding: utf-8 -*-
"""
Part A · Step 1 (v2)  —  UC (GSE38713) 应激颗粒相关基因(SGRG)差异表达分析
=====================================================================
数据源 : Gemma REST API v3, dataset id=5873 (GSE38713)
          "Expression data from intestinal mucosa of patients with UC"
          human colonic mucosal biopsy; Affymetrix HG-U133 Plus 2.0
          Gemma 已完成 RMA + quantile normalization，值为 log2 尺度。
SGRG集 : GO:0010494 (cytoplasmic stress granule) U GO:0035061 (nuclear stress granule)
          Homo sapiens；经 HGNC 官方 REST API 校验，仅保留批准的官方符号 (n=104)
对比   : (主) UC active inflamed (n=15) vs non-inflammatory control (n=13)
          (辅) all-UC (n=30) vs control (n=13)
统计   : Welch t-test on log2 expression；Benjamini-Hochberg FDR（在 SGRG 面板内校正）
输出   : results/ 下的结果表 + 火山图 + 热图 + 各组均值 + 摘要 JSON
"""
import gzip, json, os, re, sys, csv, math
from collections import Counter
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
RAW  = os.path.join(BASE, "data", "raw")
OUT  = os.path.join(BASE, "results")
os.makedirs(OUT, exist_ok=True)

MATRIX_FILE = os.path.join(RAW, "gemma_GSE38713.tsv")
SAMPLES     = os.path.join(RAW, "s_100.json")
SGRG_CLEAN  = os.path.join(OUT, "sgrg_gene_set_clean.tsv")
DS_ID       = "Gemma:5873 / GSE38713"

log = lambda *a: print(*a, flush=True)

# ---------------------------------------------------------------- SGRG 集
def load_sgrg():
    """returns {symbol: evidence_string}"""
    d = {}
    with open(SGRG_CLEAN, encoding="utf-8") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0]:
                d[p[0]] = p[1]
    log(f"[SGRG] HGNC-approved symbols = {len(d)}")
    return d

# ---------------------------------------------------------------- 矩阵
def load_matrix(path):
    header, rows = None, []
    with open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts; continue
            if len(parts) < 7:
                continue
            rows.append(parts)
    return header, rows

def col_label(c):
    m = re.search(r"Name=([^=]*)$", c)
    return (m.group(1).strip() if m else c)

def norm_label(s):
    return re.sub(r"[-+.]+$", "", (s or "").strip())

# ---------------------------------------------------------------- 分组
def load_groups(path):
    d = json.load(open(path, encoding="utf-8"))
    lab2grp, lab2gsm = {}, {}
    for it in d.get("data", []):
        lab, desc = (it.get("name") or "").strip(), (it.get("description") or "").lower()
        if "non-inflammatory control" in desc:      g = "Control"
        elif "non-involved area uc active" in desc: g = "UC_active_noninflamed"
        elif "involved area uc active" in desc:     g = "UC_active_inflamed"
        elif "in remission" in desc:                g = "UC_remission"
        else:                                       g = "Unknown"
        lab2grp[lab] = g
        lab2gsm[lab] = it.get("shortName")
    return lab2grp, lab2gsm

def bh_fdr(p):
    p = np.asarray(p, float); n = p.size
    o = np.argsort(p)
    r = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    q = np.empty(n); q[o] = np.clip(r, 0, 1)
    return q

# ================================================================ main
def main():
    sgrg = load_sgrg()
    header, rows = load_matrix(MATRIX_FILE)
    i_sym = header.index("GeneSymbol") if "GeneSymbol" in header else 2

    sample_cols = list(range(6, len(header)))
    labels = [col_label(header[j]) for j in sample_cols]
    lab2grp, lab2gsm = load_groups(SAMPLES)
    def map_label(l):
        # 精确优先；Gemma 列名把样本名末尾的 '-' 转义成了 '.'，故定向回退
        if l in lab2grp:
            return lab2grp[l]
        return lab2grp.get(re.sub(r"\.$", "-", l), "Unknown")
    groups = [map_label(l) for l in labels]
    log("[groups]", dict(Counter(groups)))
    if any(g == "Unknown" for g in groups):
        log("[warn] unmapped:", [l for l, g in zip(labels, groups) if g == "Unknown"])

    # 基因级聚合：多基因共享探针 (A|B|C) 登记到每个候选符号
    gene_vals = {}
    gene_probes = {}
    n_multi = 0
    for r in rows:
        raw = (r[i_sym] or "").strip()
        if not raw or raw == "NA":
            continue
        cands = [x.strip() for x in raw.split("|") if x.strip()]
        if len(cands) > 1:
            n_multi += 1
        try:
            v = np.array([float(r[c]) for c in sample_cols], float)
        except (ValueError, IndexError):
            continue
        for c in cands:
            gene_vals.setdefault(c, []).append(v)
            gene_probes.setdefault(c, []).append((r[0], len(cands)))
    expr = {g: np.nanmean(np.vstack(vs), axis=0) for g, vs in gene_vals.items()}
    log(f"[expr] genes = {len(expr)} (multi-gene probes expanded: {n_multi})")

    hit = sorted(set(expr) & set(sgrg))
    miss = sorted(set(sgrg) - set(expr))
    log(f"[SGRG] measured = {len(hit)} / {len(sgrg)}")

    # 分组索引
    def idx(g): return [k for k, gg in enumerate(groups) if gg == g]
    def contrast(nameA, nameB, lstA, lstB):
        res = []
        for g in hit:
            v = expr[g]; a, b = v[lstA], v[lstB]
            t, p = stats.ttest_ind(a, b, equal_var=False)
            res.append(dict(gene=g, log2FC=float(np.mean(a) - np.mean(b)), t=float(t), p=float(p),
                            mean_A=float(np.mean(a)), mean_B=float(np.mean(b)),
                            sd_A=float(np.std(a, ddof=1)), sd_B=float(np.std(b, ddof=1))))
        q = bh_fdr([r["p"] for r in res])
        for r, qq in zip(res, q): r["FDR"] = float(qq)
        res.sort(key=lambda r: r["p"])
        return res

    ct = idx("Control")
    uc_inf = idx("UC_active_inflamed")
    uc_all = idx("UC_active_inflamed") + idx("UC_active_noninflamed") + idx("UC_remission")
    log(f"[contrast] inflamed n={len(uc_inf)} | all-UC n={len(uc_all)} | control n={len(ct)}")

    main_res = contrast("UC_active_inflamed", "Control", uc_inf, ct)
    all_res  = contrast("UC_all", "Control", uc_all, ct)

    n_sig = sum(1 for r in main_res if r["FDR"] < 0.05)
    n_up = sum(1 for r in main_res if r["FDR"] < 0.05 and r["log2FC"] > 0)
    n_dn = sum(1 for r in main_res if r["FDR"] < 0.05 and r["log2FC"] < 0)
    n_nom = sum(1 for r in main_res if r["p"] < 0.05)
    log(f"[DE] inflamed vs control : p<0.05 = {n_nom} | FDR<0.05 = {n_sig} (up {n_up} / down {n_dn})")

    by_gene = {r["gene"]: r for r in all_res}

    # ---- 结果表 ----
    p1 = os.path.join(OUT, "partA_GSE38713_SGRG_DE.tsv")
    with open(p1, "w", encoding="utf-8") as f:
        f.write("\t".join(["GeneSymbol", "HR_evidence", "n_probes", "probe_quality",
                           "log2FC_inflamed_vs_ctrl", "p_inflamed",
                           "FDR_inflamed", "log2FC_allUC_vs_ctrl", "p_allUC", "FDR_allUC",
                           "mean_Control", "mean_UC_inflamed"]) + "\n")
        for r in main_res:
            a = by_gene.get(r["gene"], {})
            prs = gene_probes.get(r["gene"], [])
            pq = "unique" if prs and all(n == 1 for _, n in prs) else "multi-gene_probe"
            f.write("\t".join([r["gene"], sgrg[r["gene"]], str(len(prs)), pq,
                               f'{r["log2FC"]:.4f}', f'{r["p"]:.3e}', f'{r["FDR"]:.3e}',
                               f'{a.get("log2FC", float("nan")):.4f}', f'{a.get("p", float("nan")):.3e}',
                               f'{a.get("FDR", float("nan")):.3e}',
                               f'{r["mean_B"]:.4f}', f'{r["mean_A"]:.4f}']) + "\n")
    log("[write]", p1)

    # ---- 各组均值矩阵 ----
    gnames = ["Control", "UC_remission", "UC_active_noninflamed", "UC_active_inflamed"]
    p2 = os.path.join(OUT, "partA_GSE38713_SGRG_groupMeans.tsv")
    with open(p2, "w", encoding="utf-8") as f:
        f.write("GeneSymbol\t" + "\t".join(f"{g}(n={len(idx(g))})" for g in gnames) + "\n")
        for r in main_res:
            f.write(r["gene"] + "\t" + "\t".join(f'{np.mean(expr[r["gene"]][idx(g)]):.4f}' for g in gnames) + "\n")
    log("[write]", p2)

    # ---- 火山图 ----
    fig, ax = plt.subplots(figsize=(6.2, 5.4), dpi=300)
    xs = np.array([r["log2FC"] for r in main_res]); ys = -np.log10([r["p"] for r in main_res])
    sig = np.array([r["FDR"] < 0.05 for r in main_res])
    up, dn = sig & (xs > 0), sig & (xs < 0)
    ax.scatter(xs[~sig], ys[~sig], s=18, c="#b9bfc7", alpha=.85, lw=0)
    ax.scatter(xs[up], ys[up], s=30, c="#b42318", alpha=.9, lw=0, label=f"Up, FDR<0.05 (n={up.sum()})")
    ax.scatter(xs[dn], ys[dn], s=30, c="#0b5cab", alpha=.9, lw=0, label=f"Down, FDR<0.05 (n={dn.sum()})")
    for r in main_res[:10]:
        ax.annotate(r["gene"], (r["log2FC"], -math.log10(r["p"])), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    ax.axhline(-math.log10(0.05), ls="--", lw=.8, c="#8a9199")
    ax.axvline(0, ls="-", lw=.6, c="#c9ced4")
    ax.set_xlabel("log2 fold change  (UC active inflamed vs control)", fontsize=9)
    ax.set_ylabel("−log10 P  (Welch t-test)", fontsize=9)
    ax.set_title("GSE38713 · stress granule-related genes (n=%d)\nUC active inflamed (n=%d) vs non-inflammatory control (n=%d)"
                 % (len(hit), len(uc_inf), len(ct)), fontsize=8.8)
    ax.legend(fontsize=8, frameon=False); ax.tick_params(labelsize=8)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    fig.tight_layout()
    pv = os.path.join(OUT, "fig_volcano_SGRG_GSE38713.png"); fig.savefig(pv, dpi=300); plt.close(fig)
    log("[write]", pv)

    # ---- 热图 top30 ----
    top = [r["gene"] for r in main_res[:30]]
    order = uc_inf + idx("UC_active_noninflamed") + idx("UC_remission") + ct
    M = np.vstack([expr[g][order] for g in top])
    M = (M - M.mean(axis=1, keepdims=True)) / (M.std(axis=1, keepdims=True) + 1e-9)
    fig, ax = plt.subplots(figsize=(7.8, 7.2), dpi=300)
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top, fontsize=7)
    ax.set_xticks([])
    b = [0, len(uc_inf), len(uc_inf)+len(idx("UC_active_noninflamed")),
         len(uc_inf)+len(idx("UC_active_noninflamed"))+len(idx("UC_remission")), M.shape[1]]
    for x in b[1:-1]: ax.axvline(x - .5, color="k", lw=1.0)
    ax.set_xlabel("inflamed-UC | non-inflamed UC | remission | control", fontsize=8.5)
    ax.set_title("Top 30 SGRG by P (row z-score)", fontsize=9.5)
    cb = fig.colorbar(im, ax=ax, fraction=.03, pad=.02); cb.set_label("z-score", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.tight_layout()
    ph = os.path.join(OUT, "fig_heatmap_SGRG_GSE38713.png"); fig.savefig(ph, dpi=300); plt.close(fig)
    log("[write]", ph)

    # ---- summary ----
    summ = dict(
        dataset=DS_ID, tissue="human colonic mucosal biopsy", platform="Affymetrix HG-U133 Plus 2.0",
        sgrg_source="GO:0010494 U GO:0035061 (Homo sapiens), HGNC-validated",
        sgrg_total=len(sgrg), sgrg_measured=len(hit), sgrg_not_measured=len(miss),
        sgrg_not_measured_symbols=miss,
        group_counts=dict(Counter(groups)),
        contrast_main=f"UC_active_inflamed (n={len(uc_inf)}) vs Control (n={len(ct)})",
        contrast_secondary=f"all-UC (n={len(uc_all)}) vs Control (n={len(ct)})",
        method="Welch t-test on log2 expression; BH-FDR within SGRG panel",
        main_nominal_p05=n_nom, main_FDR05=n_sig, main_up=n_up, main_down=n_dn,
        top20=[dict(gene=r["gene"], log2FC=round(r["log2FC"], 3),
                    p=float(f'{r["p"]:.3e}'), FDR=float(f'{r["FDR"]:.3e}'),
                    evidence=sgrg[r["gene"]]) for r in main_res[:20]],
    )
    ps = os.path.join(OUT, "partA_GSE38713_SGRG_summary.json")
    json.dump(summ, open(ps, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    log("[write]", ps)
    log("[done]")

if __name__ == "__main__":
    main()
