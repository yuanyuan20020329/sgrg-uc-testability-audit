# -*- coding: utf-8 -*-
"""
Part A · Step 4 — 独立队列验证 GSE11223 (Gemma 8385)
=====================================================
发现队列 : GSE38713 (Gemma 5873)  UC active inflamed (n=15) vs control (n=13)
验证队列 : GSE11223 (Gemma 8385)  UC inflamed colon (n≈62) vs Normal uninflamed colon (n≈63)
           结肠上皮活检；比发现队列大一个量级，且完全独立。

复用 Step1 的基因级聚合与 Welch t-test + BH-FDR 流程；额外做：
  * 两队列 log2FC 相关（Pearson/Spearman）
  * 方向一致率（全体 / 名义显著子集 / 发现集 FDR 子集）
  * 验证队列内的 SGRG 面板信号（面板内 FDR<0.05 计数、上/下调）
输出 : results/partA_step4_external_validation.tsv / .json
        results/fig_step4_crosscohort_log2FC.png
        results/fig_step4_crosscohort_concordance.png
"""
import csv
import json
import math
import os
import re
from collections import Counter

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
RAW = os.path.join(BASE, "data", "raw")
TMP = os.path.join(BASE, "data", "tmp")
OUT = os.path.join(BASE, "results")

MATRIX = os.path.join(RAW, "gemma_GSE11223.tsv")
SAMPLES = os.path.join(TMP, "samp", "8385.json")
SGRG_CLEAN = os.path.join(OUT, "sgrg_gene_set_clean.tsv")
DISCOVERY = os.path.join(OUT, "partA_GSE38713_SGRG_DE.tsv")

COLON = {"sigmoid colon", "descending colon", "ascending colon"}
COHORT = "GSE11223 / Gemma:8385"

log = lambda *a: print(*a, flush=True)


def load_sgrg():
    d = {}
    with open(SGRG_CLEAN, encoding="utf-8") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0]:
                d[p[0]] = p[1]
    log("[SGRG] HGNC-approved symbols = %d" % len(d))
    return d


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def load_matrix(path):
    """Return (header, list of (symbol_candidates, np.array of sample values))."""
    header, recs = None, []
    with open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            if len(parts) < 7:
                continue
            raw = (parts[2] or "").strip()
            cands = [x.strip() for x in raw.split("|") if x.strip()] if raw and raw != "NA" else []
            vals = np.array([_f(x) for x in parts[6:]], float)
            recs.append((cands, vals))
    return header, recs


def load_sample_groups(path):
    """BioAssayId -> (disease, inflammation, location, gsm, name)"""
    d = json.load(open(path, encoding="utf-8"))
    meta = {}
    for e in d:
        ch = {}
        for c in ((e.get("sample") or {}).get("characteristics") or []):
            ch[(c.get("category") or "").lower()] = c.get("value")
        meta[int(e["id"])] = dict(disease=ch.get("disease"), infl=ch.get("inflammation status"),
                                  loc=ch.get("anatomic location"), gsm=e.get("shortName"),
                                  name=e.get("name"))
    return meta


def bh_fdr(p):
    p = np.asarray(p, float)
    n = p.size
    o = np.argsort(p)
    r = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    q = np.empty(n)
    q[o] = np.clip(r, 0, 1)
    return q


def main():
    sgrg = load_sgrg()
    header, recs = load_matrix(MATRIX)
    log("[matrix] probes parsed = %d, columns = %d" % (len(recs), len(header)))

    all_cols = list(range(6, len(header)))               # positions in header
    if not recs:
        raise SystemExit("no matrix rows parsed")
    n_expect = len(header) - 6
    good = [r for r in recs if len(r[1]) == n_expect]
    bad = len(recs) - len(good)
    if bad:
        log("[matrix] skipped %d malformed row(s) (expected %d sample values)" % (bad, n_expect))
    recs = good
    V = np.vstack([r[1] for r in recs])                  # n_probes x n_samples
    n_sam = V.shape[1]
    alive = [k for k in range(n_sam) if not np.all(np.isnan(V[:, k]))]
    dead = [k for k in range(n_sam) if k not in alive]
    if dead:
        log("[matrix] dropping all-NaN sample columns: %d (%s)"
            % (len(dead), [header[6 + k][-30:] for k in dead[:4]]))
    V = V[:, alive]
    sample_cols = [all_cols[k] for k in alive]

    meta = load_sample_groups(SAMPLES)

    # --- map each retained matrix column to a group --------------------------
    groups, matched = [], 0
    for j in sample_cols:
        m = re.search(r"BioAssayId=(\d+)", header[j])
        bid = int(m.group(1)) if m else None
        info = meta.get(bid)
        if not info:
            groups.append("Unknown")
            continue
        matched += 1
        dis, inf, loc = info["disease"], info["infl"], info["loc"]
        if dis == "UC" and inf == "Inflamed" and loc in COLON:
            g = "UC_inflamed_colon"
        elif dis == "Normal" and inf == "Uninflamed" and loc in COLON:
            g = "Normal_uninflamed_colon"
        elif loc == "terminal ileum":
            g = "ileum_excluded"
        else:
            g = "other_excluded"
        groups.append(g)
    log("[map] columns=%d matched=%d  %s" % (len(sample_cols), matched, dict(Counter(groups))))

    # --- gene-level collapse (nan-aware) -------------------------------------
    gene_rows = {}
    for (cands, _), k in zip(recs, range(len(recs))):
        for c in cands:
            gene_rows.setdefault(c, []).append(k)
    expr = {}
    for g, ridx in gene_rows.items():
        expr[g] = np.nanmean(V[ridx, :], axis=0)
    log("[expr] genes = %d" % len(expr))

    hit = sorted(set(expr) & set(sgrg))
    log("[SGRG] measured = %d / %d" % (len(hit), len(sgrg)))

    def idx(g):
        return [k for k, gg in enumerate(groups) if gg == g]

    A, B = idx("UC_inflamed_colon"), idx("Normal_uninflamed_colon")
    log("[contrast] UC inflamed colon n=%d vs Normal uninflamed colon n=%d" % (len(A), len(B)))

    res = []
    for g in hit:
        v = expr[g]
        a = v[A][~np.isnan(v[A])]
        b = v[B][~np.isnan(v[B])]
        if len(a) < 3 or len(b) < 3:
            continue
        t, p = stats.ttest_ind(a, b, equal_var=False)
        res.append(dict(gene=g, log2FC=float(np.mean(a) - np.mean(b)), t=float(t), p=float(p),
                        mean_case=float(np.mean(a)), mean_ctrl=float(np.mean(b)),
                        sd_case=float(np.std(a, ddof=1)), sd_ctrl=float(np.std(b, ddof=1))))
    q = bh_fdr([r["p"] for r in res])
    for r, qq in zip(res, q):
        r["FDR"] = float(qq)
    res.sort(key=lambda r: r["p"])
    vmap = {r["gene"]: r for r in res}

    n_sig = sum(1 for r in res if r["FDR"] < 0.05)
    n_up = sum(1 for r in res if r["FDR"] < 0.05 and r["log2FC"] > 0)
    n_dn = sum(1 for r in res if r["FDR"] < 0.05 and r["log2FC"] < 0)
    n_nom = sum(1 for r in res if r["p"] < 0.05)
    log("[DE] UC inflamed vs normal : p<0.05 = %d | FDR<0.05 = %d (up %d / down %d)"
        % (n_nom, n_sig, n_up, n_dn))

    # --- cross-cohort comparison with discovery ------------------------------
    disc = list(csv.DictReader(open(DISCOVERY, encoding="utf-8"), delimiter="\t"))
    dmap = {r["GeneSymbol"]: r for r in disc}
    disc_fdr = {r["GeneSymbol"] for r in disc if float(r["FDR_inflamed"]) < 0.05}
    log("[discovery] SGRG measured=%d  FDR<0.05=%d" % (len(disc), len(disc_fdr)))

    common = sorted(set(dmap) & set(vmap))
    dx = np.array([float(dmap[g]["log2FC_inflamed_vs_ctrl"]) for g in common])
    vy = np.array([vmap[g]["log2FC"] for g in common])
    pear = stats.pearsonr(dx, vy)
    spear = stats.spearmanr(dx, vy)
    agree = int(np.sum(np.sign(dx) * np.sign(vy) > 0))
    log("[cross] common genes=%d  Pearson r=%.3f (p=%.2g)  Spearman rho=%.3f  sign-agree=%d (%.1f%%)"
        % (len(common), pear.statistic, pear.pvalue, spear.statistic, agree, 100.0 * agree / len(common)))

    # replication among discovery-FDR-significant genes
    d_sig = [g for g in common if g in disc_fdr]
    rep_nom = [g for g in d_sig if vmap[g]["p"] < 0.05 and
               np.sign(float(dmap[g]["log2FC_inflamed_vs_ctrl"])) * np.sign(vmap[g]["log2FC"]) > 0]
    rep_fdr = [g for g in d_sig if vmap[g]["FDR"] < 0.05 and
               np.sign(float(dmap[g]["log2FC_inflamed_vs_ctrl"])) * np.sign(vmap[g]["log2FC"]) > 0]
    log("[replication] discovery FDR genes present in validation = %d" % len(d_sig))
    log("              same-direction p<0.05 = %d (%.1f%%)" % (len(rep_nom), 100.0 * len(rep_nom) / max(1, len(d_sig))))
    log("              same-direction FDR<0.05 = %d (%.1f%%)" % (len(rep_fdr), 100.0 * len(rep_fdr) / max(1, len(d_sig))))
    # expected by chance (direction) ~50%
    binom = stats.binomtest(len(rep_nom), len(d_sig), 0.5, alternative="greater") if d_sig else None
    if binom:
        log("              binomial vs 50%%: p=%.3g" % binom.pvalue)

    # --- write result table ---------------------------------------------------
    p1 = os.path.join(OUT, "partA_step4_external_validation.tsv")
    with open(p1, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["GeneSymbol", "HR_evidence", "log2FC_validation", "p_validation", "FDR_validation",
                    "mean_case", "mean_ctrl", "log2FC_discovery", "FDR_discovery", "direction_match",
                    "is_discovery_FDR"])
        for r in res:
            g = r["gene"]
            d = dmap.get(g)
            dm = ""
            if d:
                dm = "True" if np.sign(r["log2FC"]) * np.sign(float(d["log2FC_inflamed_vs_ctrl"])) > 0 else "False"
            w.writerow([g, sgrg.get(g, ""), "%.4f" % r["log2FC"], "%.3e" % r["p"], "%.3e" % r["FDR"],
                        "%.4f" % r["mean_case"], "%.4f" % r["mean_ctrl"],
                        ("%.4f" % float(d["log2FC_inflamed_vs_ctrl"])) if d else "",
                        ("%.3e" % float(d["FDR_inflamed"])) if d else "",
                        dm, "True" if g in disc_fdr else "False"])
    log("[write]", p1)

    # --- figures --------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.6, 5.4), dpi=300)
    ax.scatter(dx, vy, s=20, c="#4b5563", alpha=.7, lw=0)
    lim = max(1.5, float(np.max(np.abs(np.r_[dx, vy]))) * 1.05)
    ax.plot([-lim, lim], [-lim, lim], ls="--", lw=.8, c="#9aa1a9")
    ax.axhline(0, lw=.5, c="#c9ced4"); ax.axvline(0, lw=.5, c="#c9ced4")
    # highlight top discovery genes
    for g in sorted(disc_fdr)[:12]:
        if g in vmap:
            ax.annotate(g, (float(dmap[g]["log2FC_inflamed_vs_ctrl"]), vmap[g]["log2FC"]),
                        fontsize=6.5, xytext=(2, 2), textcoords="offset points")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("log2FC  GSE38713 (discovery)", fontsize=9)
    ax.set_ylabel("log2FC  GSE11223 (validation)", fontsize=9)
    ax.set_title("SGRG log2FC cross-cohort\nPearson r=%.2f  Spearman rho=%.2f  (n=%d)"
                 % (pear.statistic, spear.statistic, len(common)), fontsize=9)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    pv = os.path.join(OUT, "fig_step4_crosscohort_log2FC.png")
    fig.savefig(pv, dpi=300); plt.close(fig)
    log("[write]", pv)

    d_sig_present = len(d_sig)
    fig, ax = plt.subplots(figsize=(5.2, 3.2), dpi=300)
    bars = [("same-direction\np<0.05", len(rep_nom)),
            ("same-direction\nFDR<0.05", len(rep_fdr)),
            ("total discovery\nFDR genes", d_sig_present)]
    ax.bar([b[0] for b in bars], [b[1] for b in bars], color=["#0b5cab", "#b42318", "#9aa1a9"], width=.6)
    for i, b in enumerate(bars):
        ax.text(i, b[1] + .3, str(b[1]), ha="center", fontsize=9)
    ax.set_ylabel("SGRG count", fontsize=9)
    ax.set_title("GSE11223 replication of GSE38713 FDR-significant SGRG", fontsize=9)
    ax.tick_params(labelsize=8)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    pv2 = os.path.join(OUT, "fig_step4_crosscohort_concordance.png")
    fig.savefig(pv2, dpi=300); plt.close(fig)
    log("[write]", pv2)

    # --- summary json ---------------------------------------------------------
    summ = dict(
        cohort=COHORT, tissue="colonic epithelial biopsy (sigmoid/descending/ascending colon)",
        platform="Affymetrix (GPL570 &c.)",
        contrast="UC_inflamed_colon (n=%d) vs Normal_uninflamed_colon (n=%d)" % (len(A), len(B)),
        group_counts=dict(Counter(groups)),
        sgrg_measured=len(hit),
        validation_nominal_p05=n_nom, validation_FDR05=n_sig, validation_up=n_up, validation_down=n_dn,
        cross_common_genes=len(common),
        cross_pearson_r=round(float(pear.statistic), 4), cross_pearson_p=float(pear.pvalue),
        cross_spearman_rho=round(float(spear.statistic), 4),
        sign_agreement=agree, sign_agreement_pct=round(100.0 * agree / len(common), 2),
        discovery_FDR_genes_presentin_validation=d_sig_present,
        replication_same_direction_p05=len(rep_nom),
        replication_same_direction_FDR05=len(rep_fdr),
        replication_binom_p=(float(binom.pvalue) if binom else None),
        replicated_genes_nominal=sorted(rep_nom),
        replicated_genes_fdr=sorted(rep_fdr),
        top20_validation=[dict(gene=r["gene"], log2FC=round(r["log2FC"], 3),
                               p=float('%.3e' % r["p"]), FDR=float('%.3e' % r["FDR"])) for r in res[:20]],
    )
    ps = os.path.join(OUT, "partA_step4_external_validation.json")
    json.dump(summ, open(ps, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    log("[write]", ps)
    log("[done]")


if __name__ == "__main__":
    main()
