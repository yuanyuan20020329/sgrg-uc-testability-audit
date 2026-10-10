#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Part C / Stage C — two-sample MR (SGRG cis-eQTL -> ulcerative colitis) plus
colocalisation, entirely in numpy/scipy.

Design
------
  exposure : GTEx v8 cis-eQTL (Colon_Sigmoid / Colon_Transverse / Brain_Cortex /
             Whole_Blood), instrument = independent genome-wide significant
             cis-SNP, beta = GTEx NES, se reconstructed from the t-statistic.
  outcome  : de Lange 2017 ulcerative colitis GWAS, GCST004133
             (12,366 cases / 33,609 controls; N = 45,975), harmonised file
             supplied by the GWAS Catalog.
  matching : on rsID (GTEx is GRCh38, the GWAS file is GRCh37) -> no liftover.
  alleles  : both sources report forward-strand REF/ALT, so a direct allele
             identity check is unambiguous and palindromic SNPs are safe when
             the alleles match exactly.

For every (gene, tissue) instrument set we report
    Wald ratio (lead SNP), IVW (fixed + multiplicative random effects),
    MR-Egger slope, MR-Egger intercept (directional pleiotropy),
    weighted median, weighted mode, Cochran's Q, MR-PRESSO-style global test,
    leave-one-out IVW, and a Steiger directionality check.
Colocalisation uses the classical Wakefield approximate Bayes factor with
coloc-style priors (p1=p2=1e-4, p12=1e-5, sd.prior=0.15).

Outputs (results/)
------------------
partC_mr_results.tsv        one row per (gene, tissue, method)
partC_mr_geneset.tsv        one row per (gene, tissue): primary IVW + sensitivity
partC_mr_summary.json       everything, machine readable
partC_coloc_results.tsv     PP.H0..H4 per (gene, tissue)
"""
import os
import csv
import json
import math
from collections import defaultdict

import numpy as np

import mr_core as M

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
RES = os.path.join(ROOT, "results")
INS = os.path.join(RES, "partC_instruments.tsv")
CIS = os.path.join(RES, "partC_eqtl_cis.tsv")
GWASH = os.path.join(RES, "partC_gwas_hits.tsv")

GTEX_N = {"Colon_Sigmoid": 318, "Colon_Transverse": 368,
          "Brain_Cortex": 205, "Whole_Blood": 670, "Brain_Hippocampus": 165}
UC_N = 45975

OUT_MR = os.path.join(RES, "partC_mr_results.tsv")
OUT_SET = os.path.join(RES, "partC_mr_geneset.tsv")
OUT_JSON = os.path.join(RES, "partC_mr_summary.json")
OUT_COLOC = os.path.join(RES, "partC_coloc_results.tsv")


def read_tsv(p):
    out = []
    with open(p, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out.append(r)
    return out


def fnum(x):
    try:
        if x is None or x == "" or x == "NA" or x == "nan":
            return None
        return float(x)
    except Exception:
        return None


# --------------------------------------------------------------------------- #
# coloc.abf (Wakefield approximate Bayes factor)
# --------------------------------------------------------------------------- #
def log_abf(beta, se, sd_prior):
    V = se ** 2
    W = sd_prior ** 2
    r = W / (W + V)
    z2 = (beta ** 2) / V
    return 0.5 * math.log(1 - r) + 0.5 * z2 * r


def coloc_abf(bx, sx, by, sy, sd_prior_x=0.15, sd_prior_y=0.15,
              p1=1e-4, p2=1e-4, p12=1e-5):
    l1 = np.array([log_abf(b, s, sd_prior_x) for b, s in zip(bx, sx)])
    l2 = np.array([log_abf(b, s, sd_prior_y) for b, s in zip(by, sy)])
    s1 = np.exp(l1 - l1.max()) * np.exp(l1.max())
    s1 = float(np.exp(l1).sum()); s2 = float(np.exp(l2).sum())
    s12 = float(np.exp(l1 + l2).sum())
    H0 = 1.0
    H1 = p1 * s1
    H2 = p2 * s2
    H3 = p12 * max(s1 * s2 - s12, 0.0)
    H4 = p12 * s12
    tot = H0 + H1 + H2 + H3 + H4
    pp = [x / tot for x in (H0, H1, H2, H3, H4)]
    return dict(zip(["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"], pp))


def main():
    ins = read_tsv(INS)
    sets = json.load(open(os.path.join(RES, "partC_instrument_sets.json"), encoding="utf-8"))
    gwas = {}
    for r in read_tsv(GWASH):
        rs = r.get("hm_rsid") or ""
        if not rs or rs == "NA":
            continue
        gwas[rs] = r
    print("instruments:", len(ins), "| GWAS rows pulled:", len(gwas))

    cis = read_tsv(CIS)
    cis_by = defaultdict(list)
    for r in cis:
        cis_by[(r["gene"], r["tissue"])].append(r)

    ins_by = defaultdict(list)
    for r in ins:
        ins_by[(r["gene"], r["tissue"])].append(r)

    mr_rows = []
    set_rows = []
    coloc_rows = []
    details = {}

    for key, iv in sorted(sets.items()):
        gene, tissue = key.split("|")
        recs = ins_by.get((gene, tissue), [])
        if not recs:
            continue
        # ---- merge exposure (instrument) with outcome, harmonise ----
        bx, sx, by, sy, used = [], [], [], [], []
        used_recs = []
        for r in recs:
            g = gwas.get(r["snp"])
            if not g:
                continue
            b_x = fnum(r["nes"]); s_x = fnum(r["se"])
            b_y = fnum(g.get("hm_beta")) if fnum(g.get("hm_beta")) is not None else fnum(g.get("beta"))
            s_y = fnum(g.get("standard_error"))
            if None in (b_x, s_x, b_y, s_y) or s_y <= 0 or s_x <= 0:
                continue
            hx, hy, st = M.harmonise_snp(
                b_x, s_x, r["ea"], r["oa"], b_y, s_y,
                g.get("hm_effect_allele") or g.get("effect_allele"),
                g.get("hm_other_allele") or g.get("other_allele"),
                fnum(g.get("hm_effect_allele_frequency")) or fnum(g.get("effect_allele_frequency")))
            if hx is None:
                mr_rows.append({"gene": gene, "tissue": tissue, "snp": r["snp"],
                                "method": "excluded", "note": st, "b": "", "se": "",
                                "p": "", "nsnp": 1})
                continue
            bx.append(hx); sx.append(s_x); by.append(hy); sy.append(s_y)
            used.append(r["snp"]); used_recs.append(r)
        k = len(bx)
        if k == 0:
            continue
        bx = np.array(bx); sx = np.array(sx); by = np.array(by); sy = np.array(sy)
        nx = GTEX_N.get(tissue, 300)
        r2x = M.r2_from_z(bx, sx, nx)
        r2y = M.r2_from_z(by, sy, UC_N)
        st = M.steiger(r2x, r2y)

        res = {"gene": gene, "tissue": tissue, "nsnp_used": k,
               "snps_used": used, "meanF": float(np.mean(M.f_stat(bx, sx))),
               "r2_exposure": r2x, "r2_outcome": r2y,
               "steiger_correct": st["correct_direction"], "methods": {}}

        # Wald ratio on the lead SNP (most significant exposure SNP used)
        lead = max(recs, key=lambda r: (-math.log10(max(fnum(r["p"]) or 1, 1e-300))))
        li = used.index(lead["snp"]) if lead["snp"] in used else 0
        wr = {"method": "Wald ratio", "nsnp": 1, "snp": used[li],
              "b": float(by[li] / bx[li]),
              "se": float(abs(sy[li] / bx[li]))}
        wr["p"] = float(2 * M._NORM.sf(abs(wr["b"] / wr["se"])))

        methods = [wr, M.ivw(bx, sx, by, sy)]
        if k >= 3:
            methods.append(M.mr_egger(bx, sx, by, sy))
            methods.append(M.weighted_median(bx, sx, by, sy))
            methods.append(M.weighted_mode(bx, sx, by, sy))
        for m in methods:
            res["methods"][m["method"]] = m
            mr_rows.append({"gene": gene, "tissue": tissue, "snp": m.get("snp", "ALL"),
                            "method": m["method"], "note": "",
                            "b": m.get("b"), "se": m.get("se"), "p": m.get("p"),
                            "nsnp": m.get("nsnp"),
                            "Q": m.get("Q") if "Q" in m else "",
                            "Q_p": m.get("Q_p") if "Q_p" in m else "",
                            "egger_int": m.get("egger_intercept", ""),
                            "egger_int_p": m.get("egger_intercept_p", "")})

        ivw = res["methods"]["IVW"]
        egg = res["methods"].get("MR-Egger")
        wm = res["methods"].get("Weighted median")
        presso = M.mr_presso_like(bx, sx, by, sy, nsim=1000) if k >= 3 else {}
        loo, _ = M.leave_one_out(bx, sx, by, sy)
        res["presso"] = presso
        res["loo"] = loo
        details[key] = {
            "bx": bx.tolist(), "sx": sx.tolist(), "by": by.tolist(),
            "sy": sy.tolist(), "snps": used,
            "exposure": [{"snp": u, "ea": rr["ea"], "oa": rr["oa"],
                          "nes": rr["nes"], "se": rr["se"],
                          "p": rr["p"], "pos": rr["pos"], "panel": rr.get("panel", "")}
                         for u, rr in zip(used, used_recs)],
        }
        set_rows.append({
            "gene": gene, "tissue": tissue, "nsnp": k, "meanF": round(res["meanF"], 1),
            "IVW_b": round(ivw["b"], 4), "IVW_se": round(ivw["se"], 4),
            "IVW_p": ivw["p"], "IVW_Q": round(ivw["Q"], 2) if ivw["Q"] == ivw["Q"] else "",
            "IVW_Q_p": ivw["Q_p"],
            "Egger_b": round(egg["b"], 4) if egg else "", "Egger_p": egg["p"] if egg else "",
            "Egger_int_p": egg.get("egger_intercept_p", "") if egg else "",
            "WMed_b": round(wm["b"], 4) if wm else "", "WMed_p": wm["p"] if wm else "",
            "PRESSO_global_p": presso.get("global_p_sim", ""),
            "PRESSO_n_outlier": presso.get("n_outlier", ""),
            "PRESSO_corrected_b": presso.get("corrected_b", ""),
            "PRESSO_corrected_p": presso.get("corrected_p", ""),
            "r2_exposure": r2x, "r2_outcome": r2y,
            "steiger_correct": st["correct_direction"],
            "wald_b": round(wr["b"], 4), "wald_p": wr["p"],
        })

        # ---- colocalisation over the whole cis window ----
        cw = cis_by.get((gene, tissue), [])
        cbx, csx, cby, csy = [], [], [], []
        for r in cw:
            g = gwas.get(r["snp"])
            if not g:
                continue
            b_x = fnum(r["nes"]); s_x = fnum(r["se"])
            b_y = fnum(g.get("hm_beta")) if fnum(g.get("hm_beta")) is not None else fnum(g.get("beta"))
            s_y = fnum(g.get("standard_error"))
            if None in (b_x, s_x, b_y, s_y) or s_y <= 0 or s_x <= 0:
                continue
            hx, hy, s2 = M.harmonise_snp(
                b_x, s_x, r["ea"], r["oa"], b_y, s_y,
                g.get("hm_effect_allele") or g.get("effect_allele"),
                g.get("hm_other_allele") or g.get("other_allele"),
                fnum(g.get("hm_effect_allele_frequency")) or fnum(g.get("effect_allele_frequency")))
            if hx is None:
                continue
            cbx.append(hx); csx.append(s_x); cby.append(hy); csy.append(s_y)
        if len(cbx) >= 5:
            pp = coloc_abf(np.array(cbx), np.array(csx), np.array(cby), np.array(csy))
            # prior sensitivity
            alt = {}
            for sp in (0.1, 0.2, 0.3):
                alt["PP.H4_sd%.2f" % sp] = coloc_abf(
                    np.array(cbx), np.array(csx), np.array(cby), np.array(csy),
                    sd_prior_x=sp, sd_prior_y=sp)["PP.H4"]
            row = {"gene": gene, "tissue": tissue, "n_snps": len(cbx)}
            row.update({k2: round(v, 4) for k2, v in pp.items()})
            row.update({k2: round(v, 4) for k2, v in alt.items()})
            coloc_rows.append(row)
            res["coloc"] = row

    # BH-FDR on primary IVW p across (gene,tissue)
    ps = [(r["IVW_p"], i) for i, r in enumerate(set_rows) if r["IVW_p"] not in ("", None)]
    ps_sorted = sorted(ps, key=lambda t: t[0])
    n = len(ps_sorted)
    prev = 1.0
    q = [1.0] * len(set_rows)
    for rank, (p, i) in enumerate(reversed(ps_sorted), 1):
        k = n - rank + 1
        prev = min(prev, p * n / k)
        q[i] = prev
    for i, r in enumerate(set_rows):
        r["IVW_q"] = q[i]

    with open(OUT_MR, "w", newline="", encoding="utf-8") as f:
        cols = ["gene", "tissue", "snp", "method", "note", "b", "se", "p",
                "nsnp", "Q", "Q_p", "egger_int", "egger_int_p"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in mr_rows:
            w.writerow({c: r.get(c, "") for c in cols})

    with open(OUT_SET, "w", newline="", encoding="utf-8") as f:
        cols = list(set_rows[0].keys())
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in set_rows:
            w.writerow(r)

    with open(OUT_COLOC, "w", newline="", encoding="utf-8") as f:
        if coloc_rows:
            cols = list(coloc_rows[0].keys())
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
            w.writeheader()
            for r in coloc_rows:
                w.writerow(r)

    json.dump({"sets": sets, "results": {k: {kk: vv for kk, vv in v.items() if kk != "loo"}
                                         for k, v in details.items()},
               "summary": set_rows, "coloc": coloc_rows},
              open(OUT_JSON, "w", encoding="utf-8"), indent=1, default=str)

    # ---- console summary ----
    print("\n=== MR: SGRG cis-eQTL -> UC (primary IVW) ===")
    print("%-10s %-18s %3s %8s %8s %10s %8s %6s %6s" %
          ("gene", "tissue", "k", "IVW_b", "IVW_se", "IVW_p", "IVW_q", "Q_p", "stgr"))
    for r in sorted(set_rows, key=lambda x: x["IVW_p"]):
        print("%-10s %-18s %3d %8.3f %8.3f %10.2e %8.3f %6s %6s" %
              (r["gene"], r["tissue"], r["nsnp"], r["IVW_b"], r["IVW_se"],
               r["IVW_p"], r["IVW_q"],
               ("%.2f" % r["IVW_Q_p"]) if isinstance(r["IVW_Q_p"], float) else "-",
               r["steiger_correct"]))
    print("\n=== colocalisation ===")
    for r in sorted(coloc_rows, key=lambda x: -x["PP.H4"]):
        if r["PP.H4"] > 0.05:
            print("  %-10s %-18s n=%4d PP.H4=%.3f PP.H3=%.3f"
                  % (r["gene"], r["tissue"], r["n_snps"], r["PP.H4"], r["PP.H3"]))
    print("\nfiles:", OUT_MR, "|", OUT_SET, "|", OUT_COLOC, "|", OUT_JSON)


if __name__ == "__main__":
    main()
