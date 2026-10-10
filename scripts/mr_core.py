#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Part C — minimal, dependency-light two-sample Mendelian randomisation toolkit.

No R, no TwoSampleMR: every estimator below is implemented directly on
(bx, sx, by, sy) arrays with numpy/scipy, so results are fully transparent and
reproducible in this sandbox.

Estimators
----------
wald_ratio        single-instrument ratio
ivw               inverse-variance weighted (fixed + multiplicative random effects)
mr_egger          weighted linear regression with intercept (directional pleiotropy)
weighted_median   weighted median estimator (bootstrap SE)
weighted_mode     weighted mode estimator (bootstrap SE)
cochran_q         Cochran's Q heterogeneity
mr_presso_like    global + outlier-corrected IVW (approximate MR-PRESSO)
leave_one_out     LOO IVW
steiger           directionality (R^2 = z^2 / N)

Conventions
-----------
* bx = effect of the SNP on the exposure (eQTL), by = effect on the outcome (UC).
* Both betas always refer to the SAME effect allele after harmonisation.
* All betas on the standardised (log-expression / log-odds) scale.
"""
import numpy as np
from scipy import stats

_CHI2 = stats.chi2
_NORM = stats.norm
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


# --------------------------------------------------------------------------- #
# harmonisation
# --------------------------------------------------------------------------- #
def complement(a):
    return _COMP.get((a or "").upper(), None)


def harmonise_snp(bx, sx, ea_x, oa_x, by, sy, ea_y, oa_y, eaf_y=None,
                  palindromic_eaf_band=(0.42, 0.58)):
    """Return (bx, by, status). Both betas end up relative to ea_x (eQTL alt)."""
    ea_x, oa_x = (ea_x or "").upper(), (oa_x or "").upper()
    ea_y, oa_y = (ea_y or "").upper(), (oa_y or "").upper()
    if not (ea_x and oa_x and ea_y and oa_y):
        return None, None, "missing_allele"
    # direct match
    if {ea_x, oa_x} == {ea_y, oa_y}:
        if ea_x == ea_y:
            return bx, by, "match"
        return bx, -by, "flip"
    # strand flip on the outcome side
    cea, coa = complement(ea_y), complement(oa_y)
    if {ea_x, oa_x} == {cea, coa}:
        if ea_x == cea:
            return bx, by, "strand_match"
        return bx, -by, "strand_flip"
    # palindromic ambiguity (A/T or C/G) — only resolvable with EAF
    if {ea_x, oa_x} == {"A", "T"} or {ea_x, oa_x} == {"C", "G"}:
        if eaf_y is None:
            return None, None, "palindromic_no_eaf"
        e = float(eaf_y)
        if palindromic_eaf_band[0] <= e <= palindromic_eaf_band[1]:
            return None, None, "palindromic_ambiguous"
        # if the outcome effect allele is the *minor/derived* one we can infer
        # strand: for palindromic SNPs, a big deviation from 0.5 lets us assign
        # the strand when the eQTL allele frequency is assumed closer to the
        # GWAS effect-allele frequency.  Keep it simple + conservative: accept
        # only when the outcome effect allele equals the eQTL effect allele.
        if ea_x == ea_y:
            return bx, by, "palindromic_match"
        if ea_x == oa_y:
            return bx, -by, "palindromic_flip"
        return None, None, "palindromic_unresolved"
    return None, None, "allele_mismatch"


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def f_stat(bx, sx):
    """Per-SNP F = (beta/se)^2 (1 df); weak-instrument rule of thumb: mean F>10."""
    bx = np.asarray(bx, float)
    sx = np.asarray(sx, float)
    return (bx / sx) ** 2


def ratio(bx, sx, by, sy):
    b = by / bx
    se = np.abs(sy / bx)
    return b, se


def cochran_q(b, se, beta=None):
    """Cochran's Q for heterogeneity among ratio estimates (fixed effect)."""
    w = 1.0 / se ** 2
    if beta is None:
        beta = np.sum(w * b) / np.sum(w)
    Q = float(np.sum(w * (b - beta) ** 2))
    df = len(b) - 1
    p = float(_CHI2.sf(Q, df)) if df > 0 else np.nan
    return Q, df, p


# --------------------------------------------------------------------------- #
# estimators
# --------------------------------------------------------------------------- #
def ivw(bx, sx, by, sy, random_effects=True):
    b, se = ratio(bx, sx, by, sy)
    w = 1.0 / se ** 2
    beta = float(np.sum(w * b) / np.sum(w))
    se_fixed = float(np.sqrt(1.0 / np.sum(w)))
    se_used = se_fixed
    infl = 1.0
    Q, df, pQ = cochran_q(b, se, beta)
    if random_effects and df > 0 and Q > df:
        infl = float(np.sqrt(Q / df))
        se_used = se_fixed * infl
    z = beta / se_used if se_used > 0 else np.nan
    p = float(2 * _NORM.sf(abs(z))) if np.isfinite(z) else np.nan
    return {"method": "IVW", "nsnp": len(b), "b": beta, "se": se_used,
            "p": p, "beta_fixed_se": se_fixed, "infl": infl,
            "Q": Q, "Q_df": df, "Q_p": pQ}


def mr_egger(bx, sx, by, sy):
    bx = np.asarray(bx, float)
    by = np.asarray(by, float)
    sy = np.asarray(sy, float)
    # degenerate case: all instruments carry the same exposure effect (perfect
    # LD proxies) -> the weighted design matrix is singular and Egger is not
    # identifiable.
    if len(np.unique(np.round(bx, 10))) < 2:
        return {"method": "MR-Egger", "nsnp": len(bx), "b": np.nan, "se": np.nan,
                "p": np.nan, "egger_intercept": np.nan,
                "egger_intercept_se": np.nan, "egger_intercept_p": np.nan,
                "note": "singular: no variation in bx"}
    w = 1.0 / sy ** 2
    # normalise the weights so the normal equations stay well conditioned
    w = w / np.max(w)
    X = np.column_stack([np.ones_like(bx), bx])
    W = np.diag(w)
    XtWX = X.T @ W @ X
    try:
        cov = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(XtWX)
    coef = cov @ (X.T @ W @ by)
    resid = by - X @ coef
    dof = max(len(bx) - 2, 1)
    sigma2 = float((resid ** 2 * w).sum() / dof)
    se = np.sqrt(np.abs(np.diag(cov) * sigma2))
    intercept, slope = float(coef[0]), float(coef[1])
    p_slope = float(2 * _NORM.sf(abs(slope / se[1]))) if se[1] > 0 else np.nan
    p_int = float(2 * _NORM.sf(abs(intercept / se[0]))) if se[0] > 0 else np.nan
    return {"method": "MR-Egger", "nsnp": len(bx), "b": slope, "se": float(se[1]),
            "p": p_slope, "egger_intercept": intercept,
            "egger_intercept_se": float(se[0]), "egger_intercept_p": p_int}


def weighted_median(bx, sx, by, sy, nboot=10000, seed=20261008):
    b, se = ratio(bx, sx, by, sy)
    w = 1.0 / se ** 2

    def _wm(bb, ww):
        order = np.argsort(bb)
        bb, ww = bb[order], ww[order]
        cw = np.cumsum(ww) / np.sum(ww)
        k = np.searchsorted(cw, 0.5)
        k = min(max(k, 0), len(bb) - 1)
        if np.isclose(cw[k], 0.5) and k + 1 < len(bb):
            return 0.5 * (bb[k] + bb[k + 1])
        return bb[k]

    beta = float(_wm(b, w))
    rng = np.random.default_rng(seed)
    n = len(b)
    if n < 3:
        return {"method": "Weighted median", "nsnp": n, "b": beta,
                "se": np.nan, "p": np.nan}
    boots = np.empty(nboot)
    for i in range(nboot):
        idx = rng.integers(0, n, n)
        boots[i] = _wm(b[idx], w[idx])
    se = float(np.std(boots, ddof=1))
    p = float(2 * _NORM.sf(abs(beta / se))) if se > 0 else np.nan
    return {"method": "Weighted median", "nsnp": n, "b": beta, "se": se, "p": p}


def weighted_mode(bx, sx, by, sy, phi=1.0, nboot=2000, seed=20261008):
    b, se = ratio(bx, sx, by, sy)
    w = 1.0 / se ** 2
    w = w / w.sum()

    def _mode(bb, ww):
        if len(bb) == 1:
            return float(bb[0])
        bw = phi * np.std(bb, ddof=1)
        if bw <= 0:
            return float(bb[np.argmax(ww)])
        grid = np.linspace(bb.min() - 2 * bw, bb.max() + 2 * bw, 512)
        dens = np.zeros_like(grid)
        for v, wt in zip(bb, ww):
            dens += wt * np.exp(-0.5 * ((grid - v) / bw) ** 2)
        return float(grid[np.argmax(dens)])

    beta = _mode(b, w)
    rng = np.random.default_rng(seed)
    n = len(b)
    if n < 3:
        return {"method": "Weighted mode", "nsnp": n, "b": beta, "se": np.nan, "p": np.nan}
    boots = np.empty(nboot)
    for i in range(nboot):
        idx = rng.integers(0, n, n)
        ww = w[idx]
        ww = ww / ww.sum()
        boots[i] = _mode(b[idx], ww)
    se = float(np.std(boots, ddof=1))
    p = float(2 * _NORM.sf(abs(beta / se))) if se > 0 else np.nan
    return {"method": "Weighted mode", "nsnp": n, "b": beta, "se": se, "p": p}


def leave_one_out(bx, sx, by, sy):
    bx = np.asarray(bx, float); sx = np.asarray(sx, float)
    by = np.asarray(by, float); sy = np.asarray(sy, float)
    out = []
    full = ivw(bx, sx, by, sy)
    for i in range(len(bx)):
        m = np.ones(len(bx), bool); m[i] = False
        if m.sum() < 1:
            continue
        r = ivw(bx[m], sx[m], by[m], sy[m])
        out.append({"dropped": i, "b": r["b"], "se": r["se"], "p": r["p"],
                    "delta_b": r["b"] - full["b"]})
    return out, full


def mr_presso_like(bx, sx, by, sy, nsim=2000, seed=20261008, alpha=0.05):
    """Global test (observed Q vs simulated) + simple outlier-corrected IVW."""
    b, se = ratio(bx, sx, by, sy)
    Q, df, pQ = cochran_q(b, se)
    rng = np.random.default_rng(seed)
    # simulate under the null that all SNPs share one causal effect
    w = 1.0 / se ** 2
    beta0 = float(np.sum(w * b) / np.sum(w))
    sims = np.empty(nsim)
    for i in range(nsim):
        bb = rng.normal(beta0, se)
        ww = 1.0 / (se ** 2)
        bs = np.sum(ww * bb) / np.sum(ww)
        sims[i] = np.sum(ww * (bb - bs) ** 2)
    p_global = float((sims >= Q).mean())
    out = {"Q": Q, "Q_df": df, "Q_p_chi2": pQ, "global_p_sim": p_global,
           "mean_Q_sim": float(sims.mean())}
    if df >= 1 and np.isfinite(pQ):
        # outlier p-values: studentised residual of the ratio w.r.t. IVW
        resid = b - beta0
        t = np.abs(resid) / se
        p_out = 2 * stats.t.sf(t, df)
        out["outlier_p"] = p_out.tolist()
        keep = p_out > (alpha / len(b))
        out["n_outlier"] = int((~keep).sum())
        if keep.sum() >= 1 and (~keep).sum() > 0:
            r = ivw(bx[keep], sx[keep], by[keep], sy[keep])
            out["corrected_b"] = r["b"]; out["corrected_se"] = r["se"]
            out["corrected_p"] = r["p"]; out["corrected_nsnp"] = int(keep.sum())
    return out


def steiger(r2_x, r2_y):
    """Direction: exposure explains more variance than outcome."""
    return {"r2_exposure": r2_x, "r2_outcome": r2_y,
            "correct_direction": bool(r2_x > r2_y)}


def r2_from_z(beta, se, n):
    z = np.asarray(beta, float) / np.asarray(se, float)
    return float(np.sum(z ** 2) / n)
