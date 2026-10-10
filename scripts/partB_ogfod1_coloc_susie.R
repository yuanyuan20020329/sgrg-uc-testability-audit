# ============================================================
# Part B addendum — OGFOD1 locus coloc sensitivity (coloc.susie)
#
# Question: the Part C MR signal for OGFOD1 (Whole_Blood eQTL -> UC,
# IVW p=1.7e-7) — is it driven by the SAME variant as the OGFOD1 eQTL,
# or is it contaminated by the much stronger CES1P1 eQTL 631 kb away?
# Part C's coloc.abf (truncated SNP set, PP.H3=0.87) already hints at
# the latter; here we test it formally with multi-causal coloc.susie
# on both traits (OGFOD1 eQTL and CES1P1 eQTL vs the UC GWAS),
# an external LD panel (1kGP high-coverage, GRCh38, EUR n=633).
#
# Design follows the coloc-susie-external-ld skill:
#   (1) coloc.abf on the SAME SNP set  -> comparability benchmark
#   (2) coloc.susie main (p12=1e-5, matching Part C prior) + p12 grid
#   (3) L=1 limit                      -> does the panel bias single-
#       causal inference?  PP.H4 should approximate coloc.abf
#   (4) zero-noise synthetic positive control
#       z = z_peak * R[, peak]          -> must give exactly 1 CS
#       containing the true peak, else results are uninterpretable
#   (5) panel pathology: near-perfect correlations & MAF
#
# Traits:
#   d1 = GTEx v8 Whole_Blood cis-eQTL (normalised effect sizes,
#        sdY=1, N=670, type="quant")
#   d2 = UC GWAS (de Lange 2017, GCST004133; 12,366 cases /
#        33,609 controls; N=45,975, s=0.26895, type="cc")
# ============================================================
.libPaths(c("C:/Users/H-P/AppData/Local/R/win-library/4.6", .libPaths()))
suppressMessages({library(data.table); library(coloc); library(susieR)})

RES <- "D:/ld_ref_1000g/coloc_work"   # ASCII path (R fread can't read CJK paths)
REGIONS <- c("OGFOD1", "CES1P1")
EQTL_N <- 670
UC_N <- 45975
UC_S <- 12366 / UC_N
P12_MAIN <- 1e-5
P12_GRID <- c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4)

as_susie <- function(x) { if (!("susie" %in% class(x))) class(x) <- c("susie", class(x)); x }

mk_ds <- function(dt, R) {
  list(d1 = list(beta = dt$beta_eqtl, varbeta = dt$se_eqtl^2, snp = dt$snp,
                 position = dt$pos, type = "quant", N = EQTL_N, sdY = 1,
                 LD = R),
       d2 = list(beta = dt$beta_gwas, varbeta = dt$se_gwas^2, snp = dt$snp,
                 position = dt$pos, type = "cc", s = UC_S, N = UC_N,
                 LD = R))
}

run_set <- function(ds, opt) {
  a1 <- tryCatch(suppressWarnings(do.call(runsusie, c(list(d = ds$d1, suffix = 1), opt))),
                 error = function(e) NULL)
  a2 <- tryCatch(suppressWarnings(do.call(runsusie, c(list(d = ds$d2, suffix = 2), opt))),
                 error = function(e) NULL)
  list(s1 = a1, s2 = a2)
}

lead_of <- function(cs, z, snps) snps[cs[which.max(abs(z[cs]))]]

susie_sets <- list(
  main = list(maxit = 1000, repeat_until_convergence = TRUE, estimate_residual_variance = FALSE),
  erv  = list(maxit = 1000, repeat_until_convergence = TRUE, estimate_residual_variance = TRUE)
)

main_rows <- list(); cs_rows <- list(); pair_rows <- list()
abf_rows <- list(); sens_rows <- list(); set_rows <- list(); diag_rows <- list()

for (rid in REGIONS) {
  fp <- file.path(RES, sprintf("ogfod1_coloc_%s_dataset.tsv", rid))
  lp <- file.path(RES, sprintf("ogfod1_coloc_%s_ld.csv", rid))
  dt <- fread(fp)
  ld_dt <- fread(lp, header = TRUE)
  R <- as.matrix(ld_dt[, -1])
  rownames(R) <- colnames(R) <- as.character(ld_dt[[1]])
  stopifnot(identical(colnames(R), dt$snp), nrow(dt) == nrow(R))
  diag(R) <- 1
  cat(sprintf("\n=== %s | SNP=%d | min p(eQTL)=%.3g min p(GWAS)=%.3g ===\n",
              rid, nrow(dt), min(dt$p_eqtl), min(dt$p_gwas, na.rm = TRUE)))

  ds <- mk_ds(dt, R)
  z1 <- dt$beta_eqtl / dt$se_eqtl
  z2 <- dt$beta_gwas / dt$se_gwas

  # ---------- (1) coloc.abf on the same SNP set ----------
  ab <- tryCatch(suppressWarnings(coloc.abf(ds$d1, ds$d2, p12 = P12_MAIN)),
                 error = function(e) { cat("  abf ERR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(ab)) {
    abf_rows[[length(abf_rows) + 1]] <- data.table(
      region = rid, n_snp = nrow(dt),
      PP.H0 = ab$summary[["PP.H0.abf"]], PP.H1 = ab$summary[["PP.H1.abf"]],
      PP.H2 = ab$summary[["PP.H2.abf"]], PP.H3 = ab$summary[["PP.H3.abf"]],
      PP.H4 = ab$summary[["PP.H4.abf"]])
    tp <- as.data.table(ab$results)[which.max(SNP.PP.H4)]
    cat(sprintf("  [abf same-SNP-set p12=1e-5] PP.H4=%.4f PP.H3=%.4f top-SNP=%s (PP=%.3f)\n",
                ab$summary[["PP.H4.abf"]], ab$summary[["PP.H3.abf"]], tp$snp, tp$SNP.PP.H4))
  }

  # ---------- (2)+(3) coloc.susie under two settings ----------
  fits <- list()
  for (nm in names(susie_sets)) {
    ss <- run_set(ds, susie_sets[[nm]])
    if (is.null(ss$s1) || is.null(ss$s2)) { cat("  [", nm, "] susie failed\n"); next }
    n1 <- length(ss$s1$sets$cs); n2 <- length(ss$s2$sets$cs)
    r <- tryCatch(suppressWarnings(coloc.susie(as_susie(ss$s1), as_susie(ss$s2),
                                               p12 = P12_MAIN)),
                  error = function(e) NULL)
    ok <- !is.null(r) && !is.data.table(r)
    pp4 <- if (ok) max(as.data.table(r$summary)$PP.H4.abf) else NA_real_
    pp3 <- if (ok) { s <- as.data.table(r$summary); s[which.max(PP.H4.abf)]$PP.H3.abf } else NA_real_
    set_rows[[length(set_rows) + 1]] <- data.table(
      region = rid, setting = nm, n_snp = nrow(dt), n_cs1 = n1, n_cs2 = n2,
      PP.H4.best = pp4, PP.H3.best = pp3,
      lead1 = if (n1) paste(vapply(ss$s1$sets$cs, lead_of, "", z1, dt$snp), collapse = "|") else NA_character_,
      lead2 = if (n2) paste(vapply(ss$s2$sets$cs, lead_of, "", z2, dt$snp), collapse = "|") else NA_character_)
    cat(sprintf("  [susie %-4s] CS(eqtl)=%d CS(UC)=%d -> best PP.H4=%.4f PP.H3=%.4f\n",
                nm, n1, n2, ifelse(is.na(pp4), -1, pp4), ifelse(is.na(pp3), -1, pp3)))
    if (ok) fits[[nm]] <- list(sm = as.data.table(r$summary), s1 = ss$s1, s2 = ss$s2)
  }

  M <- fits[["main"]]
  if (!is.null(M)) {
    sm <- M$sm[order(-PP.H4.abf)]
    best <- sm[1]
    main_rows[[length(main_rows) + 1]] <- data.table(
      region = rid, n_snp = nrow(dt),
      n_cs_eqtl = length(M$s1$sets$cs), n_cs_uc = length(M$s2$sets$cs),
      n_pairs = nrow(sm), PP.H4.best = best$PP.H4.abf, PP.H3.best = best$PP.H3.abf,
      hit1 = best$hit1, hit2 = best$hit2)
    pair_rows[[length(pair_rows) + 1]] <-
      sm[, .(region = rid, idx1, idx2, nsnps,
             PP.H3 = PP.H3.abf, PP.H4 = PP.H4.abf, hit1, hit2)]
    for (tr in c("eqtl", "uc")) {
      obj <- if (tr == "eqtl") M$s1 else M$s2
      cs <- obj$sets$cs; z <- if (tr == "eqtl") z1 else z2
      if (length(cs)) for (i in seq_along(cs)) cs_rows[[length(cs_rows) + 1]] <- data.table(
        region = rid, trait = tr, cs = i, size = length(cs[[i]]),
        lead_snp = lead_of(cs[[i]], z, dt$snp),
        min_abs_corr = if (!is.null(obj$sets$purity)) obj$sets$purity$min.abs.corr[i] else NA_real_)
    }
    # p12 grid reusing the susie objects
    for (p12 in P12_GRID) {
      f2 <- tryCatch(suppressWarnings(coloc.susie(as_susie(M$s1), as_susie(M$s2), p12 = p12)),
                     error = function(e) NULL)
      ok2 <- !is.null(f2) && !is.data.table(f2)
      sens_rows[[length(sens_rows) + 1]] <- data.table(
        region = rid, p12 = p12,
        PP.H4.best = if (ok2) max(as.data.table(f2$summary)$PP.H4.abf) else NA_real_)
    }
  }

  # ---------- (3b) L=1 limit ----------
  l1 <- tryCatch({
    o <- list(L = 1, maxit = 1000, repeat_until_convergence = TRUE,
              estimate_residual_variance = FALSE)
    ss <- run_set(ds, o)
    if (is.null(ss$s1) || is.null(ss$s2)) NA else {
      r <- suppressWarnings(coloc.susie(as_susie(ss$s1), as_susie(ss$s2), p12 = P12_MAIN))
      if (is.data.table(r)) NA else max(as.data.table(r$summary)$PP.H4.abf)
    }
  }, error = function(e) NA_real_)
  cat(sprintf("  [L=1 limit] PP.H4=%.4f\n", ifelse(is.na(l1), -1, l1)))

  # ---------- (4) zero-noise synthetic positive control ----------
  # NOTE: the peak z is fixed at 8 (the real GWAS peak z is ~1.6, far too weak
  # for susie to detect — using it here was a design bug that produced 0 CS).
  # The control validates the code path, not the data.
  syn <- tryCatch({
    pk <- which.min(dt$p_gwas)
    zs <- 8 * R[, pk] * sign(z2[pk])
    syn_d <- list(beta = zs, varbeta = rep(1, length(zs)), snp = dt$snp,
                  position = dt$pos, type = "quant", N = UC_N, sdY = 1, LD = R)
    s <- suppressWarnings(runsusie(d = syn_d, suffix = "s",
                                   L = 10, maxit = 1000,
                                   repeat_until_convergence = TRUE,
                                   estimate_residual_variance = FALSE))
    ncs <- length(s$sets$cs)
    hit <- if (ncs) any(vapply(s$sets$cs, function(cs) dt$snp[pk] %in% dt$snp[cs], FALSE)) else FALSE
    list(n_cs = ncs, peak = dt$snp[pk], contains_peak = hit)
  }, error = function(e) list(n_cs = -1, peak = NA_character_, contains_peak = NA))
  cat(sprintf("  [synthetic z_peak=8] n_CS=%d, contains true peak (%s): %s\n",
              syn$n_cs, syn$peak, syn$contains_peak))

  # ---------- (5) panel pathology ----------
  ut <- upper.tri(R); ar <- abs(R[ut])
  ij <- which(ut, arr.ind = TRUE)
  maf <- pmin(pmin(dt$eaf_panel, 1 - dt$eaf_panel)[ij[, 1]],
              pmin(dt$eaf_panel, 1 - dt$eaf_panel)[ij[, 2]])
  hi <- which(ar > 0.98)
  diag_rows[[length(diag_rows) + 1]] <- data.table(
    region = rid, n_snp = nrow(dt), n_pairs_gt098 = length(hi),
    frac_maf_lt005 = if (length(hi)) mean(maf[hi] < 0.05) else NA_real_,
    min_maf_hipairs = if (length(hi)) min(maf[hi]) else NA_real_,
    median_abs_r = median(ar), L1_PP.H4 = l1,
    syn_n_cs = syn$n_cs, syn_contains_peak = syn$contains_peak)
}

fwrite(rbindlist(main_rows, fill = TRUE), file.path(RES, "ogfod1_coloc_susie_main.csv"))
fwrite(rbindlist(cs_rows, fill = TRUE),   file.path(RES, "ogfod1_coloc_susie_cs.csv"))
fwrite(rbindlist(pair_rows, fill = TRUE), file.path(RES, "ogfod1_coloc_susie_pairs.csv"))
fwrite(rbindlist(abf_rows, fill = TRUE),  file.path(RES, "ogfod1_coloc_abf_samesnpset.csv"))
fwrite(rbindlist(sens_rows, fill = TRUE), file.path(RES, "ogfod1_coloc_p12_sensitivity.csv"))
fwrite(rbindlist(set_rows, fill = TRUE),  file.path(RES, "ogfod1_coloc_settings.csv"))
fwrite(rbindlist(diag_rows, fill = TRUE), file.path(RES, "ogfod1_coloc_diagnostics.csv"))

cat("\n================ coloc.abf (same SNP set) ================\n")
print(rbindlist(abf_rows))
cat("\n================ coloc.susie main ================\n")
if (length(main_rows)) print(rbindlist(main_rows, fill = TRUE)) else cat("  (no CS pairs — the UC GWAS has no credible set in either region)\n")
cat("\n================ CS leads ================\n")
if (length(cs_rows)) print(rbindlist(cs_rows, fill = TRUE)) else cat("  (none)\n")
cat("\n================ diagnostics ================\n")
print(rbindlist(diag_rows, fill = TRUE))
cat("\n================ p12 sensitivity ================\n")
if (length(sens_rows)) print(dcast(rbindlist(sens_rows), region ~ p12, value.var = "PP.H4.best")) else cat("  (not estimable — coloc.susie requires >=1 CS in both traits)\n")
