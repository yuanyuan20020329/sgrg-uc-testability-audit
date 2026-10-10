# Reproduction guide

The analysis is organised as four parts (A, A2, B, C, D) plus two independent
validation cases (A, B). Each part is a self-contained pipeline: retrieve →
analyse → write a result table into `results/`. Nothing here needs a cluster; the
largest single job is the streaming pass over ds4 (320,979 cells), which is
chunked by `h5ad` row block.

## 0. Environment

```bash
python -m pip install -r env/requirements.txt
Rscript -e 'install.packages(c("data.table","coloc","susieR"))'
```

Interpreter used throughout: CPython 3.13.12 (Windows x86-64); R 4.6.0.
All scripts resolve paths relative to the repository root and write into
`results/`. Raw data go into `data/raw/` (not shipped) and are created on first
run by the corresponding retrieval script.

> **Windows note.** R cannot read paths containing non-ASCII characters on the
> machine used here, so the colocalisation scripts stage files through an
> ASCII-only working directory. Adjust the `.libPaths()` / staging path at the
> top of `partB_ogfod1_coloc_susie.R` if you run elsewhere.

## 1. Gene set (shared by all parts)

```bash
python scripts/validate_sgrg_hgnc.py       # GO:0010494 ∪ GO:0035061 → HGNC filter
python scripts/partD_get_genesets.py       # Hallmark + random calibration sets
```

Produces `results/sgrg_gene_set_clean.json`, `results/sgrg_rejected.tsv`.

## 2. Part A — bulk discovery and validation

```bash
python scripts/partA_step1_sgrg_de.py              # GSE38713 differential expression
python scripts/partA_step2_sgrg_model.py           # panel classifier
python scripts/partA_step2b_panel_specificity.py   # vs 200 random gene sets  (→ L4)
python scripts/partA_step2c_global_structure.py    # PC1 / global axis          (→ L3)
python scripts/partA_step2d_batch_corrected.py     # batch-centring collapse    (→ L2)
python scripts/partA_step3_paired_validation.py    # within-patient check
python scripts/partA_step4_external_validation.py  # GSE11223                   (→ L2)
```

Key outputs: `partA_GSE38713_SGRG_summary.json`,
`partA_step2b_panel_specificity.json`, `partA_step2c_global_structure.json`,
`partA_step2d_batch_corrected.json`, `partA_step4_external_validation.json`.

## 3. Part A2 — cohort landscape (Gemma)

```bash
python scripts/partA2_probe.py           # enumerate cohorts + curation metadata
python scripts/partA2_analyse.py         # batch-artifact prevalence
python scripts/partA2_check_svd.py       # batch-only AUC (self-consistency)
python scripts/partA2_mesh_verify.py     # Europe PMC MeSH confirmation
```

> **Interpretation guard.** `partA2_check_svd.py` compares Gemma's own
> `batchEffect` annotation against Gemma's own SVD — this is a *self-consistency*
> check, **not** an independent validation. The audit口径 must use PC1–3 (the
> curation paper's own criterion), not PC1 alone.

## 4. Part B — UC ↔ MDD bridge

```bash
python scripts/partB_fetch_mdd_de.py       # Gemma differential result sets
python scripts/partB_bridge_analysis.py    # enrichment + sign concordance
python scripts/partB_genomewide_19738.py   # genome-wide DE, GSE19738
python scripts/partB_cohort_counts.py
python scripts/partB_ogfod1_coloc_align.py # variant-space alignment
Rscript scripts/partB_ogfod1_coloc_susie.R # coloc.susie on external LD panel
```

Part B is the **negative that passes all five rungs**; it is the control that
shows the ladder is not a machine for producing failures.

## 5. Part C — Mendelian randomisation

```bash
python scripts/partC_get_eqtl.py           # GTEx v8 cis-eQTLs
python scripts/partC_a_instruments.py      # instrument selection
python scripts/partC_c_mr.py               # IVW / median / mode / Egger
python scripts/partC_ld_audit_fetch.py     # 1000G LD via BGZF + HTTP Range
python scripts/partC_ld_audit_analyse.py   # effective independent instruments (→ L2')
python scripts/partC_e_neighbour.py        # neighbour-gene colocalisation
```

> **Interpretation guard.** Distance pruning is **not** LD clumping. Report the
> LD matrix and the effective number of independent instruments with every MR
> estimate. Where `coloc.susie` returns no credible-set pair, the correct
> statement is "one end has no signal", not "colocalisation refuted".

## 6. Part D — single-cell atlases

```bash
python scripts/dl_parallel.py              # fetch the four h5ad assets
python scripts/partD_technical_inventory.py  # ds1–ds4 technical axis inventory (→ L2/L3)
python scripts/partD_batch_audit.py          # batch–disease nesting            (→ L2)
python scripts/partD_analyse_stream.py       # per-dataset pseudobulk + scores
python scripts/partD_random_control.py       # random-set null                  (→ L4)
python scripts/partD_provenance_analysis.py  # control-source audit             (→ L1)
python scripts/partD_gse133060_fetch.py      # the powered design
python scripts/partD_gse133060_test.py       # paired test + calibration        (→ L5)
```

Key outputs: `partD_technical_inventory.tsv`, `partD_batch_audit_summary.json`,
`partD_ds4_prov_summary.json`, `partD_ds4_tearle_*_summary.json`,
`partD_gse133060_verdict.json`.

> **Interpretation guard.** Rarefaction is mandatory before any between-group
> score comparison if the score correlates with `nUMI`. Always carry the random
> calibration sets through the *same* rarefaction — otherwise a reviewer can
> argue that rarefaction erased the signal.

## 7. Part D — library-wide resource audit (metadata only, no h5ad)

```bash
python scripts/partD_resfeas_cxg.py        # enumerate CELLxGENE collections
python scripts/partD_resfeas_libwide.py    # 397 collections / 2,238 datasets
python scripts/partD_resfeas_geo.py        # 521 human IBD series
python scripts/partD_resfeas_ledger.py     # four-gate ledger
```

Outputs: `partD_resfeas_libwide_{collections,datasets}.tsv`,
`partD_resfeas_libwide_summary.json`, `partD_resfeas_verdict.json`.
The 397 collection-detail JSONs are cached, so `--no-fetch` re-runs in seconds.

## 8. Independent cases

```bash
# Case A — ferroptosis × HCC
python scripts/caseA_fetch_data.py
python scripts/caseA_step34.py
python scripts/caseA_L4_random.py
python scripts/caseA_GSE65372_paired.py
python scripts/caseA_GSE65372_pseudorep.py

# Case B — pyroptosis × PD
python scripts/caseB_detail_PD.py
python scripts/caseB_step34.py
python scripts/caseB_cxg_neuro_deep.py
```

Both cases were **frozen before analysis** in a registered table
(`results/caseA_论文登记表.md`, `results/caseB_论文登记表.md` — Chinese;
the English write-ups are `manuscript/caseA_report.md` and
`manuscript/caseB_report.md`).

## 9. Self-checks

Two scripts enforce the manuscript's internal consistency. They are the reason
every number in the English text can be traced:

```bash
python scripts/reconcile_numbers.py   # every number in the EN body must exist in the CN appendix §B
python scripts/audit_citations.py     # every defined reference must be cited; no orphans
```

`reconcile_numbers.py` is the machine version of the rule *no number enters the
English draft unless it appears in the CN appendix §B source table*.

## 10. Figure generation

```bash
python scripts/partAD_figures.py          # Parts A/D main figures
python scripts/partB_figures.py
python scripts/partD_figures.py
python scripts/make_fig1_ladder.py        # Fig. 1 (the ladder schematic)
```

Fig. 1 is rendered with matplotlib at 17.0 cm width, all glyphs ≥ 6 pt, Arial.
Figure-number registration lives in `results/partAD_figure_numbers.json`.

## Assembly

`scripts/build_archive.py` rebuilds this repository's `scripts/`, `results/`,
`manuscript/`, `env/` and `zenodo/manifest.json` from the working project:

```bash
python scripts/build_archive.py --check   # report only
python scripts/build_archive.py           # rebuild release/
```
