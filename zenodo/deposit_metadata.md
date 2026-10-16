# Zenodo deposit metadata (paste-ready)

## Title

A five-rung testability audit for public omics resources: code and key results

## Resource type

Software (with supplementary dataset tables) — *or* "Dataset" if you prefer to
lead with the result tables

## Version

v1.0.2

## Publication date

2026-10-16

## Creators

* Yuan, Yike — ORCID 0009-0009-8391-4134 — Hubei University of Chinese Medicine, Wuhan, Hubei, China — yyk20020329@126.com *(first author)*
* Chen, Hongci — ORCID 0009-0000-3335-8998 — Hubei Provincial Hospital of Traditional Chinese Medicine, Wuhan, Hubei, China — chenhongci@hbhtcm.com *(corresponding author)*

## Description

Analysis code and key result tables for an ordered five-rung audit of whether
public omics resources can support the claims made from them. The rungs are
ordered by the cost of the check and each has an explicit failure signature:

1. **Control provenance** — is the control arm what the depositor says it is?
   (`is_primary_data` composition; share of control cells from case-free studies)
2. **Batch–disease nesting** — is the design variable separable from batch?
   (Cramér's V, batch-only AUC)
2′. **Instrument independence** — effective number of independent instruments
   from the LD matrix
3. **Technical-load confounding** — is the score a function of sequencing depth?
   (ρ(score, nUMI); R²(depth ~ disease); post-rarefaction effect)
4. **Null calibration** — is the null valid? (empirical p from
   expression-stratified random gene sets)
5. **Power floor** — could the design have detected a known-positive effect?
   (minimum attainable p = 2/2ⁿ; positive-control recovery)

The ladder is applied to a stress-granule-related gene set × ulcerative colitis
claim across bulk transcriptomes, summary-statistic Mendelian randomisation and
four single-cell atlases, and then validated on two independent published
claims: ferroptosis × hepatocellular carcinoma (reproducible but non-specific)
and pyroptosis × Parkinson's disease (not recoverable from public curated data
under the published fixed weights). A whole-library provenance scan covers all
397 CELLxGENE Discover collections (2,238 datasets, 292.8M cells).

**No primary data are redistributed.** Every dataset is public and identified by
accession in `DATA_ACCESS.md`, with the retrieval script that fetches it.

## Keywords

testability; reproducibility; data provenance; batch confounding;
pseudoreplication; single-cell genomics; Mendelian randomisation;
ulcerative colitis; stress granules; data re-use

## Licence

* Code: MIT
* Result tables and manuscript: CC BY 4.0

## Related identifiers

* This archive (v1.0.2): 10.5281/zenodo.23412772
* Concept DOI (resolves to the newest version): 10.5281/zenodo.23291758
* Previous version (v1.0.1): 10.5281/zenodo.23298914 (published and immutable)
* Previous version (v1.0.0): 10.5281/zenodo.23291759 (published and immutable)
* Manuscript: <DOI / preprint URL — fill in>
* GitHub repository: https://github.com/yuanyuan20020329/sgrg-uc-testability-audit
* Target paper, case A: 10.1155/2023/1624580
* Target paper, case B: 10.1038/s41598-024-80185-9
* Primary data: see `DATA_ACCESS.md` (GEO, CELLxGENE Discover, GTEx,
  1000 Genomes, GWAS Catalog)

## Notes for the depositor

* Zenodo's GitHub integration archives a tagged release automatically; prefer
  that over a manual upload so the DOI tracks the tag.
* ✅ Done (v1.0.0, 2026-10-11): version DOI `10.5281/zenodo.23291759` backfilled into
  `CITATION.cff`, `README.md` and the manuscript's *Data and code availability* statement.
* ✅ Done (v1.0.1, 2026-10-11): version DOI `10.5281/zenodo.23298914` backfilled into
  `resource_paper/code_doi.txt`, `CITATION.cff`, `README.md`, this file and the manuscript's
  *Data and code availability* statement. The data payload cites the **v1.0.1** code DOI,
  not the v1.0.0 one — `scripts/audit_data_descriptor.py --pre-publish` fails otherwise.
* ✅ Done (v1.0.2, 2026-10-16): version DOI `10.5281/zenodo.23412772` backfilled into
  `resource_paper/code_doi.txt`, `CITATION.cff`, `README.md` and this file. The v1.0.2 data
  payload cites the **v1.0.2** code DOI. The manuscripts are updated in a later pass —
  this release is data-dictionary synchronisation only.
