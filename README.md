# Testability audit of public omics resources — code and key results

This repository contains the analysis code, the key result tables, and the
environment specifications for the manuscript

> **A five-rung testability audit for public omics resources**
> (working title; the manuscript source is in [`manuscript/`](manuscript/))

The project asks a question that sits *before* reproducibility: once a dataset
has been obtained and a pipeline has been reproduced, **could that dataset have
returned the answer it returned?** We formalise the question as an ordered
ladder of five zero-cost checks, apply it to a real published-looking claim
(stress-granule-related gene set × ulcerative colitis), and then validate the
ladder on two further independent published claims.

## Release notes

**v1.0.2** — no change to any analysis, any number, or any conclusion. The deposited
dataset payload's **documentation** was synchronised with the Data Descriptor that cites it:

1. **`G3_powered_control` is described as a preset design-eligibility gate, not a "power
   gate".** The v1.0.1 `datapackage.json` used the word *power*, which contradicted the
   semantics the manuscript explains. The field itself is unchanged; only its description is.
2. **The warnings now ship with the data.** `n_paired_donors` / `n_paired_same_site` are
   marked *candidate* counts that need manual verification, and the known control-vocabulary
   misclassification (`normal pressure hydrocephalus`) plus the three per-table UC-selection
   conventions are written into `README.md` and the affected field descriptions. Field names,
   table structure, the four gates, `USABLE`, every numeric value and the row order are
   unchanged. `scripts/test_documentation_sync.py` asserts exactly that by diffing the deposit
   against the frozen v1.0.1 payload file by file, and
   `scripts/audit_data_descriptor.py --pre-publish` runs it as a hard pre-publication gate.
3. **This code archive was re-released** so that it is the archive that built the v1.0.2
   payload. New in v1.0.2: `scripts/test_documentation_sync.py` and the updated
   `scripts/build_data_deposit.py` templates.

**v1.0.1** — no change to any analysis, any number, or any conclusion. Two things moved:

1. **The deposited dataset payload was normalised to English.** The per-row `reason` string
   of the four-gate ledger, the field glossary and `caveat` note of the summary JSON, one
   example string, and two schema descriptions were written in Chinese in v1.0.0 and are now
   English, so a reader of the deposit can see *why a resource failed* without a translation
   layer. Field names, table structure, the four gates, `USABLE`, every numeric value and the
   row order are unchanged. `scripts/test_language_normalisation.py` asserts exactly that by
   diffing the deposit against the frozen v1.0.0 payload file by file, and
   `scripts/audit_data_descriptor.py --pre-publish` runs it as a hard pre-publication gate.
2. **This code archive was re-released** so that the scripts in it are the scripts that
   produced the v1.0.1 data. New in v1.0.1: `resfeas_reason_i18n.py` (the bilingual strings,
   one source of truth), `i18n_normalise_reason.py` (applies them to already-generated
   results without re-running the plots), `test_language_normalisation.py` (the byte-level
   language-only claim) and `results/partD_resfeas_ledger_reason_i18n.tsv` (the v1.0.0 →
   v1.0.1 reason mapping, row by row).

v1.0.0 remains published and immutable, with its own version DOI.

## The audit ladder

| Rung | Question | Quantity to report | Cost |
|------|----------|--------------------|------|
| **L1** | Is the control arm what the depositor says it is? | `is_primary_data` composition of the control arm; share of control cells from case-free studies | one `obs` column |
| **L2** | Is the design variable separable from batch? | Cramér's V(batch, group); batch-only AUC | two metadata columns |
| **L2'** | Are the instruments independent? | effective number of independent instruments from the LD matrix | one reference panel |
| **L3** | Is the score confounded by technical load? | ρ(score, nUMI); R²(depth ~ disease); post-rarefaction effect | one integer column |
| **L4** | Is the null the right null? | empirical p from expression-stratified random gene sets | N permuted sets |
| **L5** | Could the design have detected the effect at all? | minimum attainable p = 2/2ⁿ; positive-control recovery | design arithmetic |

## Repository layout

```
scripts/      analysis code (Python 3.13 + R 4.6)
results/      key result tables backing every number in the manuscript
manuscript/   manuscript source (EN) and the CN appendix containing the
              number → source mapping tables
env/          pinned Python and R environments
zenodo/       manifest with SHA-256 checksums for the archived release
```

## What is *not* here

**No raw public data are redistributed.** Every primary dataset is publicly
available and is identified by accession in
[`DATA_ACCESS.md`](DATA_ACCESS.md), together with the script that downloads it.
Only derived, aggregated result tables are shipped.

## Reproducing

See [`REPRODUCE.md`](REPRODUCE.md) for the step-by-step run order, and
[`env/`](env/) for the exact interpreter and package versions.

## Data availability

All primary data are public. Accessions, download routes and the corresponding
retrieval scripts are listed in [`DATA_ACCESS.md`](DATA_ACCESS.md).

## Citation

See [`CITATION.cff`](CITATION.cff). The concept DOI for this code archive is
[10.5281/zenodo.23291758](https://doi.org/10.5281/zenodo.23291758), which resolves to the
newest version. A version DOI for this release is minted via Zenodo upon publication of the
GitHub release and is recorded in `CITATION.cff` and
[`zenodo/deposit_metadata.md`](zenodo/deposit_metadata.md).

## Licence

Code: MIT (see [`LICENSE`](LICENSE)). Derived result tables: CC BY 4.0.
