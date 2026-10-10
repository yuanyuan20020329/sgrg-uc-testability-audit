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

See [`CITATION.cff`](CITATION.cff). A DOI for the archived release is minted via
Zenodo upon publication.

## Licence

Code: MIT (see [`LICENSE`](LICENSE)). Derived result tables: CC BY 4.0.
