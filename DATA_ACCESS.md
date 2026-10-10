# Data access

**No primary data are redistributed in this repository.** Every dataset below is
public and is reachable from the accession given. The corresponding *retrieval
script* in `scripts/` downloads (or queries) it and writes into `data/raw/`,
which is not shipped.

Two channels deserve a note:

* **Gemma** (https://gemma.msl.ubc.ca) is used both to retrieve curated
  expression matrices and — separately — to retrieve cohort-level curation
  metadata (`batchEffect`, `batchEffectStatistics`, `batchConfound`,
  `geeq.publicQualityScore`). Gemma dataset identifiers are given as
  `Gemma:<id>` alongside the GEO accession.
* **CELLxGENE Discover** (https://cellxgene.cziscience.com) is used at two
  levels: collection-detail JSON (metadata only — `cell_count`,
  `primary_cell_count`, `is_primary_data`, disease labels) and, for the four
  analysed datasets, the `.h5ad` asset. The library-wide provenance scan
  (§2.7) uses **metadata only** and never downloads an `h5ad`.

---

## Part A — bulk transcriptomics (ulcerative colitis)

| Role | Accession | Description | Route | Script |
|---|---|---|---|---|
| Discovery | **GSE38713** / Gemma:5873 | 43 human colonic mucosal biopsies (15 UC active inflamed, 8 UC remission, 7 UC active non-inflamed, 13 controls) | GEO / Gemma REST | `partA_step1_sgrg_de.py` |
| Validation | **GSE11223** / Gemma:8385 | 61 UC inflamed vs 62 normal colonic epithelial biopsies | GEO / Gemma REST | `partA_step4_external_validation.py` |
| Panel specificity | — | 200 size-matched random gene sets drawn from the 18,630-gene measured pool | derived | `partA_step2b_panel_specificity.py` |
| External panel test | **GSE22619** | independent colonic biopsy cohort | GEO / Gemma REST | `partA_step2*` |

## Part A2 — cohort landscape and curation-metadata audit

| Role | Source | Description | Route | Script |
|---|---|---|---|---|
| Cohort landscape | **Gemma** `query=ulcerative colitis` | 160 public UC cohorts with machine-readable curation metadata | Gemma REST `/datasets?query=` | `partA2_harvest_meta.py` |
| MeSH confirmation | **Europe PMC** | PMID → MeSH disease attribution for the retrieved cohorts | Europe PMC REST | `partA2_mesh_verify.py` |
| Batch separability | **Gemma** SVD output | batch-only AUC per cohort (self-consistency check, not independent validation) | Gemma REST | `partA2_check_svd.py` |

## Part B — UC ↔ MDD shared-bridge test

| Role | Accession | Description | Route | Script |
|---|---|---|---|---|
| MDD discovery | **GSE19738** | 66 MDD vs 66 controls, genome-wide; brain | GEO / Gemma REST | `partB_fetch_mdd_de.py` |
| Brain replication | **GSE53987**, **GSE80655**, **GSE22619** | independent brain-region cohorts | GEO / Gemma REST | `partB_bridge_analysis.py` |
| MDD differential expression | **Gemma** differential analyses | curated DE result sets (no matrix download) | Gemma REST `/datasets/{id}/analyses/differential` | `partB_fetch_mdd_de.py` |
| Colocalisation | **GTEx v8** cis-eQTLs | OGFOD1 / CES1P1 cis regions; 1000G EUR LD panel | GTEx API v2 + 1000G FTP | `partB_ogfod1_coloc_*` |

## Part C — summary-statistic Mendelian randomisation

| Role | Accession / source | Description | Route | Script |
|---|---|---|---|---|
| Exposure | **GTEx v8** | cis-eQTLs for SGRG genes in whole blood, colon sigmoid, colon transverse, brain cortex | GTEx Portal API v2 `/association/singleTissueEqtlByLocation` | `partC_get_eqtl*.py` |
| Outcome | **GCST004133** | de Lange 2017 UC GWAS (12,366 cases / 33,609 controls), harmonised, GRCh38 | GWAS Catalog | `partC_extract_gwas.py` |
| LD reference | **1000 Genomes** high-coverage GRCh38 phased panel (2022-04-22 release, 3,202 samples) | LD matrices for instrument clumping and effective-instrument count | EBI 1000G FTP | `ld_fetch_1kgp38.py`, `ld_bgzf_tabix.py` |
| Colocalisation | GTEx + GCST004133 cis windows | `coloc.abf` and `coloc.susie` (R) | derived | `partB_ogfod1_coloc_susie.R` |

EUR sample panel: `20130606_g1k_3202_samples_ped_population.txt` (space-delimited;
population column index 6).

## Part D — single-cell atlases (ulcerative colitis mucosa)

| Label | CELLxGENE collection | Dataset ID | Cells (deposited) | Source publication |
|---|---|---|---|---|
| **ds1** | `7c7bd6c2-925b-4034-baab-620ef1b760e1` | `9bfecd44-e6b5-4072-9b64-cb7f0a117f88` | 46,700 | Garrido-Trigo et al. *Nat Commun* 2023;14:4506 |
| **ds2** | `f11cb29c-b546-4738-9bd8-66ea621a7bd5` | `9c235282-2b6f-4f74-8f16-7d3ac1b51371` (*Extended – Myeloid cells*) | 50,570 | Oliver et al. *Nature* 2024;635:699–707 |
| **ds3** | `f11cb29c-b546-4738-9bd8-66ea621a7bd5` | `763d1d88-b38d-4ab0-b00f-72df816b4b08` (*Extended – Large Intestine (adult/pediatric)*) | 96,675 | Oliver et al. *Nature* 2024;635:699–707 |
| **ds4** | `75e35754-bd21-4562-a4f8-39eaa4dfce36` | `56bbabcb-3968-4206-8e0b-4a7d6bda5195` | 320,979 | Tearle et al. *Nat Commun* 2026;17:10107 (raw data E-MTAB-14955, E-MTAB-16401) |

> **Deposited cell counts are not analysis cell counts.** ds2/ds3 are read from
> the datasets listed above; ds4's 320,979 deposited cells are the full dataset
> including re-used cells. Analysis-level cell counts (after donor × cell-type
> filters) are recorded in `results/partD_technical_inventory.tsv`.

Retrieval: `dl_h5ad.py`, `dl_parallel.py` (CELLxGENE asset URLs; multi-part
download). Reading: `partD_analyse_stream.py` (h5py chunked streaming for
>100k cells).

### Library-wide provenance scan (metadata only)

| Role | Source | Description | Script |
|---|---|---|---|
| Cell-level provenance | **CELLxGENE Discover** — all 397 collection-detail endpoints | `cell_count`, `primary_cell_count`, `is_primary_data`, disease labels; 2,238 datasets, 292.8M cells | `partD_resfeas_libwide.py` |
| GEO feasibility | **GEO / NCBI eutils** | 521 human IBD × gut-segment series; modality from `!Series_type` SOFT field | `partD_resfeas_geo.py` |

Cached collection JSON (397 files, `data/raw/cxg_coll_detail/`) is re-fetchable
with `--no-fetch` disabled; the scan itself runs in seconds from the cache.

## Part D — powered bulk design

| Role | Accession | Description | Route | Script |
|---|---|---|---|---|
| Powered test | **GSE133060** | 8 UC patients × (inflamed vs non-inflamed mucosa); Agilent two-colour GPL4133; 16 arrays, 19,749 genes after probe aggregation | GEO | `partD_gse133060_fetch.py`, `partD_gse133060_test.py` |

## Case A — ferroptosis × hepatocellular carcinoma

| Role | Accession | Description | Route |
|---|---|---|---|
| Target paper | Yi S et al. *J Oncol* 2023;2023:1624580 (PMC9981290) | reported 11-gene panel, AUC = 1.000 | Europe PMC (open access) |
| Target's discovery cohort | **GSE65372** | 39 HCC + 15 non-tumour (54 samples from 18 patients) | GEO |
| Target's validation cohort | **GSE84402** | 14 + 14 — *not downloaded* (data "available upon request" upstream) | — |
| Independent replication | **GSE36376** / Gemma:7673 | 433 samples (240 tumour, 193 adjacent non-tumour) | GEO / Gemma REST |
| Gene set | **WP_FERROPTOSIS** (MSigDB C2:CP v2023.1, n = 64) | 62/64 measurable in GSE36376 | MSigDB |

## Case B — pyroptosis × Parkinson's disease

| Role | Accession | Description | Route |
|---|---|---|---|
| Target paper | Wang L et al. *Sci Rep* 2024;14:28548 (PMC11574289) | two-gene score, AUC = 0.823 | Europe PMC |
| Bulk training | **GSE7621**, **GSE20164**, **GSE20163**, **GSE8397** | 100 samples (59 PD / 41 control) | GEO / Gemma REST |
| Bulk validation | **GSE49036** | 23 samples | GEO / Gemma REST |
| scRNA "validation" arm | **GSE140231**, **GSE157783**, **GSE204796** | 5 PD / 27 healthy (as used by the target paper) | GEO (identity confirmed by literature cross-check) |
| Gene set | derived by the target paper | 3,088 genes at p < 0.05 from DESeq2 — **no curated pyroptosis set cited** | — |

## Gene-set sources

| Set | Source | Used in |
|---|---|---|
| SGRG panel | **QuickGO** GO:0010494 ∪ GO:0035061 → 156 symbols → 104 HGNC-approved → 97 measured | Parts A, B, C, D |
| Hallmark inflammation (5 programmes) | **MSigDB** Hallmark (`HALLMARK_INFLAMMATORY_RESPONSE`, `HALLMARK_TNFA_SIGNALING_VIA_NFKB`, `HALLMARK_INTERFERON_GAMMA_RESPONSE`, `HALLMARK_HYPOXIA`, `HALLMARK_UNFOLDED_PROTEIN_RESPONSE`) | Part D positive controls |
| WP_FERROPTOSIS | **MSigDB** C2:CP (WikiPathways) | Case A |
| Random calibration sets | drawn in-house, expression-stratified, panel genes excluded | all rungs (L4) |

---

## Environment

See [`env/`](env/): `requirements.txt` (Python 3.13), `R_dependencies.txt`
(R 4.6), `R_sessionInfo.txt`, `R_packages_snapshot.csv`.

## Notes on reachability

Two hosts were unreachable from the machine on which the analysis was run and
are **not** required to reproduce any shipped table: Ensembl REST and UCSC.
The NCBI web domain was intermittently blocked; the eutils API and the GEO FTP
were used instead, and every gzip payload was verified by decompression rather
than by HTTP status (a blocked request returns a reCAPTCHA page with status
200).
