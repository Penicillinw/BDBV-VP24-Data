# Analysis data — Bundibugyo virus VP24 four-position configuration and cell-type interferon restriction

This repository holds the **data behind the figures and analyses** of the accompanying manuscript:
per-figure source data, primary sequence and structural evidence, host reference tables,
external dataset registrations, the derived analysis tables, and the analysis code that
generates every number reported in the paper.

It is intended to serve as the core pillar of the manuscript's Data availability statement.

- **Manuscript**: *Bundibugyo virus carries a genus-exclusive four-position VP24 configuration that orders cell-type interferon restriction*
- **Author**: Wenhao Huang, Nanning No.3 High School, Nanning, Guangxi, China
- **Repository version**: 1.1.0 (2026-10-02)
- **Aligned to**: manuscript **v16** (the v15→v16 change is figure-citation formatting only; no data changed)
- **Contents**: 506 standardised data files (523 files in the repository), ~49 MB

## What this repository contains

| # | Directory | Contents |
|---|---|---|
| 00 | `00_index_and_metadata/` | file manifest, name map, checksums, data dictionary, licence and access terms, reproduction notes |
| 01 | `01_source_data_by_figure/` | **the numbers plotted in every panel** of Fig. 1–3 and Fig. S1–S4, one table per panel, plus a note per panel |
| 02 | `02_primary_sequence_data/` | BDBV VP24 record census, six-species cross-database screen, 2026 outbreak genomes, reference sequences, alignments and residue-call tables |
| 03 | `03_primary_structure_and_interface/` | PDB 4U2X coordinates, interface contact and distance tables, BDBV VP24 homology model coordinates |
| 04 | `04_host_reference_and_atlas/` | Human Protein Atlas v25.1 derived cell-type tables, Cell Ontology mapping, independent-atlas agreement tables, compartment membership |
| 05 | `05_external_public_datasets/` | per-series derived results for every GEO series used, plus `INDEX.tsv` |
| 06 | `06_analysis_data/` | derived tables for the IFN-tone score, the importin-α axis, the proxy tests, the tissue extension and the coverage audit |
| 07 | `07_analysis_code_and_environment/` | the analysis scripts that generate the tables, and the environment/seed record |
| 08 | `08_out_of_scope_and_superseded/` | ledger of everything deliberately **not** included, with the reason |

> **Corrections.** Building this deposit surfaced five caption-versus-data disagreements and one
> figure cell contradicted by its own voucher. All were recomputed from source and fixed in the
> manuscript (v15) and in the tables here. Every change, with its evidence file, is listed in
> `00_index_and_metadata/CORRECTIONS.md` — including the one quantity that could not be recomputed
> and is therefore marked as such rather than estimated.


## Figure → source data

| panel | source data file |
|---|---|
| Fig. 1a | `01_source_data_by_figure/Fig1/bdbv-vp24_vp24-config_table_fig1a_v01_20261002.tsv` |
| Fig. 1b | `01_source_data_by_figure/Fig1/bdbv-vp24_vp24-config_table_fig1b_v01_20261002.tsv` |
| Fig. 2a | `01_source_data_by_figure/Fig2/bdbv-vp24_ifn-tone_table_fig2a_v01_20261002.tsv` |
| Fig. 2b | `01_source_data_by_figure/Fig2/bdbv-vp24_ifn-tone_table_fig2b_v01_20261002.tsv` |
| Fig. 2c | `01_source_data_by_figure/Fig2/bdbv-vp24_ifn-tone_table_fig2c_v01_20261002.tsv` |
| Fig. 2d | `01_source_data_by_figure/Fig2/bdbv-vp24_ifn-tone_table_fig2d_v01_20261002.tsv` |
| Fig. 3a | `01_source_data_by_figure/Fig3/bdbv-vp24_importin-axis_table_fig3a_v01_20261002.tsv` |
| Fig. 3b | `01_source_data_by_figure/Fig3/bdbv-vp24_importin-axis_table_fig3b_v01_20261002.tsv` |
| Fig. 3c | `01_source_data_by_figure/Fig3/bdbv-vp24_importin-axis_table_fig3c_v01_20261002.tsv` |
| Fig. S1a–d | `01_source_data_by_figure/FigS1/bdbv-vp24_ifn-tone_table_figs1{a,b,c,d}_v01_20261002.tsv` |
| Fig. S2a–d | `01_source_data_by_figure/FigS2/bdbv-vp24_proxy-test_table_figs2{a,b,c,d}_v01_20261002.tsv` |
| Fig. S3a,b | `01_source_data_by_figure/FigS3/bdbv-vp24_vp24-interface_table_figs3{a,b}_v01_20261002.tsv` |
| Fig. S4a–d | `01_source_data_by_figure/FigS4/bdbv-vp24_epithelium-panel_table_figs4{a,b,c,d}_v01_20261002.tsv` |

Each panel table has a companion `.md` note recording where the values came from and how they were
checked against the manuscript. Two files in `01_source_data_by_figure/_shared/` complete the set:

- `bdbv-vp24_metadata_table_shared-figure-source-map_v01_20261002.tsv` — machine-readable panel → file map.
- `bdbv-vp24_metadata_doc_shared-w1-source-data-report_v01_20261002.md` — the full assembly report,
  including the five caption-versus-table discrepancies that were found and left unreconciled.

## Data sources

### Primary sequence and structure

- **NCBI Protein** — *Bundibugyo virus*[Organism] AND VP24[All Fields], retrieved 2026-09-26 and re-read 2026-09-30.
- **NCBI Nucleotide** — 2026 outbreak complete genomes; references NC_002549.1 (EBOV), NC_014373.1 (BDBV 2007),
  KC545393.1 (BDBV 2012), FJ217161 (BDBV 2007 isolate).
- **UniProt** — Q05322 (EBOV VP24), B8XCN4 / R4QRB9 (BDBV VP24); **BV-BRC** and **ENA** for the cross-database screen.
- **Pathoplexus LAPIS** — real-time surveillance instance aggregating deposited genomes.
- **RCSB PDB** — entry 4U2X (VP24–KPNA5), public domain, redistributed here.

### Host reference

- **Human Protein Atlas v25.1** — single-cell-type RNA reference table (154 human cell types, nCPM) and the
  Deep Visual Proteomics cell-type table. *Redistribution status pending — see LICENSE_AND_ACCESS.md.*
- **Cell Ontology** — 18,914 terms, used for the cross-atlas label mapping.
- **CELLxGENE Discover** — platform on which 1,629 public human datasets were screened for coverage auditing.

### External datasets (derived results only)

| accession | what it is |
|---|---|
| GSE21158 | 10 human cancer cell lines, control vs IFN-α-2a, 24 h, triplicate |
| GSE46599 | 9 immortalised lines + primary macrophages + primary CD4 T cells, ± type-I IFN, 24 h |
| GSE327707 | 6 donors × 5 sorted populations, ± IFN-α, 48 h |
| GSE306664 | 5 donors × 4 populations, IFN-α/β/γ/λ1, 21 h, single cell |
| GSE342661 | primary human nasal epithelial ALI cultures, IFN-β1 and IFN-λ1–4, 6/24/72 h |
| GSE309699 | normal human keratinocytes, IFN pre-treatment then EBOV-ΔVP30 or rVSV/EBOV-GP |
| GSE100839 | ARPE-19 cells, mock vs EBOV, 24 h |
| GSE300073 / GSE298600 | human iPSC-derived intestinal organoids, EBOV and MARV, 1 and 3 days |
| GSE107549 | clinical IFN trial, pre-treatment IFN-stimulated expression vs later control |
| GSE114905 | five-filovirus infection series in Huh7 cells, viral mRNA counts |
| Normandin et al. 2023 | rhesus macaque natural-history series, 21 animals across 17 tissues |

Per-series file counts, sizes and the number of registered link-only raw files are in
`05_external_public_datasets/INDEX.tsv`.

## Reproducing the analysis

See `00_index_and_metadata/REPRODUCE.md` for the software versions, the entry point of each script,
the random seeds and the permutation counts. Two points matter in practice:

1. The analysis scripts were written to run inside the original project tree and read their inputs by
   fixed relative paths. To run them from this repository, point `ROOT` at a checkout and use the
   packaging wrapper described in `REPRODUCE.md`.
2. Inputs that are `link_only` must be re-downloaded from their accessions first; they are not in
   this repository by design.

## Provenance and naming

Every delivered copy is renamed to a machine-readable, self-describing form:

```
bdbv-vp24_<topic>_<layer>_<unit>_<version>_<YYYYMMDD>.<ext>
```

`00_index_and_metadata/NAME_MAP.tsv` maps every delivered name back to its original path in the
project tree, and `00_index_and_metadata/FILE_MANIFEST.tsv` adds, for each file, its category,
role, originating dataset, size and SHA-256. `00_index_and_metadata/CHECKSUMS.sha256` covers every
file in the repository.

Validation performed before release: naming-grammar check (0 violations), per-file checksum
coverage (100%), per-panel source-data completeness (23/23 panels), and a scan confirming that no
figure, caption, manuscript or third-party raw file is present.

## Licence and third-party terms

Author-created content (derived tables, analysis code, metadata, documentation) is released under
**CC BY 4.0** — see `LICENSE`. Third-party data are not redistributed and remain subject to their
original terms; see `00_index_and_metadata/LICENSE_AND_ACCESS.md`.

## How to cite

Please cite the manuscript and this repository. Machine-readable metadata is in `CITATION.cff`.

## Contact

Corresponding author: Wenhao Huang. Please open an issue in this repository for questions about
the data or the code.
