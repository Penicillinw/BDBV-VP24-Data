# Reproducing the analysis

## 1. Software environment

The analysis was run on Windows 10 (10.0.19045) with:

| component | version |
|---|---|
| Python | 3.14.6 |
| numpy | 2.5.2 |
| pandas | 3.0.5 |
| scipy | 1.18.0 |
| scikit-learn | 1.9.0 |
| h5py | 3.16.0 |
| anndata | 0.13.2 |
| matplotlib | 3.11.1 |
| seaborn | 0.13.2 |
| biopython | 1.88 |
| requests | 2.34.2 |
| PyMOL | 3.1.1 (structure figures only) |

Machine-readable copy: `07_analysis_code_and_environment/environment/`.

## 2. Random seeds

| seed | where it is used |
|---|---|
| `20260919` | recorded bootstrap seed in the environment record |
| `20260926` | cross-database census, per-cell-line and primary-immune resampling |
| `20260930` | stratified bootstrap for the Fig. 3c AUC intervals; composite-axis rank drift |

Permutation counts are stated with each result: 4,000 label shuffles for the cross-atlas test,
5,000 stratified bootstrap draws for the AUC intervals, 2,000 expression-matched random modules
and 2,000 label permutations for the nasal-epithelium panel, 200,000 Monte Carlo draws for the
GSE46599 contrasts where the permutation space could not be enumerated, and full enumeration of
all 126 label assignments for the paired AUC difference.

## 3. Running the code

The scripts were written to run inside the original project tree and read their inputs by fixed
relative paths (for example `ROOT/analysis/g6c_primary_immune_20260926/...`). They are deposited
here unchanged, with a standardised name.

To run them from this repository:

1. Download the `link_only` inputs you need. `00_index_and_metadata/LINK_ONLY_RAW.tsv` lists each
   one with its dataset and size; each accession's download route is in
   `00_index_and_metadata/LICENSE_AND_ACCESS.md`.
2. Recreate the expected input layout (or edit a copy of `ROOT` in the script — do not edit the
   deposited file if you intend to keep checksums valid).
3. Entry points, by analysis block:

| block | script |
|---|---|
| VP24 census and cross-database screen | `..._vp24-config_code_g1-crossdb-census...py`, `..._g1b-crossdb-screen...py`, `..._ncbi-bdbv-seq-survey...py` |
| 2026 outbreak reads | the `bdbv_2026` / `vp24` scripts in `02_primary_sequence_data` provenance |
| interface geometry | `..._vp24-interface_code_b1-independent-interface...py`, `..._t4-vp24-struct-seqmap...py` |
| IFN-tone score Δ | `..._ifn-tone_code_g9-delta-robustness...py`, `..._a1c-components-and-kpna...py` |
| cross-atlas agreement | `..._cross-atlas_code_b5-*...py`, `..._b5b-*...py` |
| importin-α axis | `..._importin-axis_code_g7-kpna-axis...py`, `..._g7-kpna-axis-null...py` |
| probe-aggregation grid | `..._proxy-test_code_a1-level-vs-fold...py`, `..._a1b-functional-form...py`, `..._g6e-analyse-gse46599...py` |
| tissue extension | `..._tissue-extension_code_a3-deconvolution...py`, `..._t4-g8-tissue-ifn-env...py` |
| nasal epithelium | `..._epithelium-panel_code_g6f-gse342661-probe...py` |
| coverage audit | `..._coverage_code_compartment-census...py`, `..._b6-compartment-coverage...py` |

4. **Before running any script**, point `ROOT` at your own checkout: the deposited scripts retain
   their original `ROOT` constant (a local absolute path) so that the code reads as it was run.
   Replace that constant, or wrap the call in a small runner that sets `ROOT` from an environment
   variable. All other delivered files have local absolute paths masked to `<local-path>`.

## 4. What can and cannot be reproduced from this repository alone

- **Fully reproducible from this repository**: every per-figure source table, the derived analysis
  tables, and the numbers they contain. Each panel table is accompanied by a note stating its origin.
- **Requires re-downloading `link_only` inputs**: anything that starts from a GEO series matrix, an
  H5 single-cell matrix, an HPA table or an NCBI/UniProt/ENA/BV-BRC query. These are registered with
  accession, size and checksum but are not redistributed.
- **Not reproducible in full**: the interactive structural figure rendering (PyMOL sessions are not
  part of the data package) and the DOCX/PDF assembly of the manuscript. Neither affects any number
  reported in the paper.

## 5. Known issues carried in the record

Five discrepancies between figure captions and the underlying tables were found while assembling
this package, and were **recomputed and resolved** in manuscript v15. The full assembly report is
deposited as
`01_source_data_by_figure/_shared/bdbv-vp24_metadata_doc_shared-w1-source-data-report_v01_20261002.md`
and the adjudication with replacement wording as
`01_source_data_by_figure/_shared/bdbv-vp24_metadata_doc_x-caption-decisions_v01_20261002.md`.
Every change is summarised in `CORRECTIONS.md`.

| # | discrepancy | resolution | evidence in this deposit |
|---|---|---|---|
| 1 | Fig. 1b caption distance for P83 to the second VP24 copy (55.0 Å) could not be reproduced under any of 66 definitions tried | clause deleted from the caption | `03_…/bdbv_vp24_model/…x1-p83-distance-sweep…tsv` |
| 2 | Fig. 1b contact count: 6 contacts ≤ 3.6 Å in the panel versus 7 contacts ≤ 4.0 Å in a superseded render input | caption now states the exact criterion and attribution; the panel's 6 contacts are confirmed | `03_…/pdb_4u2x_experimental/…x2-contacts-authoritative…tsv` |
| 3 | Fig. 3 caption quadrant numerals (I/III) did not match the frozen ledger labels (Q1/Q4) | caption now names the quadrants by definition | `01_source_data_by_figure/Fig3/…fig3a…tsv` |
| 4 | Fig. S2b p = 0.06 is a rounding of 0.0584 | caption now gives p = 0.058 | `01_source_data_by_figure/FigS2/…figs2b…tsv` |
| 5 | Fig. S2a's 0.50–0.77 span pools the two negative-group definitions | caption now says the span is taken over both | `01_source_data_by_figure/FigS2/…figs2a…tsv` |

One further correction was made: the Fig. S2d matrix cell for GSE342661 "Publicly released" was
refuted by its voucher and changed 0 → 1, so that dataset meets 4/7 rather than 3/7. The figure was
redrawn and the caption updated accordingly.

### Quantities that cannot be recomputed

41 values reported in the pilot-MD `openmm_results.json` (local RMSD, native-contact fraction,
minimised energy, void volume) cannot be recomputed from this deposit, because the workspace holds
only minimised coordinates and no trajectories. They are marked `not_recomputable` in
`03_primary_structure_and_interface/md_pilot/…x4-md-recomputed-metrics…tsv`. The manuscript reports
the pilot qualitatively and draws no number from it, so no published value depends on them.
