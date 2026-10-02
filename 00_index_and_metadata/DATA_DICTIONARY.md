# Data dictionary

## 1. Conventions used across the package

- **Encoding**: UTF-8, no BOM. **Delimiter**: tab for `.tsv`, comma for `.csv`.
- **Missing values**: left empty. The manuscript's rule is carried through the data: a value that
  could not be read or called is recorded as *not determinable* and stays in the denominator; it is
  never silently converted into a negative. Where a table distinguishes the two, the column name
  says so (`cant_call`, `unresolved`, `n_unresolved`, `not_determinable`).
- **Directional labels**: `delta` / `Δ` is the frozen baseline IFN-tone score (equal-weight mean of
  three z-scored components; see §3.1). `KPNA_axis` is the equal-weight mean z-score of KPNA1,
  KPNA5 and KPNA6. Higher = more of the quantity named.
- **Ranks are 1-based and descending** (rank 1 = highest value) unless a table says otherwise.
- **`AUC`** values are in-sample areas under the ROC curve for "restricted" vs "non-restricted"
  units; 0.5 is chance.
- **`p_perm`** is a permutation p value; the number of shuffles is in `n_perm` / `p_perm_n` /
  `n_permutations` where present, and the seed is in `REPRODUCE.md`.
- **Cell-type names** follow the Human Protein Atlas v25.1 single-cell-type reference labels.

## 2. `FILE_MANIFEST.tsv`

| column | meaning |
|---|---|
| `file_id` | stable identifier within this release (`F0001`…) |
| `standardized_name` | delivered file name |
| `package_path` | path inside this repository |
| `original_path` | path in the author's project tree, before standardisation |
| `category_l1` / `category_l2` | one of the nine top-level categories / its subfolder |
| `role` | the file's layer: `table` (242), `code` (146), `doc` (24), `calls` (6), `meta` (1) |
| `topic` | manuscript topic the file belongs to (`vp24-config`, `ifn-tone`, `importin-axis`, `cross-atlas`, `proxy-test`, `filovirus-anchor`, `tissue-extension`, `epithelium-panel`, `spleen`, `coverage`, `md-pilot`, `metadata`) |
| `unit` | the dataset, figure panel or entity the file is about |
| `size_bytes` | file size |
| `sha256` | SHA-256 of the delivered file |

Related files: `NAME_MAP.tsv` (original → delivered name), `CHECKSUMS.sha256` (SHA-256 for every
tracked file except itself),
`LINK_ONLY_RAW.tsv` (third-party raw files that are registered but not redistributed).

## 3. `01_source_data_by_figure/`

One table per figure panel: **the numbers that are plotted**. A companion `.md` note per panel
records the origin of each value and how it was checked against the manuscript.

### 3.1 Fig. 1 — the four-position VP24 configuration

`..._table_fig1a_...tsv` (43 rows)

| column | meaning |
|---|---|
| `species`, `species_label` | orthoebolavirus species |
| `sequences_assessed`, `n_partial_skipped` | denominator construction: fragments outside the four-position window are excluded from the denominator, never scored as negative |
| `observation` | `residue_call` / `motif` / `coverage` |
| `position` | EBOV-numbered VP24 position (83, 135, 140, 141, 184, 217) |
| `residue`, `n_records`, `is_major_call` | the call and how many records carry it |
| `is_bdbv_defining_residue`, `is_bombali_call` | figure colouring flags |

`..._table_fig1b_...tsv` (10 rows) — per-contact geometry for the VP24–KPNA5 interface model:
`viral_residue`/`viral_atom` → `partner_residue`/`partner_atom`, `distance_A`, `structure`, and
`caption_value_A`/`matches_caption` for the figure-caption distances.

### 3.2 Fig. 2 — the cell-type IFN-tone axis Δ

| file | rows | key columns |
|---|---|---|
| `..._fig2a_...tsv` | 154 | `rank`, `cell_type`, `delta_restriction`, `set_membership` (Amplification / Sanctuary), `coverage_limited` |
| `..._fig2b_...tsv` | 17 | `IFN_I_capacity`, `IFN_III_capacity`, `ISG_priming` — the three Δ components |
| `..._fig2c_...tsv` | 26 | `compartment`, `n_members`, `median_delta`, `min_delta`, `max_delta`, `no_member_in_atlas`, `coverage_limited_marker` |
| `..._fig2d_...tsv` | 162 | `atlas`, `hpa_cell_type`, `cl_id`, `n_cells`, `hpa_delta`, `atlas_delta_all_labels`, `hpa_delta_matched_only`, `atlas_delta_matched_only` — one row per matched pair |

### 3.3 Fig. 3 — the importin-α axis

| file | rows | key columns |
|---|---|---|
| `..._fig3a_...tsv` | 154 | `delta_restriction`, `delta_rank`, `KPNA_axis`, `quadrant`, `KPNA1_z`/`KPNA5_z`/`KPNA6_z`, median-split columns |
| `..._fig3b_...tsv` | 26 | `n_Q1_highD_highK`…`n_Q4_lowD_lowK`, `mean_delta`, `p_exact_multinomial`, `q_bh`, `marker` |
| `..._fig3c_...tsv` | 18 | `negative_class_definition`, `probe_rule`, `candidate`, `n_positive`, `n_negative`, `auc_insample`, `ci_lo2.5`, `ci_hi97.5` |

### 3.4 Fig. S1–S4

| file | rows | key columns |
|---|---|---|
| `figs1a` / `figs1b` | 11 each | `ISG_baseline_log2` vs `ISG_IFN_level_log2` / `ISG_induction_log2` |
| `figs1c` | 11 | `cohort`, `unit`, `cell_type`, `baseline_isg_log2`, `additional_response_log2` |
| `figs1d` | 4 | donor-to-donor spread before and after IFN |
| `figs2a` | 42 | `probe_rule`, `negative_group_mode`, `quantity`, `auc`, `p_exact` |
| `figs2b` | 24 | exact null distribution of the paired AUC difference over 126 label assignments |
| `figs2c` | 15 | viral mRNA reads per virus per day (GSE114905) |
| `figs2d` | 35 | which dataset satisfies each of the seven requirements of the decisive test |
| `figs3a` | 29 | per-residue heavy-atom contacts in PDB 4U2X, with the six-species residues side by side |
| `figs3b` | 8 | informative records carrying the BDBV motif, per source database |
| `figs4a` / `figs4d` | 36 / 72 | module level per treatment × route × time (GSE342661), ISG17 and control modules |
| `figs4b` | 30 | induction above each sample's own paired control |
| `figs4c` | 10 | IFN-β1 minus type III, exact sign-flip test, 95% CI |

`..._shared/...figure-source-map...tsv` is the machine-readable panel → file map (23 rows) with
`source_origin`, `n_rows`, `key_columns`, `matches_manuscript` and `notes`.

## 4. `02_primary_sequence_data/`

- `call_tables/` — `pan_species_vp24_positions.tsv` (per-record calls at the four positions across
  the six species), `bdbv_full_length_breakdown.tsv`, `vp24_census_summary.json`,
  `routeG/routeN/routeP` records and the cross-check summaries.
- `genus_vp24_screen/` — `g1_summary.json` (NCBI/UniProt/BV-BRC/ENA census) and
  `g1b_summary.json` (exhaustive cross-database screen, with the per-species coverage denominators).
- `outbreak_2026_genomes/` — 2026 outbreak genome records, VP24 variant and QC tables, the
  reference CDS and the survey summary.
- `alignments/` — VP24 Clustal W alignment and the VP24 FASTA set. **Only VP24 was retained**;
  alignments for other genes in the project tree belong to analyses not reported here.

## 5. `03_primary_structure_and_interface/`

- `4U2X.pdb` — experimental VP24–KPNA5 complex (public domain, redistributed).
- `interface_contacts.tsv` — inter-chain heavy-atom contacts < 4.5 Å, per residue and per copy.
- `vp24_interface_residues.tsv`, `interface_species.tsv`, `chain_summary.tsv`, `summary.json`.
- `bdbv_vp24_model/` — homology-model coordinates and metadata, plus the two template-transfer models.

## 6. `04_host_reference_and_atlas/`

- `hpa_v25_1/` — HPA single-cell-type tables (154 cell types, nCPM), the Δ recompute and its
  per-component drop-out variants, and the DVP mRNA–protein paired tables.
- `independent_atlases/` — Cell Ontology mapping, label coverage, agreement tables, and the
  matched-pair tables used for Fig. 2d. Regenerate with `b5*`/`b5b*` scripts.
- `compartment_membership/` — compartment member lists, quadrant composition and coverage audit.

## 7. `05_external_public_datasets/`

One folder per accession (see `INDEX.tsv`). Each contains the **derived** tables produced from that
series — never the raw series files. Column names follow the conventions in §1; per-series design
detail (treatment arms, time points, donors) is in the manuscript's Supplementary Methods.

## 8. `06_analysis_data/`

| folder | contents |
|---|---|
| `ifn_tone_score/` | Δ robustness variants (`hpa_delta_variants.tsv`), atlas probe and matched pairs |
| `importin_alpha_axis/` | `kpna_axis.tsv` (154 cell types with all three gene z-scores), rank-drift tables, the AUC table, and the protein-layer cross-check |
| `proxy_tests/` | the probe-aggregation contrast grid (`module × group_rule × theta_deg → auc`), the functional-form fits and the level-vs-fold comparison |
| `coverage_audit/` | CELLxGENE Discover scan results, eye and ureter atlas contrast tables, reference coverage scans |
| `tissue_extension/` | NHP multi-organ deconvolution, composition, sensitivity and fixed-effects tables |
| `filovirus_anchor/` | GSE309699 module levels and the independent re-derivation |
| `spleen`, `epithelium_panel`, `sequence_census`, `md_pilot/` | as named |

## 9. `07_analysis_code_and_environment/`

`analysis/` holds the scripts that generate the tables; the file name carries the topic and the
original script stem. `environment/environment-versions...json` records the interpreter and package
versions and the bootstrap seed. See `REPRODUCE.md`.
