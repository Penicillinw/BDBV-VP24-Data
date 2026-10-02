# W1 源数据报告 — v13 七图数值源数据重建

本轮把 v13 主稿与补充稿中 7 幅图（含子面板，共 23 个 panel）的数值源数据重建为 tidy 数表。下面每一段给出：来源产物、生成方式、与 v13 正文/图注的核对结果。所有表格只读既有产物，产物目录为 `_coord/_round_W_20261002/analysis/`。

**统一口径**：只收数表，不收成图、绘图脚本和图注；每个数字都能指回一个工作区文件，或由这些文件按 `build_source_tables.py` 的对应函数复算。图注文字只用于交叉核对，未作为源数据来源（唯一例外是 Fig. S2d 这一需求矩阵，它本身不是实验测量，其取值来自 S 轮冻结的图脚本网格）。

## 复现

```powershell
cd "<local-path>"
python build_source_tables.py        # 重建 23 张数表 + 映射表 + 说明
python _w_validate.py                # pandas 读入、面板齐全性与抽查
```

## 逐图说明

### Fig. 1a — `Fig1a_source.tsv`（43 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::fig1a()`，只读既有产物，不做数值改动；关键列 `species;position;residue;n_records;observation`。
- **与 v13 核对**：Retrieved records 3,522, fragments excluded 32, sequences assessed 3,490, callable at all four positions 3,471, non-Bundibugyo 3,440; all four figures are identical to the v13 Fig. 1 caption. Bundibugyo is the only species with SQHA (31 of 31).
- **说明**：Rows with residue X or - are retained and flagged is_major_call=False so that an uncalled position is never read as a negative.

### Fig. 1b — `Fig1b_source.tsv`（10 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::fig1b()`，只读既有产物，不做数值改动；关键列 `observation;distance_A;structure`。
- **与 v13 核对**：The six inter-chain polar contacts are 2.73-3.54 A, all below the 3.6 A stated in the caption. P83 to the nearest KPNA5 atom is 11.66 A in the model (caption 11.7) and 11.36 A in 4U2X (caption 11.4).
- **说明**：The caption's 55.0 A 'second VP24 copy' distance could not be reproduced under any natural definition; see the report's inconsistency list.

### Fig. 2a — `Fig2a_source.tsv`（154 行）

- **来源**：`<local-path> <local-path>
- **方法**：`build_source_tables.py::fig2()`，只读既有产物，不做数值改动；关键列 `rank;cell_type;delta_restriction;set_membership`。
- **与 v13 核对**：All 154 ranks and Delta values are identical to the round-U ledger; 12 amplification, 5 sanctuary, 137 other members.
- **说明**：Colour names differ between the round-S caption (vermillion/sky blue) and the v13 caption (salmon/periwinkle); no numeric change.

### Fig. 2b — `Fig2b_source.tsv`（17 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::fig2()`，只读既有产物，不做数值改动；关键列 `cell_type;IFN_I_capacity;IFN_III_capacity;ISG_priming`。
- **与 v13 核对**：17 representative cell types; each component is a z score across 154 cell types, as the caption states.
- **说明**：17 representative cell types, panel-b display order

### Fig. 2c — `Fig2c_source.tsv`（26 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::fig2()`，只读既有产物，不做数值改动；关键列 `compartment;median_delta;min_delta;max_delta;n_members`。
- **与 v13 核对**：26 compartments ranked by member median; the spleen carries 0 members (grey x) and urothelium / ventricular system carry the coverage-limited asterisk, as the caption states.
- **说明**：26 compartments (25 with members + spleen as no-member); ranked by member median

### Fig. 2d — `Fig2d_source.tsv`（162 行）

- **来源**：`<local-path> <local-path>
- **方法**：`build_source_tables.py::fig2()`，只读既有产物，不做数值改动；关键列 `atlas;hpa_cell_type;hpa_delta;atlas_delta_all_labels`。
- **与 v13 核对**：162 matched pairs over 18 atlases; hpa_delta is byte-identical to the round-U ledger. Caption statistics rho 0.519, p 0.0002, matched-only 0.452 are carried by the frozen analysis, not restated here.
- **说明**：The plotted y in panel d is atlas_delta_all_labels; both standardisations are provided as columns.

### Fig. 3a — `Fig3a_source.tsv`（154 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::fig3()`，只读既有产物，不做数值改动；关键列 `cell_type;delta_restriction;KPNA_axis;quadrant`。
- **与 v13 核对**：154 cell types; quadrant sizes 30/47/47/30 and the two axis medians (Delta 0.067159, KPNA 0.088719) match the g7 summary and the v13 caption; the six named quadrant extremes are carried per cell.
- **说明**：154 cell types; quadrant sizes {'Q3_lowD_highK': 47, 'Q4_lowD_lowK': 30, 'Q1_highD_highK': 30, 'Q2_highD_lowK': 47}; medians 0.067159/0.088719; 6 named quadrant extremes

### Fig. 3b — `Fig3b_source.tsv`（26 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::fig3()`，只读既有产物，不做数值改动；关键列 `compartment;n_members;n_Q1..n_Q4;mean_delta;q_bh;marker`。
- **与 v13 核对**：26 compartments with Q1-Q4 stacks, member mean Delta, exact multinomial p and BH q; asterisks in the panel follow q < 0.05.
- **说明**：26 compartments; exact multinomial GOF vs the atlas-wide 30/47/47/30 split, BH over the 25 with members

### Fig. 3c — `Fig3c_source.tsv`（18 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::fig3()`，只读既有产物，不做数值改动；关键列 `candidate;auc_insample;ci_lo2.5;ci_hi97.5`。
- **与 v13 核对**：Delta 0.867 (95% CI 0.60-1.00), importin-alpha mean(z) 0.467 (0.10-0.83) and z(mean) 0.50 (0.13-0.87) reproduce the v13 caption; intervals are 2.5th-97.5th percentiles of 5,000 stratified bootstrap draws, seed 20260930.
- **说明**：The reported AUC is the in-sample value; the def1 strict-n9 variant of the same table is not part of the plotted panel.

### Supplementary Fig. S1a — `FigS1a_source.tsv`（11 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::figs1()`，只读既有产物，不做数值改动；关键列 `cell_type;ISG_baseline_log2;ISG_IFN_level_log2`。
- **与 v13 核对**：Nine immortalised lines give rho +0.883 (exact p 0.0031), matching the caption's +0.88 / 0.003; the two primary cell types are shown but not entered.
- **说明**：9 immortalised lines (rho +0.88, p 0.003) plus 2 primary cell types shown but not entered

### Supplementary Fig. S1b — `FigS1b_source.tsv`（11 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::figs1()`，只读既有产物，不做数值改动；关键列 `cell_type;ISG_baseline_log2;ISG_induction_log2`。
- **与 v13 核对**：Baseline against own induction gives rho -0.083 (p 0.8432) and the housekeeping control -0.367, matching the caption.
- **说明**：same 11 cell types; housekeeping control rho -0.37 on the same axis

### Supplementary Fig. S1c — `FigS1c_source.tsv`（11 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::figs1()`，只读既有产物，不做数值改动；关键列 `cohort;unit;baseline_isg_log2;additional_response_log2`。
- **与 v13 核对**：CD14+ monocytes give rho -1.00 in both cohorts (exact p 0.003 and 0.017); the pooled five-population value is -0.954 (caption -0.95).
- **说明**：CD14+ monocytes rho = -1.00 in both cohorts (exact p 0.003 and 0.017); pooled five-population rho = -0.95

### Supplementary Fig. S1d — `FigS1d_source.tsv`（4 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::figs1()`，只读既有产物，不做数值改动；关键列 `cohort;phase;spread_across_donors_log2`。
- **与 v13 核对**：Spread falls from 4.707 to 0.779 log2 (GSE327707) and from 2.413 to 0.363 (GSE306664), the convergence the caption describes.
- **说明**：donor spread before vs after IFN: 4.707 -> 0.779 (GSE327707) and 2.413 -> 0.363 (GSE306664)

### Supplementary Fig. S2a — `FigS2a_source.tsv`（42 行）

- **来源**：`<local-path> GEO in <local-path>
- **方法**：`build_source_tables.py::figs2()`，只读既有产物，不做数值改动；关键列 `probe_rule;negative_group_mode;quantity;auc;p_exact`。
- **与 v13 核对**：Level AUC is 0.90 in mode A and 0.8333 in mode B under every one of the seven probe-collapsing rules; fold induction spans 0.50-0.77 and baseline tone 0.67-0.80 across the two modes.
- **说明**：The 0.50-0.77 fold span pools modes A and B; within mode A alone the span is 0.50-0.75.

### Supplementary Fig. S2b — `FigS2b_source.tsv`（24 行）

- **来源**：`<local-path> (exact_perm)`；`raw GEO GSE46599`
- **方法**：`build_source_tables.py::figs2()`，只读既有产物，不做数值改动；关键列 `negative_group_mode;diff_value;n_assignments_at_value;proportion`。
- **与 v13 核对**：All 126 mode-A and 462 mode-B label assignments were enumerated; the observed dAUC is +0.100 with p 0.222 (mode A) and +0.133 with p 0.058 (mode B), matching the caption's +0.10 / 0.22 and +0.13 / 0.06.
- **说明**：exact paired-difference null over all 126 (mode A) and 462 (mode B) label assignments; observed dAUC +0.10 p=0.22 and +0.13 p=0.06

### Supplementary Fig. S2c — `FigS2c_source.tsv`（15 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::figs2()`，只读既有产物，不做数值改动；关键列 `quantity;virus;day;value`。
- **与 v13 核对**：BDBV/EBOV ratios 0.0781, 0.0819 and 0.3370 match the caption's 0.08, 0.08 and 0.34; per-virus per-day totals are the deposited sums.
- **说明**：One library per virus per day, so this is a consistency check and not a host-response measurement.

### Supplementary Fig. S2d — `FigS2d_source.tsv`（35 行）

- **来源**：`<local-path> in figures_captions_revised.md`
- **方法**：`build_source_tables.py::figs2()`，只读既有产物，不做数值改动；关键列 `dataset;requirement;satisfied`。
- **与 v13 核对**：Counts 4/7, 5/7, 5/7, 4/7 and 3/7 requirements met; no dataset satisfies all seven.
- **说明**：This panel is a documentation matrix, not an experimental measurement; its cells are the binary requirement grid frozen in the round-S figure script.

### Supplementary Fig. S3a — `FigS3a_source.tsv`（29 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::figs3()`，只读既有产物，不做数值改动；关键列 `ebov_full_pos;contacts_heavy_atom;bdbv_differs;plotted`。
- **与 v13 核对**：29 interface residues with heavy-atom contact sums over the three cognate copies; Q184R has the second-largest footprint (78 atoms) and position 83 is flagged as the no-contact x.
- **说明**：The panel plots 15 of the 29 residues (plus position 83); the 14 omitted residues are retained here with plotted=false.

### Supplementary Fig. S3b — `FigS3b_source.tsv`（8 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::figs3()`，只读既有产物，不做数值改动；关键列 `source_db;species_group;n_with_all_four_bdbv_type;n_informative_records_examined`。
- **与 v13 核对**：BDBV 1/1, 2/3, 9/9 and 31/31; all other orthoebolaviruses 0/39, 0/62, 0/99 and 0/3459; identical to the v13 caption.
- **说明**：The BV-BRC denominator is three retrievable BDBV sequences because one of them could not be read at all four positions.

### Supplementary Fig. S4a — `FigS4a_source.tsv`（36 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::figs4()`，只读既有产物，不做数值改动；关键列 `route;time;treatment;mean_ISG17_log2_tpm1;sem`。
- **与 v13 核对**：36 treatment-by-time-by-route cells with five replicates each; the module is the mean log2(TPM+1) of the seventeen-gene set.
- **说明**：36 treatment-by-time-by-route cells, 5 replicates each; seventeen-gene module = mean log2(TPM+1)

### Supplementary Fig. S4b — `FigS4b_source.tsv`（30 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::figs4()`，只读既有产物，不做数值改动；关键列 `route;time;treatment;mean_fold_ISG17_log2`。
- **与 v13 核对**：30 treated cells; induction is the treated value minus its own paired untreated control, in log2 units.
- **说明**：30 treated cells; induction = treated minus its own paired untreated control, log2 units

### Supplementary Fig. S4c — `FigS4c_source.tsv`（10 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::figs4()`，只读既有产物，不做数值改动；关键列 `route;isoform;mean_diff;sd_diff;sign_flip_p`。
- **与 v13 核对**：Type I minus type III pooled over three time points (n = 15 pairs) per route; the 95% interval is mean +/- 2.145*se, and asterisks follow the exact sign-flip p, as the caption states.
- **说明**：type I minus type III pooled over three time points (n = 15 pairs); 95% CI = mean +/- 2.145*se; IFNL123_mean pools the three type III isoforms

### Supplementary Fig. S4d — `FigS4d_source.tsv`（72 行）

- **来源**：`<local-path>
- **方法**：`build_source_tables.py::figs4()`，只读既有产物，不做数值改动；关键列 `module;route;time;treatment;mean_module_log2_tpm1`。
- **与 v13 核对**：72 rows: the six-gene housekeeping module and the four-gene receptor module over the same 36 cells, both flat, as the caption states.
- **说明**：negative controls over the same 36 cells: six-gene housekeeping module (HK) and four-gene receptor module (RECEPTOR)

## 不一致清单

下列条目是“图注数字 ↔ 工作区产物”的差异。按要求不改数，只登记。

**1. Fig1b | caption '55.0 A from the second VP24 copy'**

The v13 caption states that P83 lies 55.0 A from the second VP24 copy. No reproduced definition gives 55.0 A: P83 to the nearest atom of the second VP24 chain (chain B) is 22.90 A, the P83 CA-CA distance between the two copies is 80.17 A, P83 to the chain-B centroid is 58.05 A and the whole-chain A-B centroid distance is 36.20 A. The nearest values found in the workspace are 56.62 A (4U2X A-C chain centroids) and 58.04 A (model chain B-Z centroids), neither of which is 'P83 to the second VP24 copy'. The 11.7 A (model) and 11.4 A (4U2X) distances reproduce exactly. Requested action: the author should either supply the calculation that produces 55.0 A or drop/replace the number.

**2. Fig1b | two workspace products disagree on the drawn contact count**

The shipped panel geometry used by Fig1_configuration.py (figures/main_v1_20260926/fig1b_struct/Fig1b_geometry.json) holds 6 inter-chain polar contacts, all 2.73-3.54 A, none involving H140 - this is what the caption's 'six ... below 3.6 A' describes. The PyMOL render script figures/main_v1_20260926/_fig1b/render_fig1b.py lists 7 bonds over a 4.0 A cutoff, including A140-NH2 ... D475-O. The deposited figure and caption follow the 6-contact table; the 7-bond list should be treated as a superseded render input.

**3. Fig3 | caption quadrant numbering differs from the frozen ledger**

The v13 Fig. 3 caption names extremes in 'quadrant I (high Delta/high importin-alpha) and III (low Delta/low importin-alpha)', while the frozen source uses the identifiers Q1_highD_highK and Q4_lowD_lowK for the same two groups (the second-largest of the four low/low quadrants is not called 'III' anywhere in the products). No value changes; only the label convention needs to agree between caption and source.

**4. FigS2b | p value rounding**

The caption gives p = 0.06 for the mode-B paired difference; the exact enumerated value is 0.0584 over 462 assignments. This is rounding, not a discrepancy, and is recorded only so a reader comparing the two figures does not flag it.

**5. FigS2a | fold-induction span pools two negative-group definitions**

The caption's 'fold induction spans 0.50-0.77' is only true when the two negative-group definitions are pooled (mode A 0.50-0.75, mode B 0.60-0.7667). The same applies to the '0.67-0.80' baseline span. Recommended: state that the span is taken over both pre-specified negative-group definitions.

## 仍无法补齐的项

- 23 个 panel 全部给出了数表，**没有一个 panel 完全缺源数据**。
- 唯一的局部缺口是 Fig. 1b 图注中的 **55.0 Å（P83 到第二个 VP24 拷贝）**：该数字在全部工作区产物中找不到，也无法用任何自然定义复算（见不一致清单第 1 条）。因此 `Fig1b_source.tsv` 中该行只给出可复算的替代量（22.90 Å 最近原子、80.17 Å CA-CA），并置 `matches_caption=no`。补齐需要作者提供该距离的计算口径。
- Fig. S2d 是需求矩阵而非测量值；若要把它写成“数据集事实”，需要为每个 1/0 单元格补一条带凭证的注释。本轮按冻结图脚本与图注给出，并在同目录 `.md` 中标注其性质。

## 交付物

- 23 个 `Fig*_source.tsv` 与 23 个同名 `Fig*_source.md`。
- `figure_source_map.tsv`：23 行，列 `figure_panel｜source_file｜source_origin｜n_rows｜key_columns｜matches_manuscript｜notes`。
- `W1_source_data_report.md`（本文件）、`W1_console.txt`（命令与关键输出）。
