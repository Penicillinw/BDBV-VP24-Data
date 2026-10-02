# Corrections applied while building this deposit

This deposit is aligned to **v15** of the manuscript. Two rounds of independent checks ran
against the deposited data and found five places where a figure caption and the underlying
tables did not agree, plus one place where a figure cell was contradicted by its own voucher.
Everything below was recomputed from source before being changed; nothing was adjusted by eye.

## 1. Fig. 1b — removal of an unreproducible distance

**Was**: "…P83 lies 11.7 Å from the nearest KPNA5 atom in the model, **55.0 Å from the second VP24
copy** and 11.4 Å in the experimental coordinates…"

**Finding**: an exhaustive sweep of 66 distance definitions — three copy pairings in PDB 4U2X,
the shipped homology model, both template-transfer models, and metrics covering minimum heavy-atom
distance, Cα–Cα, Cα–centroid, residue-centroid, whole-chain centroid and backbone-atom pairings —
produced no value within 55.0 ± 0.1 Å. The literal "P83 to the second copy" distances are 22.9 Å
(nearest heavy atom) and 80.2 Å (Cα–Cα). The single value in the whole workspace that lands on
55.0 Å is the distance from the **second** copy's own P83 to the nearest KPNA5 atom (54.98 Å),
which is a different quantity.

**Now**: the clause is deleted; the sentence keeps the two reproducible distances (11.7 Å in the
model, 11.4 Å in the experimental coordinates) and the eight-ångström exclusion. Evidence:
`03_primary_structure_and_interface/bdbv_vp24_model/…x1-p83-distance-sweep…tsv`.

## 2. Fig. 1b — the contact count and its criterion

**Was**: "dashed yellow lines mark six inter-chain polar contacts below 3.6 Å" (the count was
right, but the criterion was left implicit and a superseded render input listed seven contacts).

**Finding**: recomputed from scratch on the structure the panel actually draws, the four variant
positions make exactly **six** inter-chain N/O···N/O contacts below 3.6 Å — two from Q135 and four
from R184; H140 and A141 make none. The seven-bond alternative came from a different structure
(experimental EBOV side chains, where position 140 is arginine rather than histidine) under a wider
4.0 Å cutoff, and was not even self-consistent for that structure.

**Now**: the caption states the criterion and the attribution explicitly. Evidence:
`03_primary_structure_and_interface/pdb_4u2x_experimental/…x2-contacts-authoritative…tsv`.

## 3. Fig. 3a — quadrant labels

**Was**: "the three most extreme members of **quadrant I** (high Δ/high importin-α) and **III**
(low Δ/low importin-α) are named", while the frozen source table labels the same groups
`Q1_highD_highK` and `Q4_lowD_lowK`.

**Now**: the caption names the two quadrants by their definition rather than by Roman numerals, so
caption and source agree. No value changed.

## 4. Fig. S2d — a figure cell corrected against its voucher

**Was**: the decisive-test requirement matrix scored 342661 at **3/7** requirements met, with
"Publicly released" recorded as not satisfied.

**Finding**: the voucher contradicts the frozen cell. GSE342661's processed TPM matrix
(`GSE342661_counts_tpm.matrix.gz`, 10,325,209 B) was released on GEO on **25 September 2026** and
is present in the analysis workspace; the project's own note of that date records that the
"not obtainable" statement had become factually wrong. All 35 matrix cells were re-vouchered; 34
were confirmed and this one was refuted.

**Now**: the cell is 1, the dataset meets **4/7**, and Fig. S2 panel d has been redrawn to match.
The caption also records the exact `p = 0.058` (previously rounded to 0.06) and states that the
0.50–0.77 fold-induction span is taken over both pre-specified negative-group definitions.
Evidence: `01_source_data_by_figure/FigS2/…figs2d…tsv`, `06_analysis_data/coverage_audit/…x3-figs2d-provenance…tsv`.

## 5. Fig. 1a — denominator chain (no change)

The caption's denominator chain (3,522 retrieved → 32 fragments excluded → 3,490 assessed → 19 with
an uncalled position → 3,471 called at all four positions → 3,440 non-Bundibugyo) was recomputed
species by species and matches digit for digit. Nothing was changed.
Evidence: `02_primary_sequence_data/call_tables/…x5-fig1a-denominator-recompute…tsv`.

## 6. Material that was previously excluded and is now deposited

The inclusion rule for this repository is: *delete the file — can the reported numbers still be
recomputed?* Applying that rule to the exclusion ledger brought four groups back in, and one
group's recomputation had to be honestly bounded:

| previously excluded | now | why |
|---|---|---|
| pilot-MD minimised coordinates (**49 PDB files, 12.1 MB**; nine of them are empty chain-D placeholders of a few bytes) | deposited under `03_primary_structure_and_interface/md_pilot/` | the interface-contact integers reported in `openmm_results.json` are reproducible from these coordinates digit for digit; without them those eight integers cannot be recomputed |
| `protein_layer_scan` derived tables | deposited under `04_host_reference_and_atlas/hpa_v25_1/`; its 84 MB of raw downloads stay `link_only` | these tables back the mRNA–protein layer comparison reported in Methods |
| Fig. 1a denominator provenance | deposited under `02_primary_sequence_data/` | the caption's denominator chain previously had no deposited table |
| round-X recomputations and the adjudication document | deposited under `06_analysis_data/`, `03_…/`, `02_…/` and `01_source_data_by_figure/_shared/` | the evidence for corrections 1–5 |

**Honest limit**: 41 quantities reported by the stability pilot (`pilot_run.json`,
`pilot_bench.json`: local RMSD, native-contact fraction, minimised energy, void volume) **cannot be
recomputed** from this deposit, because **this pilot's own 1 ns trajectories were not retained in
the workspace** — only its minimised coordinates are. (Other, unrelated GROMACS campaigns in the
project do have trajectories; they are not the runs behind these numbers.) They are marked
`not_recomputable` in `03_primary_structure_and_interface/md_pilot/…x4-md-recomputed-metrics…tsv`
rather than being approximated. The manuscript reports the pilot qualitatively and draws no
number from it, so no reported value depends on them. Any future claim that needs these numbers
must first deposit the trajectories.
