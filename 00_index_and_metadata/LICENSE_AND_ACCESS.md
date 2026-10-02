# Licence and access terms

## 1. Author-created content — CC BY 4.0

All content created by the author — derived data tables, analysis code, metadata files and
documentation — is licensed under the Creative Commons Attribution 4.0 International Licence
(CC BY 4.0). You may share and adapt it for any purpose provided you give appropriate credit,
link to the licence and indicate whether changes were made.

Full licence text: https://creativecommons.org/licenses/by/4.0/

## 2. Third-party data are not redistributed

This repository contains **derived results only** for third-party datasets. Raw third-party files
are not redistributed; they remain subject to their original terms and must be obtained from the
original source. Every such file is registered in `LINK_ONLY_RAW.tsv` with its path, size and the
dataset it belongs to (509 files).

| source | what was used | access route | terms |
|---|---|---|---|
| NCBI Protein / Nucleotide | BDBV VP24 protein records; 2026 outbreak complete genomes; reference sequences | https://www.ncbi.nlm.nih.gov/ | US Government public domain (NCBI) |
| UniProt | VP24 references Q05322, B8XCN4, R4QRB9 | https://www.uniprot.org/ | CC BY 4.0 |
| BV-BRC | VP24 records for the six orthoebolavirus species | https://www.bv-brc.org/ | see BV-BRC terms of use |
| ENA | VP24 records for the six orthoebolavirus species | https://www.ebi.ac.uk/ena/ | EMBL-EBI terms of use |
| RCSB PDB | entry 4U2X (VP24–KPNA5) — **redistributed here** | https://www.rcsb.org/ | public domain (CC0) |
| Human Protein Atlas v25.1 | single-cell-type RNA reference (154 cell types), DVP cell-type table | https://www.proteinatlas.org/ | **pending verification — see §3** |
| Cell Ontology | 18,914-term ontology used for label mapping | http://obofoundry.org/ontology/cl.html | CC BY 4.0 |
| CELLxGENE Discover | platform screened for 1,629 public human datasets | https://cellxgene.cziscience.com/ | per-dataset terms |
| NCBI GEO | GSE21158, GSE46599, GSE114905, GSE327707, GSE306664, GSE342661, GSE309699, GSE100839, GSE300073, GSE298600, GSE107549 | https://www.ncbi.nlm.nih.gov/geo/ | per-series terms; raw files not redistributed |
| Normandin et al. 2023 | rhesus macaque multi-organ series | https://doi.org/10.1016/j.xgen.2023.100440 | per-publisher terms |
| WHO Disease Outbreak News | epidemiological figures quoted in the Introduction | https://www.who.int/ | WHO terms of use; cited, not redistributed |
| Research Square preprint rs-10801608 | concurrent report on VP24 variation | https://doi.org/10.21203/rs.3.rs-10801608/v1 | not peer reviewed; cited only |

## 3. Open items requiring author confirmation

1. **Human Protein Atlas v25.1.** Confirm the licence that applies to the specific tables used
   (single-cell-type RNA reference and the Deep Visual Proteomics cell-type table) and whether
   redistribution of the derived tables included here is permitted. Until confirmed, the tables in
   `04_host_reference_and_atlas/hpa_v25_1/` should be treated as *access pending*.
2. **BV-BRC and ENA screens.** Confirm that the per-record tables in
   `02_primary_sequence_data/genus_vp24_screen/` may be redistributed as derived data.
3. **Normandin et al. multi-organ series.** Confirm the publisher's policy on redistributing
   per-sample derived scores.

## 4. Personal and sensitive information

This repository contains no personal data, no human participant data, no patient-identifiable
information and no credentials. Analysis was performed on publicly available data only.

Local absolute paths from the author's machine have been masked to `<local-path>` in all text-like
delivered copies. The analysis scripts in `07_analysis_code_and_environment/analysis/` are the one
exception: they retain their original `ROOT` constant so that the code remains readable as
deposited. See `REPRODUCE.md` §3.4.
