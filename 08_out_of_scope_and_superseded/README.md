# Out-of-scope and superseded material

Nothing in this folder is redistributed. It exists so that a reader can see **what was deliberately
left out and why**.

## Files

| file | rows | meaning |
|---|---|---|
| `OUT_OF_SCOPE_FILES.tsv` | 1,792 | every candidate file excluded from the package, with its size, the decision (`no` or `pending`) and the reason |
| `OUT_OF_SCOPE_summary.tsv` | 31 | the same ledger aggregated by decision, category and reason |

## Categories of exclusion

1. **Rendering and manuscript artefacts** — figures (`.png/.pdf/.tiff/.svg/.eps`), figure-plotting
   scripts and style files, captions, manuscript and supplementary DOCX/PDF/text, DOCX assembly and
   export scripts, and rendering QA output. These are outputs, not data: deleting them changes no
   reported number.
2. **Third-party raw data** — GEO series matrices and single-cell `.h5` files, HPA downloads and
   NCBI/UniProt/ENA/BV-BRC query responses. Registered in
   `../00_index_and_metadata/LINK_ONLY_RAW.tsv` with accession and checksum, but not redistributed.
3. **Unrelated analyses in the project tree** — structural-prediction campaigns (AlphaFold,
   AlphaFold3, Boltz), molecular docking (HADDOCK, ClusPro, LightDock and others), VP30/VP35/VP40/GP
   sequence and glycosylation work, FoldX/APBS/NMA/electrostatics calculations, and binding-free-energy
   or interface-scoring output (MM/PBSA, RBFE, PRODIGY, shape complementarity, decoy scoring). These
   belong to analyses that are not reported in the accompanying manuscript; the manuscript's
   Methods state that no cross-structure energy comparison was performed.
4. **Superseded and backup copies** — earlier manuscript versions, `_backup_*`, `_superseded_*`,
   `_stage_*` and render scratch directories.
5. **Process and QA files** — alignment/collision audit JSON, verification logs, round bookkeeping
   and agent working notes. Regenerable and not data.
6. **Pending** — pilot molecular-dynamics coordinate sets. The manuscript reports the pilot as
   uninformative and states that it supports no claim; only the summary JSON is included
   (`03_primary_structure_and_interface/`). The coordinates remain available on request.

## Reinstating an entry

If a journal or a reader requires one of these files, the decision can be revisited file by file:
each row carries the original path, so the artefact can be recovered from the author's project tree
without re-running anything.
