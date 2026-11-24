# Fixed and variable sugar contacts at the subsites of starch-degrading enzymes: a structural basis for processing functionalised starch

Manuscript, supplementary material, code and data of the article submitted to the *International Journal of
Biological Macromolecules* (Elsevier). Authors: Nedis I. Acevedo-Viloria (Programa de Doctorado en Ingeniería, Universidad del Magdalena; Departamento de Ingeniería Agroindustrial, Universidad de Sucre), Yurimar Ruiz-Rocha (Departamento de Física, Facultad de Educación y Ciencias, Universidad de Sucre;
Doctorado en Ciencias Físicas, SUE Caribe) and Aldo F. Combariza (corresponding; Grupo de Investigación *in silico*,
Departamento de Biología y Química, Facultad de Educación y Ciencias, Universidad de Sucre, Sincelejo, Colombia; Doctorado en Ciencias Físicas, SUE Caribe).

## What the paper is

An experimental database of 637 PDB entries of 14 starch-active CAZy families was built; every structure of a family
was superposed on one reference, every bound sugar unit was assigned a subsite (−3 … +3) from the catalytic pair, and
every protein residue within 4.0 Å of a unit was named by its position in the family reference frame and classed by
the kind of contact it makes (side-chain hydrogen bond H, aromatic apolar Ar, aliphatic apolar Al). Across
independent lineages (one vote per group of proteins at ≥ 50 % sequence identity), the side chains that hydrogen-bond
the rings at subsites −1 and +1 are fixed, aromatic residues keep an aromatic ring, and aliphatic residues and
residues of every kind at the outer subsites vary (Tables 1 and 2). The Supplementary material gives ring distortion
by subsite (S1), conserved water sites at the catalytic residues (S2), PROPKA pKa of the catalytic residues (S3) and
an MM/GBSA decomposition of four inhibitor complexes of human lysosomal α-glucosidase (S4).

## Layout

| Folder / file | Content |
|---|---|
| `manuscript/` | `ijbm_main.tex/.pdf` (main text, elsarticle), `ijbm_SI.tex/.pdf` (supplementary material), `ijbm.bib`, `.bbl`; `figures/` holds the figure files the `.tex` include (none at present: the figures are under review in `figures/candidates/`) |
| `submission/` | files uploaded besides the manuscript: highlights (`ijbm_highlights.*`), `cover_letter.tex`, `declarations.md` (CRediT, competing interests, generative AI, funding, data availability), graphical-abstract and suggested-reviewer placeholders |
| `figures/` | figure candidates (`candidates/`), earlier figure versions not included in the current text, and `references/REFERENCES.md` (list of third-party figures consulted for style; the images themselves are not distributed) |
| `code/` | every script the numbers, tables and figures depend on (layout below) |
| `data/` | the database and the result tables those scripts read and write (layout below) |
| `environment.txt` | Python and package versions, PyMOL, ChimeraX, AmberTools |
| `sync_from_project.sh` | maintenance: refreshes `code/` and `data/` from the author's working folder (every source → destination pair listed; applies the path adaptation described below) |

`code/` and `data/` keep the layout of the working project, with `scripts/` → `code/` and `resultados/` → `data/`:

```
code/   part1_table.py part2_structures.py part2_pdb_details.py ligands_inventory.py part3_measure.py
        census_bound_linkages.py build_database.py add_sifts.py add_chembl.py   (database builders)
        download_mmcif.sh                                                       (fetch the deposited mmCIF files)
        calc/      common.py c03_frames.py c01_propka.py c01_collect.py c09_series.py c09_series_summary.py
        analysis/  a33_subsites.py a40_contact_types.py a41_control_noncatalytic.py a42_ring_distortion_subsite.py
                   a43_mmgbsa_subsite.py a44a_build.py a44b_tests.py a44c_tables.py a45_hydroxyl_map.py
                   rcsb_q.py rcsb_q2.py   grip_independence_check/step1-5.py
        si/        si_tables_s2_s3.py
        figures/   fig1_protein_*.py fig2_v2_*.py fig2_style_*.py fig2_zoom.py fig_contact_classes.py
        obsolete_build/   scripts that converted an earlier draft into manuscript/ (do not run: they overwrite it)
data/   db/enzimas_experimental.sqlite (41 MB, the whole database) + schema.sql + summaries
        pdb_entries.tsv (the 637 entries: family, method, resolution, release date)
        census/cif_entries.txt (the 645 mmCIF files read, incl. the 8 entries later removed) , bound_linkages*.tsv
        calc/03_frames  calc/01_propka  calc/09_mmgbsa/series
        analysis/33_subsites 40_contact_types 41_control_noncatalytic 42_ring_distortion 43_mmgbsa_subsite
                 44_rule_test 45_hydroxyl_map grip_independence_check  + the check and literature reports (*.md)
        figures/fig2_v2_work fig2_style_work fig2_work (figure intermediates)
        part2_*.tsv part1_*.tsv cazy_*.tsv ligands/ part3/ (inventory inputs of the builders)
```

Path adaptation: the copies in `code/` differ from the working scripts only in their paths (`resultados/` → `data/`,
`scripts/` → `code/`, figure outputs → `figures/`), and local or cluster home paths and account names are replaced by
placeholders (`<PROJECT_ROOT>`, `~/`, `<cluster-account>`) in `code/` and `data/`. `code/si/si_tables_s2_s3.py` was
written for this package; it reproduces numbers that were computed interactively.

## Reproducing the results

Run from the repository root with the Python environment of `environment.txt`. Every script reads and writes under
`data/` and overwrites the shipped table with the same content.

| # | Item | Command | Reads → writes |
|---|---|---|---|
| 0 | Deposited coordinates (needed from step 2 on) | `bash code/download_mmcif.sh` (`--validation` also fetches the wwPDB validation reports) | `data/census/cif_entries.txt` → `data/census/cif/` |
| 1 | Database (optional: shipped) | `python code/part1_table.py`, `part2_structures.py`, `part2_pdb_details.py`, `ligands_inventory.py`, `part3_measure.py`, `census_bound_linkages.py`, `build_database.py`, `add_sifts.py`, `add_chembl.py`, in this order | CAZy, RCSB, SIFTS, ChEMBL, M-CSA (network) → `data/db/enzimas_experimental.sqlite` |
| 2 | Family frames, catalytic residues, conserved waters | `python code/calc/c03_frames.py` | db, mmCIF → `data/calc/03_frames/` |
| 3 | Subsites and mapped contacts | `python code/analysis/a33_subsites.py` | db, frames → `data/analysis/33_subsites/` |
| 4 | Kind of contact (Methods; GH14/GH15 in Sect. 3.4) | `python code/analysis/a40_contact_types.py`; `python code/analysis/a40_contact_types.py GH14 GH15` | → `data/analysis/40_contact_types/` (`all_*` for the second) |
| 5 | Protein identity groups (Sect. 3.2) | `code/analysis/grip_independence_check/step1.py … step5.py <workdir>` | see note 2 → `data/analysis/grip_independence_check/` |
| 6 | Table 1 and Table 2 (rule test) | `python code/analysis/a44a_build.py`; `python code/analysis/a44b_tests.py`; `python code/analysis/a44c_tables.py` | → `data/analysis/44_rule_test/` (`C_identity_T20.tsv`, `D_tests.tsv`, `tab_rule.tex`, `tab_tests.tex`, `REPORT.md`) |
| 7 | Control at non-catalytic sites (Sect. 3.3) | `python code/analysis/a41_control_noncatalytic.py` | → `data/analysis/41_control_noncatalytic/` |
| 8 | Hydroxyl map, subsite × sugar oxygen (supporting) | `python code/analysis/a45_hydroxyl_map.py` | → `data/analysis/45_hydroxyl_map/` |
| 9 | SI S1, ring distortion | `python code/analysis/a42_ring_distortion_subsite.py` | db, `census/bound_linkages_checked.tsv` → `data/analysis/42_ring_distortion/summary.tsv` |
| 10 | SI S2 (water sites, Table S2) and S3 (PROPKA, Table S3) | `python code/calc/c01_propka.py <PDB …>`; `python code/calc/c01_collect.py`; then `python code/si/si_tables_s2_s3.py` | → `data/calc/01_propka/propka_all.tsv`; printed tables |
| 11 | SI S4, MM/GBSA (Table S4) | `AMBER_ENV=<AmberTools prefix> python code/calc/c09_series.py 5NN5 NOJ` (and 5NN6 MIG, 9GTC A1IO7, 9GTW A1IPC); `python code/calc/c09_series_summary.py`; `python code/analysis/a43_mmgbsa_subsite.py` | → `data/calc/09_mmgbsa/series/`, `data/analysis/43_mmgbsa_subsite/` |
| 12 | Figure candidates | `python code/figures/fig2_v2_data.py`; `fig2_v2_layouts.py`; `fig2_zoom.py`; `fig2_style_[1-4]_*.py`; `fig1_protein_struct.py`; `fig1_protein_frame.py`; `fig1_protein_layouts.py`; `fig_contact_classes.py` (PyMOL and ChimeraX on the PATH) | → `figures/candidates/`, `data/figures/` |
| 13 | Manuscript | `cd manuscript && pdflatex ijbm_main && bibtex ijbm_main && pdflatex ijbm_main && pdflatex ijbm_main` (same for `ijbm_SI`) | |

Steps 6 (a44c), 9, 10 (`si_tables_s2_s3.py`) and 11 (`c09_series_summary.py`, a43) run from the shipped `data/` alone
and were checked to reproduce the shipped tables and the values printed in the manuscript. Steps 2–8 and 12 need the
mmCIF files of step 0.

Notes.
1. Literature and survey statements rest on `data/analysis/substituted_starch_hydrolysis_literature.md`,
   `carbpol_computational_survey.md` (RCSB queries: `code/analysis/rcsb_q.py`, `rcsb_q2.py`),
   `hydroxyl_ts_literature_check.md` and `mutagenesis_literature_check.md`; the other `*.md` files in
   `data/analysis/` are the checks behind individual statements (curation of GH63, sign convention of subsites,
   resolution dependence of the water sites, PROPKA order, clan comparison, protein independence).
2. Steps 1 and 2 of `grip_independence_check` read UniProt entries (`https://rest.uniprot.org/uniprotkb/<ACC>.json`)
   from `<workdir>/dl/`, which are not shipped; their result, `ident.json` and `entries_minus1.tsv`, is in
   `data/analysis/grip_independence_check/` and is what `a44a_build.py` reads.
3. Per-entry PROPKA runs (601 MB) and MM/GBSA topologies, coordinates and trajectories are not shipped; they are
   regenerated by steps 10 and 11.

## Data availability (text for the article)

The code and the data needed to reproduce every number, table and figure of this article and its Supplementary
material are available at https://github.com/In-Silico-RG/Starch_Enzyme_Subsites_IJBM: the experimental database
(SQLite), the per-position, per-unit and per-test result tables, and the scripts that produce them. Atomic coordinates
are those deposited in the Protein Data Bank under the accession codes listed in `data/pdb_entries.tsv`; they are
retrieved from the RCSB PDB with `code/download_mmcif.sh`. [PENDIENTE: archived release and DOI]

## Licence

[PENDIENTE]
