# A fixed core and a variable rim: hydroxyl contacts of α-amylases and α-glucosidases set the tolerance to starch substituents

Manuscript, supplementary material, code and data of the article submitted to the *International Journal of
Biological Macromolecules* (Elsevier). Authors: Nedis I. Acevedo-Viloria (Programa de Doctorado en Ingeniería, Universidad del Magdalena; Departamento de Ingeniería Agroindustrial, Universidad de Sucre), Yurimar Ruiz-Rocha (Departamento de Física, Facultad de Educación y Ciencias, Universidad de Sucre;
Doctorado en Ciencias Físicas, SUE Caribe) and Aldo F. Combariza (corresponding; Grupo de Investigación *in silico*,
Departamento de Biología y Química, Facultad de Educación y Ciencias, Universidad de Sucre, Sincelejo, Colombia; Doctorado en Ciencias Físicas, SUE Caribe).

## What the paper is

From an experimental database of 637 PDB entries of 14 starch-active CAZy families (every structure of a family
superposed on one reference, every density-supported sugar unit assigned a subsite −3 … +3 from the catalytic pair,
every residue within 4.0 Å named by its position in the family reference frame and classed by kind of contact), the
paper asks which glucose oxygens each subsite holds.

1. **Residue test (Sect. 3.1, Tables 1–2).** In GH13 (79 entries, 29 groups at ≥ 50 % identity) and GH31 (18 entries,
   6 groups), 11 of the 14 hydrogen-bonding positions at subsites −1/+1 carry one residue type in every group;
   aliphatic positions and positions of every kind at the outer subsites vary (H − Al 0.33, two-sided permutation
   p = 0.004; 0.018 within location). Variants that weaken the result are reported (Table 2, Table S8).
2. **Subsite × oxygen map (Sect. 3.2, Figs 3–4).** Over 114 entries of GH13, GH31, GH14 and GH15, each subsite ×
   oxygen cell gets a contact state (INV fixed, VAR variable, MC main chain, EXP free exposed, BUR free buried).
3. **Map against experiment (Table 3, Tables S6, S9).** 12 of 14 INV cells agree with published products and kinetics
   of substituted substrates; VAR and EXP cells are compatible (19 of 21 overall); a trivial map scores 17 of 21 and
   maps with permuted states 0.59 of their cells; the evidence rests mainly on phosphorylated glucan.
4. **Free volume (Sect. 3.3, Table 4, Table S4)**, a static steric screen of eight substituents; clean agreement with
   the literature in 16 of 36 pairs.
5. **Hypotheses for processing substituted starch (Sect. 3.4).**

The Supplementary material also holds ring distortion, conserved waters, PROPKA, MM/GBSA (S1–S4), the molecular
dynamics consistency check (S6; unfinished runs), the sequence groups with UniProt/EC/GH13 subfamily (Table S10) and the
numbering concordance with TAKA-amylase and porcine pancreatic α-amylase (Table S11).

## Layout

| Folder / file | Content |
|---|---|
| `manuscript/` | `ijbm_main.tex/.pdf` (main text, elsarticle), `ijbm_SI.tex/.pdf` (supplementary material), `ijbm.bib`, `.bbl`; `figures/` holds the figure files the `.tex` include (Figs 1–4, Figs S1–S2) |
| `submission/` | files uploaded besides the manuscript: highlights (`ijbm_highlights.*`), `cover_letter.tex`, `declarations.md` (CRediT, competing interests, generative AI, funding, data availability), graphical abstract (`graphical_abstract.*`) and `suggested_reviewers.md` |
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
                   a43_mmgbsa_subsite.py a44a_build.py a44b_tests.py a44c_tables.py a44d_revision.py
                   a45_hydroxyl_map.py a46_substituent_fit.py a47_md_hydroxyl_contacts.py
                   rcsb_q.py rcsb_q2.py   grip_independence_check/step1-5.py
        si/        si_tables_s2_s3.py
        figures/   fig1_protein_*.py (Fig. 1, Fig. S1)  fig2_v3_*.py (Fig. 2)  fig3_map.py (Fig. 3)
                   figA_*.py (Fig. 4)  figB_fixed_variable.py (Fig. S2)  figdesc_common.py  graphical_abstract.py
                   fig2_v2_*.py fig2_style_*.py fig2_zoom.py fig_contact_classes.py (earlier figure candidates)
        obsolete_build/   scripts that converted an earlier draft into manuscript/ (do not run: they overwrite it)
data/   db/enzimas_experimental.sqlite (41 MB, the whole database) + schema.sql + summaries
        pdb_entries.tsv (the 637 entries: family, method, resolution, release date)
        census/cif_entries.txt (the 645 mmCIF files read, incl. the 8 entries later removed) , bound_linkages*.tsv
        calc/03_frames  calc/01_propka  calc/09_mmgbsa/series  calc/08_md (MD selection and build inputs)
        analysis/33_subsites 40_contact_types 41_control_noncatalytic 42_ring_distortion 43_mmgbsa_subsite
                 44_rule_test (+ revision/) 45_hydroxyl_map 46_substituent_fit 47_md_hydroxyl grip_independence_check
                 + the check and literature reports (*.md; substituted_starch_hydrolysis_literature_v2.md is Table S6)
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
| 4 | Kind of contact (Sect. 2.3; GH14/GH15 in Section S9) | `python code/analysis/a40_contact_types.py`; `python code/analysis/a40_contact_types.py GH14 GH15` | → `data/analysis/40_contact_types/` (`all_*` for the second) |
| 5 | Protein identity groups (Sect. 2.2) | `code/analysis/grip_independence_check/step1.py … step5.py <workdir>` | see note 2 → `data/analysis/grip_independence_check/` |
| 6 | Tables 1 and 2 (residue test) | `python code/analysis/a44a_build.py`; `python code/analysis/a44b_tests.py`; `python code/analysis/a44c_tables.py` | → `data/analysis/44_rule_test/` (`C_identity_T20.tsv`, `D_tests.tsv`, …) |
| 7 | Revision analyses: inner/outer test, location strata, H−Al variants, m vs m_map, map baselines (Tables 2–3, S7–S9), sequence groups with GH13 subfamily (Table S10), numbering concordance (Table S11) | `python code/analysis/a44d_revision.py` (part D reads UniProt once, cached in `revision/uniprot_ncbi_cache.json`, and the CAZy data file `cazy_data.zip`, not shipped) | → `data/analysis/44_rule_test/revision/` |
| 8 | Control at non-catalytic sites (Section S8) | `python code/analysis/a41_control_noncatalytic.py` | → `data/analysis/41_control_noncatalytic/` |
| 9 | Subsite × oxygen map (Sect. 3.2, Figs 3–4) | `python code/analysis/a45_hydroxyl_map.py` | → `data/analysis/45_hydroxyl_map/` |
| 10 | Free volume (Sect. 3.3, Table 4, Table S4) | `python code/analysis/a46_substituent_fit.py` | → `data/analysis/46_substituent_fit/` |
| 11 | MD consistency check (Section S6, Table S5) | systems built from `data/calc/08_md/systems/*/build.leap` (AmberTools 24, GROMACS 2024.4; protocol in Section S12.3); contacts: `python code/analysis/a47_md_hydroxyl_contacts.py` | pocket trajectories `data/analysis/47_md_hydroxyl/traj/` → `data/analysis/47_md_hydroxyl/` |
| 12 | SI S1–S4 (ring distortion, waters, PROPKA, MM/GBSA) | `a42_ring_distortion_subsite.py`; `c01_propka.py`, `c01_collect.py`, `si/si_tables_s2_s3.py`; `c09_series.py`, `c09_series_summary.py`, `a43_mmgbsa_subsite.py` | → `data/analysis/42_*`, `data/calc/01_propka/`, `data/calc/09_mmgbsa/`, `data/analysis/43_*` |
| 12b | Figures | `python code/figures/fig1_protein_layouts.py` (Fig. 1, S1); `fig2_v3_figure.py` (Fig. 2); `fig3_map.py` (Fig. 3); `figA_map_on_molecule.py` (Fig. 4); `figB_fixed_variable.py` (Fig. S2); `graphical_abstract.py` (ChimeraX on the PATH; renders cached, `--rerender` to redo) | → `figures/candidates/`, `manuscript/figures/` |
| 13 | Manuscript | `cd manuscript && pdflatex ijbm_main && bibtex ijbm_main && pdflatex ijbm_main && pdflatex ijbm_main` (same for `ijbm_SI`) | |

Steps 6 (a44c), 7 (parts A–C, E) and 12 (`si_tables_s2_s3.py`, `c09_series_summary.py`, a43) run from the shipped
`data/` alone. Steps 2–6, 8–10 and the structure figures need the mmCIF files of step 0.

Notes.
1. Literature and survey statements rest on `data/analysis/substituted_starch_hydrolysis_literature_v2.md` (Table S6; v1 kept),
   `carbpol_computational_survey.md` (RCSB queries: `code/analysis/rcsb_q.py`, `rcsb_q2.py`),
   `hydroxyl_ts_literature_check.md` and `mutagenesis_literature_check.md`; the other `*.md` files in
   `data/analysis/` are the checks behind individual statements (curation of GH63, sign convention of subsites,
   resolution dependence of the water sites, PROPKA order, clan comparison, protein independence).
2. Steps 1 and 2 of `grip_independence_check` read UniProt entries (`https://rest.uniprot.org/uniprotkb/<ACC>.json`)
   from `<workdir>/dl/`, which are not shipped; their result, `ident.json` and `entries_minus1.tsv`, is in
   `data/analysis/grip_independence_check/` and is what `a44a_build.py` reads.
3. Per-entry PROPKA runs (601 MB), MM/GBSA topologies and coordinates, and the solvated MD topologies, coordinates and
   full trajectories (about 500 MB) are not shipped; they are regenerated from the shipped inputs.

## Data availability (text for the article)

The code and the data needed to reproduce every number, table and figure of this article and its Supplementary
material are available at https://github.com/In-Silico-RG/Starch_Enzyme_Subsites_IJBM: the experimental database
(SQLite), the per-position, per-unit and per-test result tables, and the scripts that produce them. Atomic coordinates
are those deposited in the Protein Data Bank under the accession codes listed in `data/pdb_entries.tsv`; they are
retrieved from the RCSB PDB with `code/download_mmcif.sh`. [PENDIENTE: archived release and DOI]

## Licence

[PENDIENTE]
