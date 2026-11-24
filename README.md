# A structural map of glucose-hydroxyl recognition in α-amylases and α-glucosidases and its implications for the enzymatic hydrolysis of chemically modified starch

Manuscript, supplementary material, code and data of the article submitted to the *International Journal of
Biological Macromolecules* (Elsevier). Authors: Nedis I. Acevedo-Viloria (Programa de Doctorado en Ingeniería, Universidad del Magdalena; Departamento de Ingeniería Agroindustrial, Universidad de Sucre), Yurimar Ruiz-Rocha (Departamento de Física, Facultad de Educación y Ciencias, Universidad de Sucre;
Doctorado en Ciencias Físicas, SUE Caribe) and Aldo F. Combariza (corresponding; Grupo de Investigación *in silico*,
Departamento de Biología y Química, Facultad de Educación y Ciencias, Universidad de Sucre, Sincelejo, Colombia; Doctorado en Ciencias Físicas, SUE Caribe).

## What the paper is

From an experimental database of 637 PDB entries of 14 starch-active CAZy families (every structure of a family
superposed on one reference, every density-supported sugar unit assigned a subsite −3 … +3 from the catalytic pair,
every residue within 4.0 Å named by its position in the family reference frame and classed by kind of contact), the
paper asks which glucose oxygens each subsite holds.

1. **Residue test (Sect. 3.1, Table 1, Table S1).** In GH13 (79 entries, 29 groups at ≥ 50 % identity) and GH31 (18
   entries, 6 groups), 11 of the 14 hydrogen-bonding positions at subsites −1/+1 carry one residue type in every group
   (mean m 0.94 against 0.57 for the 3 outer-only ones, two-sided permutation p = 0.009); aliphatic positions vary
   (pooled H − Al 0.33, p = 0.005; 100 000 permutations, one seed). Settings that weaken the result and a model of m
   on location, family and class are reported (Table 1, Tables S6–S7).
2. **Subsite × oxygen map (Sect. 3.2, Figs 3–4).** Over 114 entries of GH13, GH31, GH14 and GH15, each subsite ×
   oxygen cell gets a contact state (INV fixed, VAR variable, MC main chain, EXP free exposed, BUR free buried).
3. **Map against experiment (Table 2, Tables S4, S8).** 12 of 14 INV cells agree with published products and kinetics
   of substituted substrates (GH13 4/4, GH14 4/5, GH15 4/5); a trivial map that blocks every oxygen at −1/+1 agrees in
   10 of 12 of its blocked cells and differs from the published map at two GH14 cells only (McNemar p = 0.50); the
   evidence rests mainly on phosphorylated glucan; 17 of 18 cells agree without the three abstract-only cells.
4. **Free volume (Sect. 3.3, Table 3, Table S2)**, a static steric screen of eight substituents; clean agreement with
   the literature in 16 of 36 pairs.
5. **Hypotheses for processing substituted starch (Sect. 3.4).**

The Supplementary material also holds the contact positions (Table S1), unfinished molecular dynamics of three
glucose-only complexes (S3; not used in the main text), the GH13 map over the 13 α-amylase groups (Table S9), the
distinct enzymes with UniProt/EC/GH13 subfamily (Table S10) and the numbering concordance with TAKA-amylase and porcine
pancreatic α-amylase (Table S11). Ring distortion, conserved waters, PROPKA and MM/GBSA, which the article does not
use, were removed from the Supplementary material in the round-2 revision; their scripts and results remain in `code/`
and `data/` (steps 12 and 2).

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
                   a43_mmgbsa_subsite.py a44a_build.py a44b_tests.py a44c_tables.py a44d_revision.py a44e_round2.py
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
                 44_rule_test (+ revision/, revision2/, superseded/) 45_hydroxyl_map 46_substituent_fit 47_md_hydroxyl grip_independence_check
                 + the check and literature reports (*.md; substituted_starch_hydrolysis_literature_v2.md is Table S4)
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
| 4 | Kind of contact (Sect. 2.3; GH14/GH15 in Section S4) | `python code/analysis/a40_contact_types.py`; `python code/analysis/a40_contact_types.py GH14 GH15` | → `data/analysis/40_contact_types/` (`all_*` for the second) |
| 5 | Protein identity groups (Sect. 2.2) | `code/analysis/grip_independence_check/step1.py … step5.py <workdir>` | see note 2 → `data/analysis/grip_independence_check/` |
| 6 | Table S1, earlier test runs (residue test) | `python code/analysis/a44a_build.py`; `python code/analysis/a44b_tests.py`; `python code/analysis/a44c_tables.py` | → `data/analysis/44_rule_test/` (`C_identity_T20.tsv`, `D_tests.tsv`, …; the pre-revision `REPORT.md`, `rule_section.tex`, `tab_tests.tex`, `E_inner_vs_outer_H.tsv` are in `superseded/`, see its `NOTE.md`) |
| 7 | Revision analyses: m vs m_map (Table S5), map baselines (Table 2, Table S8), distinct enzymes with GH13 subfamily (Table S10), numbering concordance (Table S11) | `python code/analysis/a44d_revision.py` (part D reads UniProt once, cached in `revision/uniprot_ncbi_cache.json`, and the CAZy data file `cazy_data.zip`, not shipped) | → `data/analysis/44_rule_test/revision/` |
| 7b | Round-2 analyses: every reported permutation test with 100 000 permutations and one seed (Table 1, Table S6), model of m (Table S7), McNemar, per-family and abstract-only scores (Table 2, Section S8), α-amylase-only map (Table S9), MD numbers (Section S3) | `python code/analysis/a44e_round2.py` (runs from the shipped `data/` alone) | → `data/analysis/44_rule_test/revision2/` |
| 8 | Control at non-catalytic sites (not in the article since round 2) | `python code/analysis/a41_control_noncatalytic.py` | → `data/analysis/41_control_noncatalytic/` |
| 9 | Subsite × oxygen map (Sect. 3.2, Figs 3–4) | `python code/analysis/a45_hydroxyl_map.py` | → `data/analysis/45_hydroxyl_map/` |
| 10 | Free volume (Sect. 3.3, Table 3, Table S2) | `python code/analysis/a46_substituent_fit.py` | → `data/analysis/46_substituent_fit/` |
| 11 | MD of three complexes (Section S3, Table S3; not used in the main text) | systems built from `data/calc/08_md/systems/*/build.leap` (AmberTools 24, GROMACS 2024.4; protocol in Section S5.3); contacts: `python code/analysis/a47_md_hydroxyl_contacts.py` | pocket trajectories `data/analysis/47_md_hydroxyl/traj/` → `data/analysis/47_md_hydroxyl/` |
| 12 | Ring distortion, waters, PROPKA, MM/GBSA (not in the article since round 2) | `a42_ring_distortion_subsite.py`; `c01_propka.py`, `c01_collect.py`, `si/si_tables_s2_s3.py`; `c09_series.py`, `c09_series_summary.py`, `a43_mmgbsa_subsite.py` | → `data/analysis/42_*`, `data/calc/01_propka/`, `data/calc/09_mmgbsa/`, `data/analysis/43_*` |
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
   full trajectories (about 500 MB) are not shipped; they are regenerated from the shipped inputs. The MD numbers of
   Section S3 are read from the shipped pocket trajectories (`data/analysis/47_md_hydroxyl/traj/`) and the shipped
   per-cell table; the MD is not used in the main text.

## Data availability (text for the article)

The code and the data needed to reproduce every number, table and figure of this article and its Supplementary
material are available at https://github.com/In-Silico-RG/Starch_Enzyme_Subsites_IJBM: the experimental database
(SQLite), the per-position, per-unit and per-test result tables, and the scripts that produce them. Atomic coordinates
are those deposited in the Protein Data Bank under the accession codes listed in `data/pdb_entries.tsv`; they are
retrieved from the RCSB PDB with `code/download_mmcif.sh`. [PENDIENTE: archived release and DOI]

## Licence

Code (`code/`, `sync_from_project.sh`) under the MIT License; data, result tables, figures and text (`data/`,
`figures/`, `manuscript/`, `submission/`) under CC BY 4.0. See `LICENSE`. PDB coordinates and the records of CAZy,
UniProt, SIFTS and ChEMBL keep their providers' terms.
