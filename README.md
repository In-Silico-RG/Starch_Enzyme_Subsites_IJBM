# Fixed and variable sugar contacts at the subsites of starch-degrading enzymes: a structural basis for processing functionalised starch

Manuscript for the *International Journal of Biological Macromolecules* (Elsevier). Carbohydrate Polymers is the
stretch option. Author: Aldo F. Combariza (Grupo de Investigación *in silico*, Universidad de Sucre; Doctorado en
Ciencias Físicas, SUE Caribe).

This repository holds the manuscript. It is synced from the folder `Proyecto_F_Enzimas_Amiloliticas/manuscript_ijbm/`
of the project repository In-Silico-RG/Enzimas_Amiloliticas, where the database, scripts and results also live.

## Files
- `ijbm_main.tex` / `.pdf`: main text (elsarticle, numbered references); `ijbm_SI.tex` / `.pdf`: supplementary
  material; `ijbm_highlights.*`: highlights; `ijbm.bib`: references.
- `figures/`: figures in the manuscript; `figures/candidates/`: figure candidates under review.
- `build/`: the scripts that converted the earlier draft. They are obsolete: running them overwrites the
  manuscript.

Build: `pdflatex ijbm_main && bibtex ijbm_main && pdflatex ijbm_main && pdflatex ijbm_main` (the same for `ijbm_SI`).

## Where each result comes from (paths in Enzimas_Amiloliticas)

| Section | Scripts | Results |
|---|---|---|
| Database, superposition, subsites | `scripts/analysis/a3*`, `scripts/calc/` | `resultados/analysis/33_subsites/`, `resultados/calc/03_frames/` |
| Kind of contact | `scripts/analysis/a40_contact_types.py` | `resultados/analysis/40_contact_types/` |
| Conservation by class (Table 1, Table 2) | `scripts/analysis/a44a_build.py`, `a44b_tests.py`, `a44c_tables.py` | `resultados/analysis/44_rule_test/` |
| Control at non-catalytic sites | `scripts/analysis/a41_control_noncatalytic.py` | `resultados/analysis/41_control_noncatalytic/` |
| Hydroxyl map (subsite × oxygen) | `scripts/analysis/a45_hydroxyl_map.py` | `resultados/analysis/45_hydroxyl_map/` |
| Ring distortion (SI S1) | `scripts/analysis/a42_ring_distortion_subsite.py` | `resultados/analysis/42_ring_distortion/` |
| PROPKA (SI S3) | — | `resultados/calc/01_propka/` |
| MM/GBSA (SI S4) | `scripts/calc/c09_series.py`, `scripts/analysis/a43_mmgbsa_subsite.py` | `resultados/calc/09_mmgbsa/series/`, `resultados/analysis/43_mmgbsa_subsite/` |
| Literature on substituted starch | — | `resultados/analysis/substituted_starch_hydrolysis_literature.md` |
| Figures | `scripts/figures/` | `resultados/figures/` |

Reviews of the draft, with their adjudication, are in `docs/reviews/` of Enzimas_Amiloliticas, and the process
log is `docs/LOG.md` of the Amilosa repository.
