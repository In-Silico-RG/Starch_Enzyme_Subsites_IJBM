# Superseded files (moved 2026-11-24, IJBM round-2 revision)

These four files were written by `a44b_tests.py` (`E_inner_vs_outer_H.tsv`) and `a44c_tables.py` (`REPORT.md`,
`rule_section.tex`, `tab_tests.tex`) before the editorial revisions of 2026-11-24. They carry numbers and wording that
the article no longer uses: 13 vs 4 hydrogen-bonding positions at -1/+1 vs outer (the subsites re-derived from contacts
at -1/+1 only, which put GH13 Lys182 among the outer positions), p = 0.0015 for that test, and the Ar(arom)-Al
comparison, whose permutation null is not valid (Section S7 of the Supplementary material).

They are kept for the record only. The numbers of the article come from:
- `C_identity_T20.tsv`, `D_tests.tsv` (positions, m, earlier test runs) — a44b_tests.py;
- `revision/` — a44d_revision.py (14 vs 3 inner/outer test, baselines, groups, numbering);
- `revision2/` — a44e_round2.py (every reported permutation test with 100 000 permutations and one seed, model of m,
  McNemar, per-family and abstract-only scores, alpha-amylase-only map, MD numbers).

Re-running a44b_tests.py or a44c_tables.py writes these files again at the folder top level; those copies are
superseded in the same way.
