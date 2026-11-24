# Ring puckering by subsite (exploratory check, 2026-11-23)

**Problem.** Glycosidase mechanism textbooks predict that the enzyme distorts the sugar ring away from
its relaxed chair specifically at the subsite next to the catalytic nucleophile/base (subsite −1), to help
reach the transition state. Do our own data show this?

**Data.** Census ring puckering θ (Cremer–Pople; θ ≈ 0° is the ⁴C₁ chair) from
`resultados/census/bound_linkages_checked.tsv`, joined to subsite labels from
`resultados/analysis/33_subsites/units.tsv`, matched through the deposited Glc-α(1→4)-Glc covalent bonds in
the database table `bond` (chain/residue-number keys, not string parsing of the census labels).

**Method.** For each of the two rings of a linkage (non-reducing/C1 side, reducing/O4 side), look up which
subsite the corresponding residue occupies in its family's frame, and group θ into "subsite −1" vs "other
negative" vs "other positive" subsites. θ is folded to its distance from the nearest chair,
min(θ, 180° − θ) — a ring near the ⁴C₁ chair (θ≈0°) and one near the ¹C4 chair (θ≈180°) are both undistorted;
only θ near 90° is a boat/skew-boat. **Correction (2026-11-23):** the first pass used raw θ and would have
mislabelled a ¹C4 chair as "distorted"; checked against the whole census, 10 of 2 688 ring values are
θ>140° (true ¹C4 chairs). Refolding changed the numbers below only slightly; the finding stands.

**Result.** Only 208 of 1 344 census linkages could be matched to a subsite at all (many families/entries
have no subsite assignment in `units.tsv`), so this is a partial, exploratory look, not a census-wide result.

| Family | Ring | At subsite −1 | Other subsites |
|---|---|---|---|
| GH13 | reducing (O4 side) | n=33, median 17.7°, 13/33 (39%) above 20° | n=66, median 7.5–9.1°, 7/66 (11%) above 20° |
| GH13 | non-reducing (C1 side) | n=3 (too few to trust) | n=95, median 5.3–7.8°, 4/95 (4%) above 20° |
| GH57 | reducing | n=14, median 11.7°, 0/14 above 20° | n=28, median 5.3–7.6°, 0/28 above 20° |
| GH31 | — | no linkage at subsite −1 matched | n=58 |

**Reading (mine, kept apart from the numbers).** In GH13, the reducing-side ring at subsite −1 is
markedly more distorted than at any other subsite — consistent with the textbook prediction. GH57 shows no
such enrichment at −1. GH31 cannot be checked this way: no linkage with subsite data lands at −1.

**Limits, stated flatly.**
- Coverage is thin (208/1344); most of the census has no subsite label yet.
- The GH13 non-reducing/−1 cell (n=3) is driven by one or two structures and is not evidence by itself —
  not identified which PDB entries yet.
- No check for whether the GH13 −1/reducing pattern is itself driven by one or two entries (e.g. the
  already-known distorted cases) rather than being spread across many proteins.
- Not controlled for resolution or density support in this pass.
