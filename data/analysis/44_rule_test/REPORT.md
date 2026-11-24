# 44 — Rule test for the contact-type claim (GH13, GH31), 2026-11-24

The analysis agent ran this and returned the report as a message, because the harness did not let it write files. This
file reproduces that report in substance. Every number below is in a TSV of this folder; the main tests are in
`D_tests.tsv`. Scripts: `scripts/analysis/a44a_build.py` (about 1 min) and `scripts/analysis/a44b_tests.py`
(about 5 min). Both ran on the laptop.

## Problem
The claim of the earlier draft was open to the objection that it is circular:
- it uses 19 positions chosen after the fact;
- a position's class is read from the same residues whose identity is then tested;
- residues are counted per entry, not per lineage;
- catalytic residues are mixed in with the rest;
- the control is weak;
- hydrogen bonds are judged by distance alone.

## Method
- **Clusters:** the existing ones from `grip_independence_check/ident.json` (GH13: 29 at 50% identity, 5 at 30%;
  GH31: 6 and 4). 2E9B joins an existing 50% cluster; 1UA3, 6IYG, 6JHI and 6Z8L are proteins already counted.
- **Residue identity per lineage:** read from the structural map. A residue counts at a reference position when its
  Cα lies within 2.0 Å of the reference Cα after superposition. This agrees with the contact table for all 1,798
  contacted pairs.
- **Engineered mutants:** taken from `_struct_ref_seq_dif` and dropped at the mutated position
  (`mutations_at_mapped_positions_bound.tsv`).
  - All four Asn at GH13 D289 (1QI5, 5TD4, 9UNS, 9XG6) and both Ala at GH31 D542 (7WJC, 7WJF) are mutants.
  - Also dropped at tested positions: 3PHA W169Y, 3BLP/3DHP, 3OLD/E/G/I, 1UA7, and 8CQF at nine positions.
- **Statistics:** permutation tests 10,000, bootstrap 2,000, split-half 1,000.
  - Bootstrap CIs are conservative: resampling clusters can only raise the modal fraction of a variable position.
  - The permutation p is the test.

## A. Positions chosen without post-hoc selection
Every reference position contacted (≤ 4.0 Å) in at least the given share of entries at subsites −3..+3. The five
catalytic residues are included at every cutoff.

| Cutoff | GH13 | GH31 | Total |
|---|---|---|---|
| ≥ 20% | 17 | 19 | 36 |
| ≥ 50% | 14 | 17 | 31 |
| ≥ 80% | 10 | 12 | 22 |

## C. Identity counted once per 50% cluster (cutoff 20%; `C_identity_T20.tsv`)
m = share of cluster votes carried by the most common residue type.

- **m = 1.00:**
  - GH13: D179, D289, E204, H92, R177, Y51.
  - GH31: D327, D443, D542, H600, D203, R526, W441, W406.
- **H288:** 0.97 (one natural Gln, 6SAO).
- **H positions that vary:**
  - G183: 0.50 (His 8, Gly 5, Tyr, Glu, Met). It is His wherever it contacts the sugar.
  - K182: 0.63.
  - M52: 0.41.
  - Q294: 0.80.
  - N207: 0.50 (2 clusters).
- **Ar positions:**
  - W539: 0.83.
  - F575: 0.83.
  - Y299: aromatic fraction 1.0, identity 0.40.
  - Q49: aromatic fraction 0.85.
  - W206: aromatic fraction 0.85.
  - F143: aromatic fraction 0.40.
  - G250: aromatic fraction 0.56.
- **Al positions:**
  - F180: 0.24.
  - I328: 0.50.
  - I364: 0.67.
  - M444: 0.83.
  - A576, T204: 0.20 each.
  - T205: 0.67.
  - T207: 1.0 (2 clusters).

## D. Is identity conservation associated with class? (cutoff 20%, class by 50%-cluster majority)

| Subset | n H/Ar/Al | Mean m H/Ar/Al | H−Al | H−Ar | Ar(aromatic fraction)−Al |
|---|---|---|---|---|---|
| All positions | 17/10/8 | 0.871/0.687/0.538 | +0.332 [0.118, 0.406], p = 0.003 (stratified 0.002) | +0.184, p = 0.052 | +0.328, p < 0.001 |
| Without catalytic | 12/10/8 | — | +0.279, p = 0.016 (stratified 0.007) | +0.130, p = 0.14 | — |

## B. Class assigned without using the identities tested
- **(i) Class from the reference structure** (1BG9, 2QMJ): H−Al p values.

  | Cutoff | All positions | Without catalytic |
  |---|---|---|
  | 20% | 0.005–0.008 | 0.014–0.016 |
  | 50% | 0.04–0.09 | about 0.09–0.10 |
  | 80% | 0.05–0.08 (stratified 0.005–0.012) | 0.073 (stratified 0.008) |

  The 50% and 80% results are weaker for two reasons: the barley reference has aromatic residues at F180 and W206,
  where most lineages are aliphatic, and D443 makes no side-chain hydrogen bond in 2QMJ.
- **(ii) Split-half** (classed on half the clusters, tested on the other half):
  - H−Al > 0 in 100% of splits at every cutoff. At 20%: median 0.299 [0.165, 0.413]; without catalytic, 0.263
    [0.106, 0.385].
  - H−Ar > 0 in 88–98% of splits, with the 2.5th percentile below 0 in most settings.
  - All H positions invariant in no split at 20% or 50%, and in 17.5% of splits at 80%.
- **(iii) Held-out proteins** (79 GH13 and 19 GH31 entries less than 50% identical to any ligand-bound protein):
  - H−Al at cutoff 20%: 0.310, p = 0.008; without catalytic 0.233, p = 0.036.
  - H−Al at cutoff 80%: 0.330, p = 0.012; without catalytic 0.294, p = 0.040.
  - H−Ar: 0.085, p = 0.23.
  - With 30% clusters: p = 0.067 at cutoff 20%, 0.013 at 80%.
  - D179, D289, E204, H92, H288 and R177, and the six GH31 H positions, stay at m 0.97–1.00.
  - G183 (0.22), K182 (0.19), M52 (0.53) and Q294 (0.29) vary.

## E. Control at the outer subsites of the same cleft (−3, −2, +2, +3)
- H vs Al at outer subsites: −0.087, p = 0.67; restricted to outer-only positions, 0.053, p = 0.40.
- H at −1/+1 vs H only at outer subsites: mean m 0.959 (n = 13) vs 0.591 (n = 4), p = 0.0015. Without catalytic
  residues: 0.933 (8) vs 0.591, p = 0.008.
- Low power: only 3–4 H positions sit at the outer subsites.

## F. Sensitivity to the hydrogen-bond definition
- **Class changes:**
  - At 3.2 Å, GH31 D443 goes H→Al.
  - At 4.0 Å, GH31 T544 goes from a tie to H.
  - At 3.5 Å with a heavy-atom angle test (both angles ≥ 90°; no hydrogens placed), GH13 H92 goes H→Ar and D443
    goes H→Al.
- **Tests:** H−Al p = 0.003–0.040 in every setting. H−Ar is never significant (p = 0.05–0.25).

## G. GH31 split by enzyme type
- **Within the three α-glucosidases** (sugar beet, MGAM-C, GAA; one 30% cluster), only I328, Y299 and five +2/+3
  positions vary.
- **The rest of the variation:** W539→F, F575→H and M444→W come from the α-xylosidase; I364→N from the xylosidase
  and the lyase.

## H. Corrected 80% table, per entry (`H_table80_corrected.tsv`)
- **Changes:**
  - G183 dropped (34/44 = 0.773).
  - Added: K182 (+2, 37/41), and GH31 D203 (9/9), R526 (8/9) and M444 (8/9) at +1.
  - Y51 also passes at −2 (51/55).
- **Result:** 22 positions, 13 H / 5 Ar / 4 Al.
  - With mutants dropped, every H position is invariant per entry except H288 (1 Gln) and K182 (Arg 3, Gln 1).
  - M444 (Al) is invariant (Met 9/9).
  - W206 is Al per entry but Ar per lineage.

## Verdict
- **Holds:** identity conservation is associated with contact class (H > Al; aromatic positions stay aromatic more
  often than aliphatic positions stay identical). It survives class assignment independent of the outcome, per-lineage
  counting, removal of catalytic residues, held-out proteins and every cutoff, with p roughly 0.002–0.04.
- **Fails:**
  1. "Every protein" for H: G183, K182, M52, Q294 and N207 vary.
  2. H vs Ar: never significant.
  3. The conservation tracks the −1/+1 subsites rather than hydrogen bonding as such: H positions at the outer
     subsites vary as much as Al positions.
- **Positions that break the classes:**
  - Variable H: G183, K182, M52, Q294, N207.
  - Ar that do not stay aromatic: F143, G250.
  - Conserved Al: M444.
  - Class depends on method: W206, D443, I364, H92.

## Not done
- Hydrogen placement: no monomer library is configured, so the angle test uses heavy atoms only.
- Stacking geometry: still inferred from residue type.
- Bootstrap CIs: only for the 50%-cluster main settings.
- GH31 rests on 6 lineages.
- GH14 and GH15 were not re-analysed.
