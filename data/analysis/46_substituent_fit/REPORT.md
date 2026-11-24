# 46 - Steric fit of starch substituents at each glucose oxygen and subsite (GH13, GH31, GH14, GH15)

Script: `scripts/analysis/a46_substituent_fit.py` (run: `~/miniconda3/bin/python scripts/analysis/a46_substituent_fit.py`;
`... post` rebuilds the tables from `_dd.pkl`; `... spot` rebuilds the three spot-check poses and images).
Environment: ~/miniconda3/bin/python (gemmi 0.7.5, RDKit 2025.03.5). Laptop only, 734 s on 7 processes.
Run: 2026-11-24. Nothing was committed. No manuscript or docs/ file was touched.

## 1. The problem

Analysis 45 says, for every subsite x oxygen cell, whether a protein N/O lies within 3.5 Å of the glucose oxygen
(held) and whether that partner is lineage-invariant (INV), variable (VAR), main-chain only (MC), or whether the
oxygen is free and exposed (EXP) or free and buried (BUR). That does not say whether a substituent on that oxygen
would fit. A held oxygen can still have room on its far side, and a free oxygen can face a wall. The question here is
purely geometric: with the deposited protein and sugar coordinates fixed, is there any orientation of the substituent
that avoids every protein heavy atom and every other sugar atom?

## 2. What the data contain

- 439 density-supported glucose-like units (45: RSCC >= 0.80, subsites -3..+3) in 114 entries, labelled C1..C6,
  O1..O6 from the bonded graph (`45_hydroxyl_map/units_labelled.tsv`). Kinds: Glc (fully glucose), Glc-modified
  (mostly acarbose-type 4-amino-4,6-dideoxy units, no O6), carba (valienamine-type cyclitol), imino (deoxynojirimycin).
- Sites tested: O2 at 439 units, O3 at 438 (TUP's O3 is F and is skipped), deposited O6 at 350 (88 units have no
  O6 and SGD's O6 is S). Modelled O6 rotamers gg, gt, tg at all 439 units (1317 sites). In total 2544 unit x
  oxygen sites x 9 groups (OH hydrogen control + 8 substituents) = 22 896 rows (`site_results.tsv`).
- 50 % sequence groups: `cl50` of `45_hydroxyl_map/atoms.tsv` (from `44_rule_test/clusters_bound.tsv` for
  GH13/GH31 and `45_hydroxyl_map/clusters_GH14_GH15.tsv`). GH14: 2 groups (soybean, *B. cereus*). GH15: 2 groups
  (*Aspergillus*/*Hormoconis* glucoamylases 1GAH, 1GAI, 6FHW; *Saccharomycopsis fibuligera* 2F6D), so every GH15
  fraction rests on two votes. GH15 has only acarbose-type units at -1 and +1, so the glc_units set has no GH15 -1/+1 cells.
- States: `45_hydroxyl_map/cell_summary.tsv` (sets all_glucoselike and glc_units).

## 3. The method

**Substituents** (heavy atoms; as attached to the glucose O):

| name | attached group | heavy atoms beyond O | SMILES used (Cprev-Cring-O-R) | stereo |
|---|---|---|---|---|
| methyl (reference) | O-CH3 | 1 | CCOC | - |
| 2-hydroxyethyl (added, for the HES literature) | O-CH2-CH2-OH | 3 | CCOCCO | - |
| 2-hydroxypropyl | O-CH2-CH(OH)-CH3 | 4 | CCOC[C@H](O)C + CCOC[C@@H](O)C | both (racemic) |
| acetyl | O-C(=O)-CH3 | 3 | CCOC(C)=O | - |
| carboxymethyl | O-CH2-COO⁻ | 4 | CCOCC(=O)[O-] | - |
| phosphate monoester | O-PO3²⁻ | 4 | CCOP(=O)([O-])[O-] | - |
| succinyl (lower bound for OSA) | O-C(=O)-CH2-CH2-COO⁻ | 7 | CCOC(=O)CCC(=O)[O-] | - |
| octenylsuccinyl | O-C(=O)-CH2-CH(COO⁻)-CH2-CH=CH-C5H11 | 15 | CCOC(=O)C[C@H](C(=O)[O-])C/C=C/CCCCC + @@ | both at the branch carbon; E alkene |

**Steps.**
1. Conformers: RDKit ETKDGv3 (seed 4646, RMS pruning 0.05 Å), each MMFF94-optimised (all converged), hydrogens
   dropped. Requested 30-400 per SMILES; kept after pruning: methyl 9, hydroxyethyl 86, hydroxypropyl 91 + 92, acetyl 16,
   carboxymethyl 67, phosphate 76, succinyl 315, octenylsuccinyl 400 + 400 (`fragment_library.tsv`). The small
   fragments have few distinct conformers, so pruning leaves fewer than 300; the rotation in step 3 multiplies each by 12.
   MMFF geometry at the attachment: C-O-X 113-116°, O-C(ether) 1.42-1.44 Å, O-C(ester) 1.36 Å, O-P 1.64 Å.
2. Placement: the fragment's Cprev, Cring, O are matched to the deposited atoms through the local frame at O
   (origin O; e1 along O→Cring; e2 toward Cprev in the Cprev-Cring-O plane; e3 = e1 × e2). The deposited O becomes
   the substituent's attachment oxygen. The fragment's own C-O-X angle and internal geometry are kept exactly.
   Cprev = C1 for O2, C2 for O3, C5 for O6.
3. Systematic Cring-O torsion: every placed conformer is rotated about the Cring-O bond in 30° steps (12 values).
   Poses per site: methyl 108, HE 1032, HP 2196, Ac 192, CM 804, P 912, Suc 3780, OSA 9600.
4. O6 rotamers: O6 is rebuilt on C6 at ω(O5-C5-C6-O6) = -60° (gg), +60° (gt), 180° (tg), with C6-O6 1.43 Å and
   C5-C6-O6 111°. The ω values were checked back from the built coordinates (deviation 0.0°). In carba units "O5" is
   the ring carbon C7. The deposited O6 is removed from the sugar set for these poses. **O6\*** = fits at any of the
   deposited O6, gg, gt, tg. This is the only O6 test possible where no O6 was deposited (acarbose-type units at
   +1 in GH13/GH31/GH15).
5. Distances per pose: the minimum distance from any substituent heavy atom to (a) protein heavy atoms (all amino-acid
   residues of the asymmetric unit, all chains, waters removed, symmetry mates not generated) and (b) sugar heavy
   atoms. Sugar atoms are the covalently connected non-protein, non-water component that contains the unit. From it
   are excluded the atoms within three bonds of the first substituent atom (the O itself, Cring, and the ring
   neighbours of Cring; for O6, C5) and any other residue with >= 3 atoms within 1.2 Å of the unit (overlapping
   alternative copies; see check (a)).
6. Fit criteria (one pose suffices):
   - **flat cut** (main): d_prot >= 3.0 Å and d_sugar >= 2.8 Å; also d_prot cut 2.8 and 3.2 Å.
   - **polar-aware** (added after the first run, see check (a)): a substituent O/N may approach a protein or sugar N/O
     to 2.5 Å (hydrogen-bond or salt-bridge distance); every other pair needs >= 3.0 Å (protein) or >= 2.8 Å (sugar).
7. Aggregation (as 45): a unit fits if any pose fits. An entry fits if >= half of its units in that cell fit. A
   50 % group fits if >= half of its entries fit. Category per cell: **F** = fits in >= 50 % of groups; **s** = fits in
   at least one unit but in < 50 % of groups; **N** = fits in no unit.

The one assumption: the deposited protein and sugar coordinates are where a substituted chain would sit. A pose
that fits is therefore possible only without protein or chain movement, and a pose that does not fit needs some movement.

## 4. What is read off

Quantities per unit x oxygen x substituent (`site_results.tsv`): `fits_<cut>` (any pose fits), `frac_<cut>` (share of
sampled poses that fit, which depends on the pose set and is not a population), `best_clear_prot` (largest d_prot
over poses that pass the sugar test; capped at 8 Å), `o_min_prot` (distance of the oxygen itself to the nearest
protein heavy atom), `held` (45). Per cell (`cell_fit_summary.tsv`): n units, entries, groups; fractions fitting;
category; 45 state.

### 4.1 Summary: fit category against 45 state (all_glucoselike, 3.0 Å flat, cells O2, O3, O6*)

Share of family x subsite x oxygen cells in which the substituent fits in >= 50 % of groups (F). Cells per state:
INV 21, VAR 7, MC 6, EXP 15, BUR 3, no 45 state (O6* where no O6 is deposited) 2.

| 45 state | Me | HE | HP | Ac | CM | P | Suc | OSA |
|---|---|---|---|---|---|---|---|---|
| INV | 0.38 | 0.19 | 0.14 | 0.05 | 0.19 | 0.00 | 0.05 | 0.05 |
| VAR | 0.43 | 0.43 | 0.43 | 0.14 | 0.14 | 0.00 | 0.14 | 0.00 |
| MC | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 | 0.17 | 0.33 | 0.00 |
| EXP | 1.00 | 1.00 | 1.00 | 0.93 | 1.00 | 0.67 | 0.87 | 0.80 |
| BUR | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.33 | 1.00 | 1.00 |

Facts: at the unit level, methyl fits at 94 % of free oxygens (held = False, n = 344) and at 29 % of held oxygens
(n = 883) (deposited O2/O3/O6, all families). By cell state: EXP 0.97 (n 270), BUR 0.89 (19), MC 0.73 (73), VAR 0.44
(349), INV 0.19 (516). Phosphate fits in no INV or VAR cell at the >= 50 %-of-groups level.

### 4.2 Combined tables

The combined tables (45 state in bold, then the category for Me HE HP Ac CM P Suc OSA) and the per-substituent tables
with numbers (fraction of groups, n fit / n groups, fraction of entries / n entries) are in the generated sections below
("Fit tables", "Combined table"), for both sets and both criteria. The cells that read most directly against the literature:

| cell | 45 state | Me | HP | Ac | CM | P | Suc | OSA |
|---|---|---|---|---|---|---|---|---|
| GH13 -1 O2 / O3 / O6* | INV | s / s / s | s / s / s | N / N / s | s / N / s | N / N / N | N / N / s | N / N / N |
| GH13 +1 O3 | INV | s | N | N | N | N | N | N |
| GH13 +1 O6* | (BUR at the 3 deposited O6) | F | F | F | F | F (9/17 groups) | F | F |
| GH13 -2 O2 / O3 | VAR | s / F | s / F | s / s | s / s | N / s | s / s | s / s |
| GH13 +2 O6 (dep.) | EXP | F | F | F | F | F (10/15) | F | F |
| GH13 -3, +3 O6 | EXP | F | F | F | F | F | F | F |
| GH31 -1 all O | INV | O3, O6 F; O2 s | N | N | N | N | N | N |
| GH31 +1 O3 | INV | N | N | N | N | N | N | N |
| GH14 -2 O6, +1 O3 | INV | N | N | N | N | N | N | N |
| GH14 -1 O2 | MC | F | F | F | F | F | F | s |
| GH15 -1 O2 | INV | N | N | N | N | N | N | N |
| GH15 -1 O6* | INV | F | N | N | N | N | N | N |
| GH15 +1 O2 | INV | F | F | F | F | N | F | F |
| GH15 +1 O3 | INV | N | N | N | N | N | N | N |
| GH15 +1 O6* | (no O6 deposited) | F | F | N | N | N | N | N |
| GH15 +2 O6 | MC | N | N | N | N | N | N | N |

## 5. The checks

**(a) Control: the hydroxyl itself.** The OH hydrogen (O-H 0.97 Å, Cring-O-H 108.5°, 12 torsions) was tested with
H...heavy >= 1.6 / 1.8 / 2.0 Å. In the first run 21 of 1227 deposited sites failed. All 21 were in 4 entries (5KZW,
5KZX, 7WJB GH31 -1; 1WDR GH14 +2), where an α and a β anomer, or two ligand copies, are deposited as separate
residues at the same place (0.03-0.8 Å apart). The sugar test counted the copy as a clash. Fix: residues with >= 3
atoms within 1.2 Å of the unit are dropped from the sugar set (7 units). After the fix the OH hydrogen fits at 1227
of 1227 sites at all three cuts. The oxygen itself is within 3.0 Å of a protein heavy atom at many sites (GH13 Glc:
-1 O2 91 %, O3 94 %; +2 O3 87 %; table under "Check (a)" below). This is why the test starts at the first substituent atom and
not at O. Methyl in GH13 Glc units fits at 95-100 % of units at -3 (EXP) and +3 and at 98 % at +2 O6 (EXP). It fits
at 6-9 % at -1 (INV), where the oxygen is held. The control therefore behaves as the task predicts: methyl fits
where the oxygen is exposed.

The first run also showed that the flat 3.0 Å rule counts an O...N/O contact at 2.5-3.0 Å (a hydrogen bond, a
phosphate-arginine salt bridge) as a clash. The polar-aware criterion was added for that reason. It changes no
methyl or HE cell and 2-7 cells per polar substituent (table under "Check (b)" below).

**(b) Clash cut-off.** Cells (O2, O3, O6, O6*; all_glucoselike; 70 per substituent) whose category changes relative
to 3.0 Å: 2.8 Å 4-13 cells (most for methyl, 13, and acetyl, 9); 3.2 Å 2-6 cells; polar-aware 0-7. Mean unit fit
fraction 2.8 / 3.0 / 3.2 Å: methyl 0.71 / 0.61 / 0.54; phosphate 0.19 / 0.16 / 0.14; OSA 0.32 / 0.28 / 0.25. Full
table under "Check (b)" below. In the literature comparison, 4 of 36 rows have a category at 2.8 or 3.2 Å that
differs from 3.0 Å: GH13 -1 O3 HE (s; N at 3.2), GH13 +1 O6* P (F; s at 3.2), GH13 +1 O2 Me (s; F at 2.8),
GH15 -1 O3 Me (s; F at 2.8, N at 3.2).

**(c) Visual spot checks** (`spot/*.pdb`, `*.pml`, `*.png`). Each PDB holds the protein and hetero residues within
10 Å of the oxygen plus the best pose (resname SUB, chain Z; the modelled O6 is resname MOD):

| cell | entry, unit | result | best pose: closest protein contact |
|---|---|---|---|
| GH13 -1 O3, phosphate | 5TD4 GLC | no pose fits (0/912) | phosphate O on Tyr62 CE2 at 1.44 Å; the best pose still lies on the Tyr62 ring |
| GH13 +2 O6, phosphate | 1MFU GLC | fits; 25 % of poses | 4.75 Å (Ser163 OG); the phosphate points into solvent |
| GH15 +1 O6gt, phosphate | 1GAH AC1 | no pose fits; best 2.44 Å | phosphate O on Trp52 CH2 at 2.44 Å (an apolar contact, so the polar-aware rule does not rescue it) |

The images agree with the numbers. They are quick renderings for checking, not figures.

## 6. Literature comparison (summary)

Rows from `substituted_starch_hydrolysis_literature.md` §7, §8, §11 (36 rows; full table under "Literature comparison"
below, also `literature_comparison.tsv`). Agreement rule (mechanical): literature tolerated ↔ F; not tolerated ↔ N; s, or a
mixed or partial literature verdict, gives "partial". Flat 3.0 Å: 16 yes, 15 partial, 5 no. Polar-aware: 17 yes, 14
partial, 5 no (GH13 -3 O3 phosphate moves s→F and becomes yes).

Agree (yes): GH13 -1 O3 and O6* phosphate (N); GH13 +1 O3 phosphate (N); GH13 +2 O6* phosphate (F, tolerated by
*B. subtilis*); GH13 -3 O6* phosphate (F); GH14 -1, -2 O6* phosphate (N); GH14 +1, -2 O3 phosphate (N); GH14 +1 O2
methyl (N); GH15 -1 O6* and O3 phosphate (N); GH15 +1 O3 phosphate (N); GH15 +2 O3 phosphate (F); GH15 -1 O2 methyl
(N); GH15 +1 O2 methyl (F; Lemieux 1996: O-methyl at +1 O2 tolerated).

Do not agree (no):
- **GH15 +1 O6*, phosphate**: literature tolerated (*A. niger* GA). Here no pose fits in any of the 4 entries at
  any rotamer. Best clearances gg 1.3-1.8, gt 2.4-2.6, tg ≤ 1.6 Å, set by Trp52 (1GAH). There is no deposited O6 at
  +1 (acarbose 6-deoxy unit), so the O6 position itself is modelled. Methyl and HP fit at gt.
- **GH14 +1 O6*, phosphate**: literature tolerated (sweet-potato β-amylase). Here N in 3 entries (1J18, 1WDR, 8ZRZ), best
  clearance per entry 2.1, 2.5, 2.2 Å. Methyl, HE, HP and CM fit (F). The structures are soybean and *B. cereus* β-amylase, not
  sweet potato.
- **GH15 -1 O6*, methyl**: literature blocks (Lemieux 1996). Here F (2/2 groups, 3/4 entries) with 3-14 % of poses
  fitting, best clearance 3.0-3.4 Å, on carba (acarbose valienamine) units whose O6 is held (45 INV).
- **GH14 -1 O2, methyl**: literature not tolerated at -2/-1 (van der Burgt 2000b; substituent position not resolved
  in that study). Here F, 18 of 19 units, best 3.2-4.4 Å; held by main chain (45 MC).
- **GH13 +1 O6*, methyl**: literature >90 % of terminal residues of *B. subtilis* limit dextrins unsubstituted
  (position mixed). Here F 17/17 groups; at +1 the O6 is free (45 BUR).

In the "no" rows where the substituent fits sterically, the oxygen is held in 45 (GH15 -1 O6, GH14 -1 O2) or the
literature position is unresolved. Steric fit does not include the hydrogen bond lost when the hydroxyl hydrogen is
replaced ("What it does not include").

Partial rows include all the GH13 -2/+2/+3 phosphate cells where enzymes within GH13 disagree in the literature
(PPA, *B. subtilis*, TAKA). In this map those cells are s (a minority of groups fit) except +2 O3 (N at 3.0 Å, s
polar-aware).

## 7. What it does not include

- Static protein and static sugar chain: no induced fit, no side-chain rotamer change, no shift of the chain
  register. A substituent that clashes by 0.5 Å with a flexible side chain is scored N.
- Substituent conformers only: no energies, no solvation, no electrostatics beyond the polar-aware distance rule. A
  fitting pose may be strained or desolvated. `frac_*` is a share of the sampled poses, not a population.
- Steric fit is not tolerance. Replacing the hydroxyl hydrogen removes its hydrogen-bond donation. Where 45 marks the
  oxygen held (INV, VAR, MC) a fitting substituent still loses that bond, and the oxygen's acceptor role remains only
  for ethers. The `held` share is reported beside every fit.
- Inhibitor units dominate -1 and +1: GH13 -1 has 33 Glc of 82 units (46 carba), GH13 +1 has 3 Glc of 50 (47
  acarbose-type), GH31 -1 has 8 Glc of 32, and GH31 +1 and GH15 -1/+1 have no Glc unit at all. Acarbose units lack O6,
  so O6 at +1 in GH13, GH31 and GH15 is tested only on modelled rotamers (no density for O6 at +1 in GH31/GH15).
- Asymmetric unit only: symmetry mates were not generated. A crystal contact near an exposed oxygen is missed, and
  a substituent pointing into the lattice can score as fitting.
- Hydrogens are absent from protein and sugar. The 3.0 Å heavy-atom cut is permissive for C...C (sum of van der Waals
  radii ≈ 3.4 Å) and strict for O...N hydrogen bonds; hence the 2.8/3.2 Å and polar-aware variants.
- GH14 and GH15 rest on two 50 % groups each. GH15 +1/-1 are four acarbose complexes.
- Literature cells are from other enzymes (porcine, *B. subtilis*, *A. oryzae* α-amylases; sweet-potato β-amylase;
  *A. niger* glucoamylase) than many of the structures. Hydroxyethyl (HES literature) and the phosphate monoester
  were modelled as given; mixed-position literature data (methyl, HE) were compared with each oxygen separately.

## Files

- `scripts/analysis/a46_substituent_fit.py`
- `resultados/analysis/46_substituent_fit/site_results.tsv`: one row per unit x oxygen site x substituent
- `.../cell_fit_summary.tsv`: per set x criterion x family x subsite x oxygen x substituent
- `.../literature_comparison.tsv`, `.../fragment_library.tsv`, `.../run.log`
- `.../spot/` three PDB + PyMOL script + PNG, `spot_checks.tsv`
- `.../_dd.pkl`: intermediate (site rows + O6* rows) used by `post`

---

The sections below are generated by the script.


## Fragment library (RDKit ETKDGv3 + MMFF94, heavy atoms only; each conformer x 12 rotations about Cring-O)

| substituent     | smiles                                |   conformers |   mmff_converged |   n_heavy_substituent |   C_O_X_angle_last |   O_X_bond_last |
|:----------------|:--------------------------------------|-------------:|-----------------:|----------------------:|-------------------:|----------------:|
| methyl          | CCOC                                  |            9 |                9 |                     1 |              113.3 |           1.421 |
| hydroxyethyl    | CCOCCO                                |           86 |               86 |                     3 |              114   |           1.434 |
| hydroxypropyl   | CCOC[C@H](O)C                         |           91 |               91 |                     4 |              113.3 |           1.426 |
| hydroxypropyl   | CCOC[C@@H](O)C                        |           92 |               92 |                     4 |              113.5 |           1.434 |
| acetyl          | CCOC(C)=O                             |           16 |               16 |                     3 |              115.7 |           1.361 |
| carboxymethyl   | CCOCC(=O)[O-]                         |           67 |               67 |                     4 |              113   |           1.444 |
| phosphate       | CCOP(=O)([O-])[O-]                    |           76 |               76 |                     4 |              115   |           1.643 |
| succinyl        | CCOC(=O)CCC(=O)[O-]                   |          315 |              315 |                     7 |              115.9 |           1.364 |
| octenylsuccinyl | CCOC(=O)C[C@H](C(=O)[O-])C/C=C/CCCCC  |          400 |              400 |                    15 |              115.2 |           1.362 |
| octenylsuccinyl | CCOC(=O)C[C@@H](C(=O)[O-])C/C=C/CCCCC |          400 |              400 |                    15 |              115.5 |           1.36  |

Pose counts per site: OH_H 12, methyl 108, hydroxyethyl 1032, hydroxypropyl 2196, acetyl 192, carboxymethyl 804, phosphate 912, succinyl 3780, octenylsuccinyl 9600.


## Fit tables at 3.0 A, set all_glucoselike

Cell: category (F = fits in >= 50 % of 50 % groups; s = fits in some unit but < 50 % of groups; N = fits in no unit), fraction of groups (n fit / n groups; fraction of entries / n entries). O6 = deposited O6 position; O6* = any of deposited, gg, gt, tg. 45 state in the first column of each oxygen.


### GH13 - methyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 0.82 (14/17; ent 0.92/38) | EXP | F 0.88 (15/17; ent 0.95/38) | EXP | F 1.00 (10/10; ent 0.95/22) | EXP | F 0.88 (15/17; ent 0.92/38) |
| -2 | VAR | s 0.36 (8/22; ent 0.46/54) | VAR | F 0.73 (16/22; ent 0.81/54) | VAR | F 0.67 (14/21; ent 0.47/53) | VAR | F 0.86 (19/22; ent 0.67/54) |
| -1 | INV | s 0.03 (1/29; ent 0.04/70) | INV | s 0.00 (0/29; ent 0.03/70) | INV | s 0.03 (1/29; ent 0.04/70) | INV | s 0.48 (14/29; ent 0.46/70) |
| +1 | VAR | s 0.35 (6/17; ent 0.28/43) | INV | s 0.06 (1/17; ent 0.09/43) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 1.00 (17/17; ent 0.98/43) |
| +2 | VAR | F 0.80 (12/15; ent 0.41/41) | VAR | s 0.47 (7/15; ent 0.24/41) | EXP | F 0.93 (14/15; ent 0.98/41) | EXP | F 1.00 (15/15; ent 1.00/41) |
| +3 | MC | F 1.00 (10/10; ent 1.00/13) | MC | F 1.00 (10/10; ent 1.00/13) | EXP | F 1.00 (10/10; ent 1.00/13) | EXP | F 1.00 (10/10; ent 1.00/13) |

### GH13 - hydroxyethyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 0.77 (13/17; ent 0.90/38) | EXP | F 0.88 (15/17; ent 0.95/38) | EXP | F 1.00 (10/10; ent 0.95/22) | EXP | F 0.88 (15/17; ent 0.90/38) |
| -2 | VAR | s 0.27 (6/22; ent 0.32/54) | VAR | F 0.64 (14/22; ent 0.76/54) | VAR | s 0.38 (8/21; ent 0.26/53) | VAR | F 0.64 (14/22; ent 0.44/54) |
| -1 | INV | s 0.00 (0/29; ent 0.01/70) | INV | s 0.00 (0/29; ent 0.01/70) | INV | s 0.03 (1/29; ent 0.03/70) | INV | s 0.28 (8/29; ent 0.29/70) |
| +1 | VAR | s 0.06 (1/17; ent 0.05/43) | INV | s 0.06 (1/17; ent 0.02/43) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 0.88 (15/17; ent 0.84/43) |
| +2 | VAR | F 0.60 (9/15; ent 0.32/41) | VAR | s 0.40 (6/15; ent 0.22/41) | EXP | F 0.87 (13/15; ent 0.95/41) | EXP | F 0.93 (14/15; ent 0.98/41) |
| +3 | MC | F 1.00 (10/10; ent 1.00/13) | MC | F 0.90 (9/10; ent 0.77/13) | EXP | F 1.00 (10/10; ent 1.00/13) | EXP | F 1.00 (10/10; ent 1.00/13) |

### GH13 - hydroxypropyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 0.77 (13/17; ent 0.90/38) | EXP | F 0.88 (15/17; ent 0.95/38) | EXP | F 1.00 (10/10; ent 0.95/22) | EXP | F 0.88 (15/17; ent 0.84/38) |
| -2 | VAR | s 0.04 (1/22; ent 0.18/54) | VAR | F 0.50 (11/22; ent 0.61/54) | VAR | s 0.29 (6/21; ent 0.19/53) | VAR | F 0.59 (13/22; ent 0.39/54) |
| -1 | INV | s 0.00 (0/29; ent 0.01/70) | INV | s 0.00 (0/29; ent 0.01/70) | INV | s 0.03 (1/29; ent 0.03/70) | INV | s 0.17 (5/29; ent 0.13/70) |
| +1 | VAR | N 0.00 (0/17; ent 0.00/43) | INV | N 0.00 (0/17; ent 0.00/43) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 0.71 (12/17; ent 0.77/43) |
| +2 | VAR | F 0.60 (9/15; ent 0.32/41) | VAR | s 0.27 (4/15; ent 0.17/41) | EXP | F 0.87 (13/15; ent 0.95/41) | EXP | F 0.93 (14/15; ent 0.98/41) |
| +3 | MC | F 1.00 (10/10; ent 1.00/13) | MC | F 0.80 (8/10; ent 0.69/13) | EXP | F 0.80 (8/10; ent 0.85/13) | EXP | F 0.80 (8/10; ent 0.85/13) |

### GH13 - acetyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 0.71 (12/17; ent 0.79/38) | EXP | F 0.82 (14/17; ent 0.92/38) | EXP | F 1.00 (10/10; ent 0.86/22) | EXP | F 0.88 (15/17; ent 0.79/38) |
| -2 | VAR | s 0.04 (1/22; ent 0.15/54) | VAR | s 0.46 (10/22; ent 0.59/54) | VAR | s 0.14 (3/21; ent 0.07/53) | VAR | s 0.41 (9/22; ent 0.26/54) |
| -1 | INV | N 0.00 (0/29; ent 0.00/70) | INV | N 0.00 (0/29; ent 0.00/70) | INV | s 0.03 (1/29; ent 0.01/70) | INV | s 0.07 (2/29; ent 0.03/70) |
| +1 | VAR | s 0.06 (1/17; ent 0.02/43) | INV | N 0.00 (0/17; ent 0.00/43) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 0.77 (13/17; ent 0.72/43) |
| +2 | VAR | F 0.53 (8/15; ent 0.29/41) | VAR | s 0.13 (2/15; ent 0.05/41) | EXP | F 0.80 (12/15; ent 0.88/41) | EXP | F 0.93 (14/15; ent 0.98/41) |
| +3 | MC | F 0.50 (5/10; ent 0.54/13) | MC | F 0.60 (6/10; ent 0.46/13) | EXP | F 0.70 (7/10; ent 0.77/13) | EXP | F 0.80 (8/10; ent 0.85/13) |

### GH13 - carboxymethyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 0.71 (12/17; ent 0.79/38) | EXP | F 0.82 (14/17; ent 0.92/38) | EXP | F 1.00 (10/10; ent 0.91/22) | EXP | F 0.88 (15/17; ent 0.82/38) |
| -2 | VAR | s 0.09 (2/22; ent 0.15/54) | VAR | s 0.41 (9/22; ent 0.56/54) | VAR | s 0.29 (6/21; ent 0.19/53) | VAR | F 0.59 (13/22; ent 0.37/54) |
| -1 | INV | s 0.00 (0/29; ent 0.01/70) | INV | N 0.00 (0/29; ent 0.00/70) | INV | s 0.03 (1/29; ent 0.03/70) | INV | s 0.14 (4/29; ent 0.07/70) |
| +1 | VAR | N 0.00 (0/17; ent 0.00/43) | INV | N 0.00 (0/17; ent 0.00/43) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 0.71 (12/17; ent 0.70/43) |
| +2 | VAR | s 0.47 (7/15; ent 0.27/41) | VAR | s 0.20 (3/15; ent 0.15/41) | EXP | F 0.87 (13/15; ent 0.95/41) | EXP | F 0.93 (14/15; ent 0.98/41) |
| +3 | MC | F 1.00 (10/10; ent 1.00/13) | MC | F 0.50 (5/10; ent 0.46/13) | EXP | F 0.80 (8/10; ent 0.77/13) | EXP | F 0.80 (8/10; ent 0.77/13) |

### GH13 - phosphate (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | s 0.29 (5/17; ent 0.21/38) | EXP | s 0.47 (8/17; ent 0.42/38) | EXP | F 0.60 (6/10; ent 0.41/22) | EXP | F 0.53 (9/17; ent 0.34/38) |
| -2 | VAR | N 0.00 (0/22; ent 0.00/54) | VAR | s 0.09 (2/22; ent 0.22/54) | VAR | s 0.00 (0/21; ent 0.02/53) | VAR | s 0.09 (2/22; ent 0.07/54) |
| -1 | INV | N 0.00 (0/29; ent 0.00/70) | INV | N 0.00 (0/29; ent 0.00/70) | INV | N 0.00 (0/29; ent 0.00/70) | INV | N 0.00 (0/29; ent 0.00/70) |
| +1 | VAR | N 0.00 (0/17; ent 0.00/43) | INV | N 0.00 (0/17; ent 0.00/43) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 0.53 (9/17; ent 0.28/43) |
| +2 | VAR | s 0.20 (3/15; ent 0.10/41) | VAR | N 0.00 (0/15; ent 0.00/41) | EXP | F 0.60 (9/15; ent 0.78/41) | EXP | F 0.67 (10/15; ent 0.83/41) |
| +3 | MC | s 0.40 (4/10; ent 0.39/13) | MC | s 0.20 (2/10; ent 0.15/13) | EXP | F 0.50 (5/10; ent 0.54/13) | EXP | F 0.60 (6/10; ent 0.61/13) |

### GH13 - succinyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 0.65 (11/17; ent 0.76/38) | EXP | F 0.88 (15/17; ent 0.95/38) | EXP | F 0.90 (9/10; ent 0.77/22) | EXP | F 0.82 (14/17; ent 0.74/38) |
| -2 | VAR | s 0.00 (0/22; ent 0.06/54) | VAR | s 0.27 (6/22; ent 0.43/54) | VAR | s 0.14 (3/21; ent 0.07/53) | VAR | s 0.41 (9/22; ent 0.26/54) |
| -1 | INV | N 0.00 (0/29; ent 0.00/70) | INV | N 0.00 (0/29; ent 0.00/70) | INV | s 0.03 (1/29; ent 0.01/70) | INV | s 0.03 (1/29; ent 0.01/70) |
| +1 | VAR | N 0.00 (0/17; ent 0.00/43) | INV | N 0.00 (0/17; ent 0.00/43) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 0.59 (10/17; ent 0.63/43) |
| +2 | VAR | F 0.53 (8/15; ent 0.29/41) | VAR | s 0.07 (1/15; ent 0.02/41) | EXP | F 0.80 (12/15; ent 0.90/41) | EXP | F 0.93 (14/15; ent 0.98/41) |
| +3 | MC | F 0.50 (5/10; ent 0.54/13) | MC | s 0.40 (4/10; ent 0.31/13) | EXP | F 0.70 (7/10; ent 0.77/13) | EXP | F 0.80 (8/10; ent 0.85/13) |

### GH13 - octenylsuccinyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 0.65 (11/17; ent 0.76/38) | EXP | F 0.88 (15/17; ent 0.95/38) | EXP | F 0.80 (8/10; ent 0.73/22) | EXP | F 0.82 (14/17; ent 0.71/38) |
| -2 | VAR | s 0.00 (0/22; ent 0.06/54) | VAR | s 0.27 (6/22; ent 0.43/54) | VAR | s 0.14 (3/21; ent 0.07/53) | VAR | s 0.36 (8/22; ent 0.24/54) |
| -1 | INV | N 0.00 (0/29; ent 0.00/70) | INV | N 0.00 (0/29; ent 0.00/70) | INV | N 0.00 (0/29; ent 0.00/70) | INV | N 0.00 (0/29; ent 0.00/70) |
| +1 | VAR | N 0.00 (0/17; ent 0.00/43) | INV | N 0.00 (0/17; ent 0.00/43) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 0.53 (9/17; ent 0.60/43) |
| +2 | VAR | s 0.47 (7/15; ent 0.22/41) | VAR | s 0.07 (1/15; ent 0.02/41) | EXP | F 0.73 (11/15; ent 0.88/41) | EXP | F 0.87 (13/15; ent 0.95/41) |
| +3 | MC | s 0.40 (4/10; ent 0.46/13) | MC | s 0.40 (4/10; ent 0.31/13) | EXP | F 0.70 (7/10; ent 0.69/13) | EXP | F 0.80 (8/10; ent 0.85/13) |

### GH31 - methyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | s 0.00 (0/6; ent 0.06/17) | INV | F 1.00 (6/6; ent 0.88/17) | INV | F 0.50 (3/6; ent 0.53/17) | INV | F 0.67 (4/6; ent 0.65/17) |
| +1 | INV | F 0.60 (3/5; ent 0.78/9) | INV | N 0.00 (0/5; ent 0.00/9) | — | — | — | F 1.00 (5/5; ent 0.78/9) |
| +2 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) |
| +3 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | F 1.00 (5/5; ent 1.00/9) | BUR | F 1.00 (5/5; ent 1.00/9) |

### GH31 - hydroxyethyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) |
| +1 | INV | F 0.60 (3/5; ent 0.78/9) | INV | N 0.00 (0/5; ent 0.00/9) | — | — | — | F 0.60 (3/5; ent 0.67/9) |
| +2 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 0.80 (4/5; ent 0.89/9) |
| +3 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | F 1.00 (5/5; ent 1.00/9) | BUR | F 1.00 (5/5; ent 1.00/9) |

### GH31 - hydroxypropyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) |
| +1 | INV | s 0.40 (2/5; ent 0.67/9) | INV | N 0.00 (0/5; ent 0.00/9) | — | — | — | F 0.60 (3/5; ent 0.67/9) |
| +2 | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 0.80 (4/5; ent 0.67/9) | EXP | F 0.80 (4/5; ent 0.67/9) |
| +3 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | F 1.00 (5/5; ent 0.89/9) | BUR | F 1.00 (5/5; ent 1.00/9) |

### GH31 - acetyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) |
| +1 | INV | s 0.40 (2/5; ent 0.67/9) | INV | N 0.00 (0/5; ent 0.00/9) | — | — | — | s 0.20 (1/5; ent 0.11/9) |
| +2 | EXP | F 0.60 (3/5; ent 0.78/9) | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | s 0.40 (2/5; ent 0.33/9) | EXP | s 0.40 (2/5; ent 0.33/9) |
| +3 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | s 0.20 (1/5; ent 0.22/9) | BUR | F 0.80 (4/5; ent 0.89/9) |

### GH31 - carboxymethyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) |
| +1 | INV | F 0.60 (3/5; ent 0.78/9) | INV | N 0.00 (0/5; ent 0.00/9) | — | — | — | s 0.20 (1/5; ent 0.22/9) |
| +2 | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 0.80 (4/5; ent 0.67/9) | EXP | F 0.80 (4/5; ent 0.67/9) |
| +3 | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | F 0.80 (4/5; ent 0.56/9) | BUR | F 1.00 (5/5; ent 0.78/9) |

### GH31 - phosphate (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) |
| +1 | INV | N 0.00 (0/5; ent 0.00/9) | INV | N 0.00 (0/5; ent 0.00/9) | — | — | — | N 0.00 (0/5; ent 0.00/9) |
| +2 | EXP | F 0.60 (3/5; ent 0.78/9) | EXP | s 0.00 (0/5; ent 0.00/9) | EXP | s 0.20 (1/5; ent 0.11/9) | EXP | s 0.20 (1/5; ent 0.11/9) |
| +3 | EXP | F 0.60 (3/5; ent 0.67/9) | EXP | F 0.80 (4/5; ent 0.44/9) | BUR | s 0.20 (1/5; ent 0.22/9) | BUR | s 0.20 (1/5; ent 0.22/9) |

### GH31 - succinyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) |
| +1 | INV | s 0.20 (1/5; ent 0.11/9) | INV | N 0.00 (0/5; ent 0.00/9) | — | — | — | N 0.00 (0/5; ent 0.00/9) |
| +2 | EXP | F 0.60 (3/5; ent 0.78/9) | EXP | s 0.40 (2/5; ent 0.67/9) | EXP | s 0.40 (2/5; ent 0.33/9) | EXP | s 0.40 (2/5; ent 0.33/9) |
| +3 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | s 0.20 (1/5; ent 0.11/9) | BUR | F 0.60 (3/5; ent 0.33/9) |

### GH31 - octenylsuccinyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) | INV | N 0.00 (0/6; ent 0.00/17) |
| +1 | INV | N 0.00 (0/5; ent 0.00/9) | INV | N 0.00 (0/5; ent 0.00/9) | — | — | — | N 0.00 (0/5; ent 0.00/9) |
| +2 | EXP | s 0.40 (2/5; ent 0.67/9) | EXP | s 0.40 (2/5; ent 0.67/9) | EXP | s 0.20 (1/5; ent 0.11/9) | EXP | s 0.20 (1/5; ent 0.11/9) |
| +3 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | s 0.20 (1/5; ent 0.11/9) | BUR | F 0.60 (3/5; ent 0.33/9) |

### GH14 - methyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | F 1.00 (2/2; ent 1.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 1.00 (2/2; ent 0.94/16) | MC | s 0.00 (0/2; ent 0.40/15) | INV | s 0.00 (0/2; ent 0.00/16) | INV | s 0.00 (0/2; ent 0.06/16) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | F 0.50 (1/2; ent 0.67/3) | INV | F 1.00 (2/2; ent 1.00/3) |
| +2 | INV | F 1.00 (2/2; ent 1.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | F 1.00 (2/2; ent 1.00/3) | BUR | F 1.00 (2/2; ent 1.00/3) |

### GH14 - hydroxyethyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 1.00 (2/2; ent 0.88/16) | MC | s 0.00 (0/2; ent 0.33/15) | INV | N 0.00 (0/2; ent 0.00/16) | INV | s 0.00 (0/2; ent 0.00/16) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | F 0.50 (1/2; ent 0.33/3) | INV | F 1.00 (2/2; ent 0.67/3) |
| +2 | INV | F 1.00 (2/2; ent 1.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | F 1.00 (2/2; ent 1.00/3) | BUR | F 1.00 (2/2; ent 1.00/3) |

### GH14 - hydroxypropyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 1.00 (2/2; ent 0.81/16) | MC | s 0.00 (0/2; ent 0.33/15) | INV | N 0.00 (0/2; ent 0.00/16) | INV | N 0.00 (0/2; ent 0.00/16) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | F 0.50 (1/2; ent 0.33/3) |
| +2 | INV | F 0.50 (1/2; ent 0.67/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | F 0.50 (1/2; ent 0.67/3) | BUR | F 0.50 (1/2; ent 0.67/3) |

### GH14 - acetyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 0.50 (1/2; ent 0.69/16) | MC | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/16) | INV | N 0.00 (0/2; ent 0.00/16) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) |
| +2 | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | N 0.00 (0/2; ent 0.00/3) | BUR | F 0.50 (1/2; ent 0.67/3) |

### GH14 - carboxymethyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 0.50 (1/2; ent 0.62/16) | MC | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/16) | INV | N 0.00 (0/2; ent 0.00/16) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | F 0.50 (1/2; ent 0.33/3) |
| +2 | INV | F 0.50 (1/2; ent 0.67/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | F 0.50 (1/2; ent 0.67/3) | BUR | F 0.50 (1/2; ent 0.67/3) |

### GH14 - phosphate (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 0.50 (1/2; ent 0.38/16) | MC | s 0.00 (0/2; ent 0.07/15) | INV | N 0.00 (0/2; ent 0.00/16) | INV | N 0.00 (0/2; ent 0.00/16) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) |
| +2 | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | N 0.00 (0/2; ent 0.00/3) | BUR | N 0.00 (0/2; ent 0.00/3) |

### GH14 - succinyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 0.50 (1/2; ent 0.69/16) | MC | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/16) | INV | N 0.00 (0/2; ent 0.00/16) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) |
| +2 | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | N 0.00 (0/2; ent 0.00/3) | BUR | F 0.50 (1/2; ent 0.67/3) |

### GH14 - octenylsuccinyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | s 0.00 (0/2; ent 0.19/16) | MC | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/16) | INV | N 0.00 (0/2; ent 0.00/16) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) |
| +2 | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | N 0.00 (0/2; ent 0.00/3) | BUR | F 0.50 (1/2; ent 0.33/3) |

### GH15 - methyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/4) | INV | s 0.00 (0/2; ent 0.25/4) | INV | F 1.00 (2/2; ent 0.75/4) | INV | F 1.00 (2/2; ent 0.75/4) |
| +1 | INV | F 1.00 (2/2; ent 1.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | — | — | — | F 1.00 (2/2; ent 1.00/4) |
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) |

### GH15 - hydroxyethyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) |
| +1 | INV | F 1.00 (2/2; ent 1.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | — | — | — | F 1.00 (2/2; ent 1.00/4) |
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) |

### GH15 - hydroxypropyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) |
| +1 | INV | F 1.00 (2/2; ent 1.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | — | — | — | F 0.50 (1/2; ent 0.75/4) |
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) |

### GH15 - acetyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) |
| +1 | INV | F 1.00 (2/2; ent 1.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | — | — | — | N 0.00 (0/2; ent 0.00/4) |
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 0.50 (1/2; ent 0.75/4) | EXP | F 0.50 (1/2; ent 0.50/4) | EXP | F 0.50 (1/2; ent 0.50/4) |

### GH15 - carboxymethyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) |
| +1 | INV | F 1.00 (2/2; ent 1.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | — | — | — | N 0.00 (0/2; ent 0.00/4) |
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 0.75/4) | EXP | F 1.00 (2/2; ent 1.00/4) |

### GH15 - phosphate (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) |
| +1 | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | — | — | — | N 0.00 (0/2; ent 0.00/4) |
| +2 | EXP | F 0.50 (1/2; ent 0.50/4) | EXP | F 0.50 (1/2; ent 0.50/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 0.50 (1/2; ent 0.50/4) | EXP | s 0.00 (0/2; ent 0.25/4) | EXP | s 0.00 (0/2; ent 0.25/4) |

### GH15 - succinyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) |
| +1 | INV | F 1.00 (2/2; ent 1.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | — | — | — | N 0.00 (0/2; ent 0.00/4) |
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 0.50 (1/2; ent 0.50/4) | EXP | F 0.50 (1/2; ent 0.50/4) |

### GH15 - octenylsuccinyl (all_glucoselike)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | INV | N 0.00 (0/2; ent 0.00/4) |
| +1 | INV | F 1.00 (2/2; ent 1.00/4) | INV | N 0.00 (0/2; ent 0.00/4) | — | — | — | N 0.00 (0/2; ent 0.00/4) |
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 0.50 (1/2; ent 0.50/4) | EXP | F 0.50 (1/2; ent 0.50/4) |

## Fit tables at 3.0 A, set glc_units

Cell: category (F = fits in >= 50 % of 50 % groups; s = fits in some unit but < 50 % of groups; N = fits in no unit), fraction of groups (n fit / n groups; fraction of entries / n entries). O6 = deposited O6 position; O6* = any of deposited, gg, gt, tg. 45 state in the first column of each oxygen.


### GH13 - methyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 1.00 (8/8; ent 1.00/20) | EXP | F 1.00 (8/8; ent 1.00/20) | EXP | F 1.00 (8/8; ent 0.95/20) | EXP | F 1.00 (8/8; ent 1.00/20) |
| -2 | VAR | s 0.38 (8/21; ent 0.47/53) | VAR | F 0.71 (15/21; ent 0.81/53) | VAR | F 0.67 (14/21; ent 0.47/53) | VAR | F 0.86 (18/21; ent 0.66/53) |
| -1 | INV | s 0.06 (1/18; ent 0.11/28) | INV | s 0.06 (1/18; ent 0.07/28) | INV | s 0.11 (2/18; ent 0.11/28) | INV | s 0.44 (8/18; ent 0.39/28) |
| +1 | VAR | F 0.67 (2/3; ent 0.67/3) | INV | s 0.33 (1/3; ent 0.33/3) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 1.00 (3/3; ent 1.00/3) |
| +2 | VAR | F 0.80 (12/15; ent 0.41/41) | VAR | s 0.47 (7/15; ent 0.24/41) | EXP | F 0.93 (14/15; ent 0.98/41) | EXP | F 1.00 (15/15; ent 1.00/41) |
| +3 | MC | F 1.00 (10/10; ent 1.00/13) | MC | F 1.00 (10/10; ent 1.00/13) | EXP | F 1.00 (10/10; ent 1.00/13) | EXP | F 1.00 (10/10; ent 1.00/13) |

### GH13 - hydroxyethyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 0.88 (7/8; ent 0.95/20) | EXP | F 1.00 (8/8; ent 1.00/20) | EXP | F 1.00 (8/8; ent 0.95/20) | EXP | F 1.00 (8/8; ent 0.95/20) |
| -2 | VAR | s 0.29 (6/21; ent 0.32/53) | VAR | F 0.67 (14/21; ent 0.77/53) | VAR | s 0.38 (8/21; ent 0.26/53) | VAR | F 0.62 (13/21; ent 0.43/53) |
| -1 | INV | s 0.00 (0/18; ent 0.04/28) | INV | s 0.00 (0/18; ent 0.04/28) | INV | s 0.06 (1/18; ent 0.07/28) | INV | s 0.28 (5/18; ent 0.21/28) |
| +1 | VAR | N 0.00 (0/3; ent 0.00/3) | INV | s 0.33 (1/3; ent 0.33/3) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 0.67 (2/3; ent 0.67/3) |
| +2 | VAR | F 0.60 (9/15; ent 0.32/41) | VAR | s 0.40 (6/15; ent 0.22/41) | EXP | F 0.87 (13/15; ent 0.95/41) | EXP | F 0.93 (14/15; ent 0.98/41) |
| +3 | MC | F 1.00 (10/10; ent 1.00/13) | MC | F 0.90 (9/10; ent 0.77/13) | EXP | F 1.00 (10/10; ent 1.00/13) | EXP | F 1.00 (10/10; ent 1.00/13) |

### GH13 - hydroxypropyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 0.88 (7/8; ent 0.95/20) | EXP | F 1.00 (8/8; ent 1.00/20) | EXP | F 1.00 (8/8; ent 0.95/20) | EXP | F 1.00 (8/8; ent 0.95/20) |
| -2 | VAR | s 0.05 (1/21; ent 0.19/53) | VAR | F 0.52 (11/21; ent 0.62/53) | VAR | s 0.29 (6/21; ent 0.19/53) | VAR | F 0.57 (12/21; ent 0.38/53) |
| -1 | INV | s 0.00 (0/18; ent 0.04/28) | INV | s 0.00 (0/18; ent 0.04/28) | INV | s 0.06 (1/18; ent 0.07/28) | INV | s 0.17 (3/18; ent 0.14/28) |
| +1 | VAR | N 0.00 (0/3; ent 0.00/3) | INV | N 0.00 (0/3; ent 0.00/3) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 0.67 (2/3; ent 0.67/3) |
| +2 | VAR | F 0.60 (9/15; ent 0.32/41) | VAR | s 0.27 (4/15; ent 0.17/41) | EXP | F 0.87 (13/15; ent 0.95/41) | EXP | F 0.93 (14/15; ent 0.98/41) |
| +3 | MC | F 1.00 (10/10; ent 1.00/13) | MC | F 0.80 (8/10; ent 0.69/13) | EXP | F 0.80 (8/10; ent 0.85/13) | EXP | F 0.80 (8/10; ent 0.85/13) |

### GH13 - acetyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 0.75 (6/8; ent 0.80/20) | EXP | F 1.00 (8/8; ent 1.00/20) | EXP | F 1.00 (8/8; ent 0.85/20) | EXP | F 1.00 (8/8; ent 0.85/20) |
| -2 | VAR | s 0.05 (1/21; ent 0.15/53) | VAR | s 0.48 (10/21; ent 0.60/53) | VAR | s 0.14 (3/21; ent 0.07/53) | VAR | s 0.38 (8/21; ent 0.24/53) |
| -1 | INV | N 0.00 (0/18; ent 0.00/28) | INV | N 0.00 (0/18; ent 0.00/28) | INV | s 0.06 (1/18; ent 0.04/28) | INV | s 0.11 (2/18; ent 0.07/28) |
| +1 | VAR | N 0.00 (0/3; ent 0.00/3) | INV | N 0.00 (0/3; ent 0.00/3) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 0.67 (2/3; ent 0.67/3) |
| +2 | VAR | F 0.53 (8/15; ent 0.29/41) | VAR | s 0.13 (2/15; ent 0.05/41) | EXP | F 0.80 (12/15; ent 0.88/41) | EXP | F 0.93 (14/15; ent 0.98/41) |
| +3 | MC | F 0.50 (5/10; ent 0.54/13) | MC | F 0.60 (6/10; ent 0.46/13) | EXP | F 0.70 (7/10; ent 0.77/13) | EXP | F 0.80 (8/10; ent 0.85/13) |

### GH13 - carboxymethyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 0.75 (6/8; ent 0.85/20) | EXP | F 1.00 (8/8; ent 1.00/20) | EXP | F 1.00 (8/8; ent 0.90/20) | EXP | F 1.00 (8/8; ent 0.90/20) |
| -2 | VAR | s 0.10 (2/21; ent 0.15/53) | VAR | s 0.43 (9/21; ent 0.57/53) | VAR | s 0.29 (6/21; ent 0.19/53) | VAR | F 0.57 (12/21; ent 0.36/53) |
| -1 | INV | s 0.00 (0/18; ent 0.04/28) | INV | N 0.00 (0/18; ent 0.00/28) | INV | s 0.06 (1/18; ent 0.07/28) | INV | s 0.17 (3/18; ent 0.14/28) |
| +1 | VAR | N 0.00 (0/3; ent 0.00/3) | INV | N 0.00 (0/3; ent 0.00/3) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 0.67 (2/3; ent 0.67/3) |
| +2 | VAR | s 0.47 (7/15; ent 0.27/41) | VAR | s 0.20 (3/15; ent 0.15/41) | EXP | F 0.87 (13/15; ent 0.95/41) | EXP | F 0.93 (14/15; ent 0.98/41) |
| +3 | MC | F 1.00 (10/10; ent 1.00/13) | MC | F 0.50 (5/10; ent 0.46/13) | EXP | F 0.80 (8/10; ent 0.77/13) | EXP | F 0.80 (8/10; ent 0.77/13) |

### GH13 - phosphate (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | s 0.12 (1/8; ent 0.15/20) | EXP | F 0.62 (5/8; ent 0.50/20) | EXP | F 0.50 (4/8; ent 0.35/20) | EXP | F 0.50 (4/8; ent 0.35/20) |
| -2 | VAR | N 0.00 (0/21; ent 0.00/53) | VAR | s 0.10 (2/21; ent 0.23/53) | VAR | s 0.00 (0/21; ent 0.02/53) | VAR | s 0.10 (2/21; ent 0.07/53) |
| -1 | INV | N 0.00 (0/18; ent 0.00/28) | INV | N 0.00 (0/18; ent 0.00/28) | INV | N 0.00 (0/18; ent 0.00/28) | INV | N 0.00 (0/18; ent 0.00/28) |
| +1 | VAR | N 0.00 (0/3; ent 0.00/3) | INV | N 0.00 (0/3; ent 0.00/3) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | s 0.33 (1/3; ent 0.33/3) |
| +2 | VAR | s 0.20 (3/15; ent 0.10/41) | VAR | N 0.00 (0/15; ent 0.00/41) | EXP | F 0.60 (9/15; ent 0.78/41) | EXP | F 0.67 (10/15; ent 0.83/41) |
| +3 | MC | s 0.40 (4/10; ent 0.39/13) | MC | s 0.20 (2/10; ent 0.15/13) | EXP | F 0.50 (5/10; ent 0.54/13) | EXP | F 0.60 (6/10; ent 0.61/13) |

### GH13 - succinyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 0.75 (6/8; ent 0.75/20) | EXP | F 1.00 (8/8; ent 1.00/20) | EXP | F 0.88 (7/8; ent 0.75/20) | EXP | F 0.88 (7/8; ent 0.80/20) |
| -2 | VAR | s 0.00 (0/21; ent 0.06/53) | VAR | s 0.29 (6/21; ent 0.43/53) | VAR | s 0.14 (3/21; ent 0.07/53) | VAR | s 0.38 (8/21; ent 0.24/53) |
| -1 | INV | N 0.00 (0/18; ent 0.00/28) | INV | N 0.00 (0/18; ent 0.00/28) | INV | s 0.06 (1/18; ent 0.04/28) | INV | s 0.06 (1/18; ent 0.04/28) |
| +1 | VAR | N 0.00 (0/3; ent 0.00/3) | INV | N 0.00 (0/3; ent 0.00/3) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 0.67 (2/3; ent 0.67/3) |
| +2 | VAR | F 0.53 (8/15; ent 0.29/41) | VAR | s 0.07 (1/15; ent 0.02/41) | EXP | F 0.80 (12/15; ent 0.90/41) | EXP | F 0.93 (14/15; ent 0.98/41) |
| +3 | MC | F 0.50 (5/10; ent 0.54/13) | MC | s 0.40 (4/10; ent 0.31/13) | EXP | F 0.70 (7/10; ent 0.77/13) | EXP | F 0.80 (8/10; ent 0.85/13) |

### GH13 - octenylsuccinyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -3 | EXP | F 0.75 (6/8; ent 0.75/20) | EXP | F 1.00 (8/8; ent 1.00/20) | EXP | F 0.75 (6/8; ent 0.70/20) | EXP | F 0.88 (7/8; ent 0.80/20) |
| -2 | VAR | s 0.00 (0/21; ent 0.06/53) | VAR | s 0.29 (6/21; ent 0.43/53) | VAR | s 0.14 (3/21; ent 0.07/53) | VAR | s 0.33 (7/21; ent 0.23/53) |
| -1 | INV | N 0.00 (0/18; ent 0.00/28) | INV | N 0.00 (0/18; ent 0.00/28) | INV | N 0.00 (0/18; ent 0.00/28) | INV | N 0.00 (0/18; ent 0.00/28) |
| +1 | VAR | N 0.00 (0/3; ent 0.00/3) | INV | N 0.00 (0/3; ent 0.00/3) | BUR | s 0.33 (1/3; ent 0.33/3) | BUR | F 0.67 (2/3; ent 0.67/3) |
| +2 | VAR | s 0.47 (7/15; ent 0.22/41) | VAR | s 0.07 (1/15; ent 0.02/41) | EXP | F 0.73 (11/15; ent 0.88/41) | EXP | F 0.87 (13/15; ent 0.95/41) |
| +3 | MC | s 0.40 (4/10; ent 0.46/13) | MC | s 0.40 (4/10; ent 0.31/13) | EXP | F 0.70 (7/10; ent 0.69/13) | EXP | F 0.80 (8/10; ent 0.85/13) |

### GH31 - methyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | s 0.00 (0/2; ent 0.20/5) | INV | F 1.00 (2/2; ent 1.00/5) | INV | F 0.50 (1/2; ent 0.40/5) | INV | F 0.50 (1/2; ent 0.40/5) |
| +2 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) |
| +3 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | F 1.00 (5/5; ent 1.00/9) | BUR | F 1.00 (5/5; ent 1.00/9) |

### GH31 - hydroxyethyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) |
| +2 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 0.80 (4/5; ent 0.89/9) |
| +3 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | F 1.00 (5/5; ent 1.00/9) | BUR | F 1.00 (5/5; ent 1.00/9) |

### GH31 - hydroxypropyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) |
| +2 | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 0.80 (4/5; ent 0.67/9) | EXP | F 0.80 (4/5; ent 0.67/9) |
| +3 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | F 1.00 (5/5; ent 0.89/9) | BUR | F 1.00 (5/5; ent 1.00/9) |

### GH31 - acetyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) |
| +2 | EXP | F 0.60 (3/5; ent 0.78/9) | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | s 0.40 (2/5; ent 0.33/9) | EXP | s 0.40 (2/5; ent 0.33/9) |
| +3 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | s 0.20 (1/5; ent 0.22/9) | BUR | F 0.80 (4/5; ent 0.89/9) |

### GH31 - carboxymethyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) |
| +2 | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 0.80 (4/5; ent 0.67/9) | EXP | F 0.80 (4/5; ent 0.67/9) |
| +3 | EXP | F 0.80 (4/5; ent 0.89/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | F 0.80 (4/5; ent 0.56/9) | BUR | F 1.00 (5/5; ent 0.78/9) |

### GH31 - phosphate (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) |
| +2 | EXP | F 0.60 (3/5; ent 0.78/9) | EXP | s 0.00 (0/5; ent 0.00/9) | EXP | s 0.20 (1/5; ent 0.11/9) | EXP | s 0.20 (1/5; ent 0.11/9) |
| +3 | EXP | F 0.60 (3/5; ent 0.67/9) | EXP | F 0.80 (4/5; ent 0.44/9) | BUR | s 0.20 (1/5; ent 0.22/9) | BUR | s 0.20 (1/5; ent 0.22/9) |

### GH31 - succinyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) |
| +2 | EXP | F 0.60 (3/5; ent 0.78/9) | EXP | s 0.40 (2/5; ent 0.67/9) | EXP | s 0.40 (2/5; ent 0.33/9) | EXP | s 0.40 (2/5; ent 0.33/9) |
| +3 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | s 0.20 (1/5; ent 0.11/9) | BUR | F 0.60 (3/5; ent 0.33/9) |

### GH31 - octenylsuccinyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -1 | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) | INV | N 0.00 (0/2; ent 0.00/5) |
| +2 | EXP | s 0.40 (2/5; ent 0.67/9) | EXP | s 0.40 (2/5; ent 0.67/9) | EXP | s 0.20 (1/5; ent 0.11/9) | EXP | s 0.20 (1/5; ent 0.11/9) |
| +3 | EXP | F 1.00 (5/5; ent 1.00/9) | EXP | F 1.00 (5/5; ent 1.00/9) | BUR | s 0.20 (1/5; ent 0.11/9) | BUR | F 0.60 (3/5; ent 0.33/9) |

### GH14 - methyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | F 1.00 (2/2; ent 1.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 1.00 (2/2; ent 0.92/13) | MC | s 0.00 (0/2; ent 0.31/13) | INV | s 0.00 (0/2; ent 0.00/13) | INV | s 0.00 (0/2; ent 0.00/13) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | F 0.50 (1/2; ent 0.67/3) | INV | F 1.00 (2/2; ent 1.00/3) |
| +2 | INV | F 1.00 (2/2; ent 1.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | F 1.00 (2/2; ent 1.00/3) | BUR | F 1.00 (2/2; ent 1.00/3) |

### GH14 - hydroxyethyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 1.00 (2/2; ent 0.85/13) | MC | s 0.00 (0/2; ent 0.23/13) | INV | N 0.00 (0/2; ent 0.00/13) | INV | s 0.00 (0/2; ent 0.00/13) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | F 0.50 (1/2; ent 0.33/3) | INV | F 1.00 (2/2; ent 0.67/3) |
| +2 | INV | F 1.00 (2/2; ent 1.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | F 1.00 (2/2; ent 1.00/3) | BUR | F 1.00 (2/2; ent 1.00/3) |

### GH14 - hydroxypropyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 1.00 (2/2; ent 0.77/13) | MC | s 0.00 (0/2; ent 0.23/13) | INV | N 0.00 (0/2; ent 0.00/13) | INV | N 0.00 (0/2; ent 0.00/13) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | F 0.50 (1/2; ent 0.33/3) |
| +2 | INV | F 0.50 (1/2; ent 0.67/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | F 0.50 (1/2; ent 0.67/3) | BUR | F 0.50 (1/2; ent 0.67/3) |

### GH14 - acetyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 0.50 (1/2; ent 0.61/13) | MC | N 0.00 (0/2; ent 0.00/13) | INV | N 0.00 (0/2; ent 0.00/13) | INV | N 0.00 (0/2; ent 0.00/13) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) |
| +2 | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | N 0.00 (0/2; ent 0.00/3) | BUR | F 0.50 (1/2; ent 0.67/3) |

### GH14 - carboxymethyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 0.50 (1/2; ent 0.54/13) | MC | s 0.00 (0/2; ent 0.00/13) | INV | N 0.00 (0/2; ent 0.00/13) | INV | N 0.00 (0/2; ent 0.00/13) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | F 0.50 (1/2; ent 0.33/3) |
| +2 | INV | F 0.50 (1/2; ent 0.67/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | F 0.50 (1/2; ent 0.67/3) | BUR | F 0.50 (1/2; ent 0.67/3) |

### GH14 - phosphate (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 0.50 (1/2; ent 0.31/13) | MC | N 0.00 (0/2; ent 0.00/13) | INV | N 0.00 (0/2; ent 0.00/13) | INV | N 0.00 (0/2; ent 0.00/13) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) |
| +2 | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | N 0.00 (0/2; ent 0.00/3) | BUR | N 0.00 (0/2; ent 0.00/3) |

### GH14 - succinyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | F 0.50 (1/2; ent 0.61/13) | MC | N 0.00 (0/2; ent 0.00/13) | INV | N 0.00 (0/2; ent 0.00/13) | INV | N 0.00 (0/2; ent 0.00/13) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) |
| +2 | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | N 0.00 (0/2; ent 0.00/3) | BUR | F 0.50 (1/2; ent 0.67/3) |

### GH14 - octenylsuccinyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| -2 | VAR | s 0.00 (0/2; ent 0.13/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) | INV | N 0.00 (0/2; ent 0.00/15) |
| -1 | MC | N 0.00 (0/2; ent 0.00/13) | MC | N 0.00 (0/2; ent 0.00/13) | INV | N 0.00 (0/2; ent 0.00/13) | INV | N 0.00 (0/2; ent 0.00/13) |
| +1 | MC | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) |
| +2 | INV | N 0.00 (0/2; ent 0.00/3) | INV | N 0.00 (0/2; ent 0.00/3) | BUR | N 0.00 (0/2; ent 0.00/3) | BUR | F 0.50 (1/2; ent 0.33/3) |

### GH15 - methyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) |

### GH15 - hydroxyethyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) |

### GH15 - hydroxypropyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) |

### GH15 - acetyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 0.50 (1/2; ent 0.75/4) | EXP | F 0.50 (1/2; ent 0.50/4) | EXP | F 0.50 (1/2; ent 0.50/4) |

### GH15 - carboxymethyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 0.75/4) | EXP | F 1.00 (2/2; ent 1.00/4) |

### GH15 - phosphate (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| +2 | EXP | F 0.50 (1/2; ent 0.50/4) | EXP | F 0.50 (1/2; ent 0.50/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 0.50 (1/2; ent 0.50/4) | EXP | s 0.00 (0/2; ent 0.25/4) | EXP | s 0.00 (0/2; ent 0.25/4) |

### GH15 - succinyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 0.50 (1/2; ent 0.50/4) | EXP | F 0.50 (1/2; ent 0.50/4) |

### GH15 - octenylsuccinyl (glc_units)

| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |
|---|---|---|---|---|---|---|---|---|
| +2 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | MC | N 0.00 (0/2; ent 0.00/4) | MC | N 0.00 (0/2; ent 0.00/4) |
| +3 | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 1.00 (2/2; ent 1.00/4) | EXP | F 0.50 (1/2; ent 0.50/4) | EXP | F 0.50 (1/2; ent 0.50/4) |

## Combined table, 3.0 A flat (all_glucoselike): 45 state, then category per substituent in the order Me HE HP Ac CM P Suc OSA


### GH13

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -3 | **EXP** F F F F F s F F | **EXP** F F F F F s F F | **EXP** F F F F F F F F | **EXP** F F F F F F F F |
| -2 | **VAR** s s s s s N s s | **VAR** F F F s s s s s | **VAR** F s s s s s s s | **VAR** F F F s F s s s |
| -1 | **INV** s s s N s N N N | **INV** s s s N N N N N | **INV** s s s s s N s N | **INV** s s s s s N s N |
| +1 | **VAR** s s N s N N N N | **INV** s s N N N N N N | **BUR** s s s s s s s s | **BUR** F F F F F F F F |
| +2 | **VAR** F F F F s s F s | **VAR** s s s s s N s s | **EXP** F F F F F F F F | **EXP** F F F F F F F F |
| +3 | **MC** F F F F F s F s | **MC** F F F F F s s s | **EXP** F F F F F F F F | **EXP** F F F F F F F F |

### GH31

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -1 | **INV** s N N N N N N N | **INV** F N N N N N N N | **INV** F N N N N N N N | **INV** F N N N N N N N |
| +1 | **INV** F F s s F N s N | **INV** N N N N N N N N | — | **—** F F F s s N N N |
| +2 | **EXP** F F F F F F F s | **EXP** F F F F F s s s | **EXP** F F F s F s s s | **EXP** F F F s F s s s |
| +3 | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **BUR** F F F s F s s s | **BUR** F F F F F s F F |

### GH14

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -2 | **VAR** s s s s s s s s | **INV** F N N N N N N N | **INV** N N N N N N N N | **INV** N N N N N N N N |
| -1 | **MC** F F F F F F F s | **MC** s s s s s s s s | **INV** s N N N N N N N | **INV** s s N N N N N N |
| +1 | **MC** N N N N N N N N | **INV** N N N N N N N N | **INV** F F N N N N N N | **INV** F F F N F N N N |
| +2 | **INV** F F F N F N N N | **INV** N N N N N N N N | **BUR** F F F N F N N N | **BUR** F F F F F N F F |

### GH15

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -1 | **INV** N N N N N N N N | **INV** s N N N N N N N | **INV** F N N N N N N N | **INV** F N N N N N N N |
| +1 | **INV** F F F F F N F F | **INV** N N N N N N N N | — | **—** F F F N N N N N |
| +2 | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **MC** N N N N N N N N | **MC** N N N N N N N N |
| +3 | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **EXP** F F F F F s F F | **EXP** F F F F F s F F |

## Combined table, 3.0 A flat (glc_units): 45 state, then category per substituent in the order Me HE HP Ac CM P Suc OSA


### GH13

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -3 | **EXP** F F F F F s F F | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **EXP** F F F F F F F F |
| -2 | **VAR** s s s s s N s s | **VAR** F F F s s s s s | **VAR** F s s s s s s s | **VAR** F F F s F s s s |
| -1 | **INV** s s s N s N N N | **INV** s s s N N N N N | **INV** s s s s s N s N | **INV** s s s s s N s N |
| +1 | **VAR** F N N N N N N N | **INV** s s N N N N N N | **BUR** s s s s s s s s | **BUR** F F F F F s F F |
| +2 | **VAR** F F F F s s F s | **VAR** s s s s s N s s | **EXP** F F F F F F F F | **EXP** F F F F F F F F |
| +3 | **MC** F F F F F s F s | **MC** F F F F F s s s | **EXP** F F F F F F F F | **EXP** F F F F F F F F |

### GH31

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -1 | **INV** s N N N N N N N | **INV** F N N N N N N N | **INV** F N N N N N N N | **INV** F N N N N N N N |
| +2 | **EXP** F F F F F F F s | **EXP** F F F F F s s s | **EXP** F F F s F s s s | **EXP** F F F s F s s s |
| +3 | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **BUR** F F F s F s s s | **BUR** F F F F F s F F |

### GH14

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -2 | **VAR** s s s s s s s s | **INV** F N N N N N N N | **INV** N N N N N N N N | **INV** N N N N N N N N |
| -1 | **MC** F F F F F F F N | **MC** s s s N s N N N | **INV** s N N N N N N N | **INV** s s N N N N N N |
| +1 | **MC** N N N N N N N N | **INV** N N N N N N N N | **INV** F F N N N N N N | **INV** F F F N F N N N |
| +2 | **INV** F F F N F N N N | **INV** N N N N N N N N | **BUR** F F F N F N N N | **BUR** F F F F F N F F |

### GH15

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| +2 | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **MC** N N N N N N N N | **MC** N N N N N N N N |
| +3 | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **EXP** F F F F F s F F | **EXP** F F F F F s F F |

## Combined table, polar-aware (O/N...N/O >= 2.5 A, other pairs >= 3.0 A) (all_glucoselike): 45 state, then category per substituent in the order Me HE HP Ac CM P Suc OSA


### GH13

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -3 | **EXP** F F F F F s F F | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **EXP** F F F F F F F F |
| -2 | **VAR** s s s s s N s s | **VAR** F F F s F s s s | **VAR** F s s s s s s s | **VAR** F F F F F s F s |
| -1 | **INV** s s s N s N N N | **INV** s s s N s N N N | **INV** s s s s s N s N | **INV** s s s s s N s s |
| +1 | **VAR** s s s s s N s N | **INV** s s N N N N N N | **BUR** s s s s s s s s | **BUR** F F F F F F F F |
| +2 | **VAR** F F F F F s F s | **VAR** s s s s s s s s | **EXP** F F F F F F F F | **EXP** F F F F F F F F |
| +3 | **MC** F F F F F F F s | **MC** F F F F F s s s | **EXP** F F F F F F F F | **EXP** F F F F F F F F |

### GH31

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -1 | **INV** s N N N N N N N | **INV** F N N N N N N N | **INV** F N N N N N N N | **INV** F N N N N N N N |
| +1 | **INV** F F F s F N s N | **INV** N N N N N N N N | — | **—** F F F s s N s s |
| +2 | **EXP** F F F F F F F s | **EXP** F F F F F F s s | **EXP** F F F s F s s s | **EXP** F F F F F s F F |
| +3 | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **BUR** F F F s F s s s | **BUR** F F F F F s F F |

### GH14

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -2 | **VAR** s s s s s s s s | **INV** F N N N N N N N | **INV** N N N N N N N N | **INV** N N N N N N N N |
| -1 | **MC** F F F F F F F F | **MC** s s s s s s s s | **INV** s N N N N N N N | **INV** s s N N N N N N |
| +1 | **MC** N N N N N N N N | **INV** N N N N N N N N | **INV** F F N N N N N N | **INV** F F F N F N N N |
| +2 | **INV** F F F F F N F F | **INV** N N N N N N N N | **BUR** F F F F F N F F | **BUR** F F F F F F F F |

### GH15

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -1 | **INV** N N N N N N N N | **INV** s N N N N N N N | **INV** F N N N N N N N | **INV** F N N N N N N N |
| +1 | **INV** F F F F F N F F | **INV** N N N N N N N N | — | **—** F F F s F N s s |
| +2 | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **MC** N N N N N N N N | **MC** N N N N N N N N |
| +3 | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **EXP** F F F F F F F F |

## Combined table, polar-aware (O/N...N/O >= 2.5 A, other pairs >= 3.0 A) (glc_units): 45 state, then category per substituent in the order Me HE HP Ac CM P Suc OSA


### GH13

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -3 | **EXP** F F F F F s F F | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **EXP** F F F F F F F F |
| -2 | **VAR** s s s s s N s s | **VAR** F F F s F s s s | **VAR** F s s s s s s s | **VAR** F F F s F s s s |
| -1 | **INV** s s s N s N N N | **INV** s s s N s N N N | **INV** s s s s s N s N | **INV** s s s s s N s N |
| +1 | **VAR** F N N N N N N N | **INV** s s N N N N N N | **BUR** s s s s s s s s | **BUR** F F F F F s F F |
| +2 | **VAR** F F F F F s F s | **VAR** s s s s s s s s | **EXP** F F F F F F F F | **EXP** F F F F F F F F |
| +3 | **MC** F F F F F F F s | **MC** F F F F F s s s | **EXP** F F F F F F F F | **EXP** F F F F F F F F |

### GH31

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -1 | **INV** s N N N N N N N | **INV** F N N N N N N N | **INV** F N N N N N N N | **INV** F N N N N N N N |
| +2 | **EXP** F F F F F F F s | **EXP** F F F F F F s s | **EXP** F F F s F s s s | **EXP** F F F F F s F F |
| +3 | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **BUR** F F F s F s s s | **BUR** F F F F F s F F |

### GH14

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| -2 | **VAR** s s s s s s s s | **INV** F N N N N N N N | **INV** N N N N N N N N | **INV** N N N N N N N N |
| -1 | **MC** F F F F F F F F | **MC** s s s N s N N N | **INV** s N N N N N N N | **INV** s s N N N N N N |
| +1 | **MC** N N N N N N N N | **INV** N N N N N N N N | **INV** F F N N N N N N | **INV** F F F N F N N N |
| +2 | **INV** F F F F F N F F | **INV** N N N N N N N N | **BUR** F F F F F N F F | **BUR** F F F F F F F F |

### GH15

| subsite | O2 | O3 | O6 (dep.) | O6* |
|---|---|---|---|---|
| +2 | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **MC** N N N N N N N N | **MC** N N N N N N N N |
| +3 | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **EXP** F F F F F F F F | **EXP** F F F F F F F F |

## Check (a): hydroxyl hydrogen control and methyl, GH13, glucose units (kind Glc)

OH-H fits at 1.8 A (H...protein heavy and H...sugar heavy >= 1.8 A); O itself: median minimum distance of the deposited O to a protein heavy atom and share of units with that distance < 3.0 A; methyl at 3.0 A.

| subsite | label | 45 state | n units | OH-H fits (units) | O-prot median A | O < 3.0 A (units) | methyl fits (units) | methyl F/s/N |
|---|---|---|---|---|---|---|---|---|
| -3 | O2 | EXP | 21 | 1.00 | 4.05 | 0.14 | 0.95 | F |
| -3 | O3 | EXP | 21 | 1.00 | 4.54 | 0.00 | 1.00 | F |
| -3 | O6 | EXP | 21 | 1.00 | 3.87 | 0.19 | 0.95 | F |
| -3 | O6* | EXP | 21 | 1.00 | 2.54 | 0.86 | 1.00 | F |
| -2 | O2 | VAR | 61 | 1.00 | 2.92 | 0.67 | 0.43 | s |
| -2 | O3 | VAR | 61 | 1.00 | 2.94 | 0.54 | 0.80 | F |
| -2 | O6 | VAR | 61 | 1.00 | 2.79 | 0.67 | 0.49 | F |
| -2 | O6* | VAR | 61 | 1.00 | 2.62 | 0.93 | 0.69 | F |
| -1 | O2 | INV | 33 | 1.00 | 2.61 | 0.91 | 0.09 | s |
| -1 | O3 | INV | 33 | 1.00 | 2.73 | 0.94 | 0.06 | s |
| -1 | O6 | INV | 33 | 1.00 | 2.85 | 0.73 | 0.09 | s |
| -1 | O6* | INV | 33 | 1.00 | 2.64 | 1.00 | 0.36 | s |
| +1 | O2 | VAR | 3 | 1.00 | 2.88 | 0.67 | 0.67 | F |
| +1 | O3 | INV | 3 | 1.00 | 2.63 | 0.67 | 0.33 | s |
| +1 | O6 | BUR | 3 | 1.00 | 3.82 | 0.00 | 0.33 | s |
| +1 | O6* | BUR | 3 | 1.00 | 2.84 | 0.67 | 1.00 | F |
| +2 | O2 | VAR | 47 | 1.00 | 2.79 | 0.74 | 0.45 | F |
| +2 | O3 | VAR | 47 | 1.00 | 2.84 | 0.87 | 0.23 | s |
| +2 | O6 | EXP | 47 | 1.00 | 4.11 | 0.02 | 0.98 | F |
| +2 | O6* | EXP | 47 | 1.00 | 3.11 | 0.32 | 1.00 | F |
| +3 | O2 | MC | 14 | 1.00 | 3.25 | 0.36 | 1.00 | F |
| +3 | O3 | MC | 14 | 1.00 | 3.48 | 0.14 | 1.00 | F |
| +3 | O6 | EXP | 14 | 1.00 | 4.58 | 0.00 | 1.00 | F |
| +3 | O6* | EXP | 14 | 1.00 | 3.39 | 0.43 | 1.00 | F |

OH-H over all units, all families (deposited O only): fits at 1.6/1.8/2.0 A in 1.000, 1.000, 1.000 of 1227 unit-oxygen sites.

Methyl fit (3.0 A) by 45 state of the cell (all families, all units, deposited O2/O3/O6):

| 45 state | n unit-sites | methyl fits | held units among them |
|---|---|---|---|
| BUR | 19 | 0.89 | 0.37 |
| EXP | 270 | 0.97 | 0.13 |
| INV | 516 | 0.19 | 0.98 |
| MC | 73 | 0.73 | 0.64 |
| VAR | 349 | 0.44 | 0.83 |

Methyl fit by the unit's own held flag (45, deposited O2/O3/O6): held=False: 0.94 of 344; held=True: 0.29 of 883.


## Check (b): sensitivity to the protein clash cut-off (2.8 / 3.0 / 3.2 A)

Number of family x subsite x oxygen cells (O2, O3, O6, O6*; all_glucoselike) in each category, and number of cells whose category differs from the 3.0 A result.

| substituent | F/s/N at 2.8 | F/s/N at 3.0 | F/s/N at 3.2 | F/s/N polar-aware | changed 2.8 vs 3.0 | changed 3.2 vs 3.0 | changed polar-aware vs 3.0 | mean unit fit fraction 2.8/3.0/3.2/polar-aware |
|---|---|---|---|---|---|---|---|---|
| methyl | 56/9/5 | 45/15/10 | 41/17/12 | 45/15/10 | 13 | 6 | 0 | 0.71/0.61/0.54/0.61 |
| hydroxyethyl | 41/12/17 | 38/13/19 | 35/12/23 | 38/13/19 | 4 | 6 | 0 | 0.55/0.50/0.44/0.50 |
| hydroxypropyl | 38/11/21 | 36/11/23 | 34/12/24 | 37/11/22 | 4 | 3 | 2 | 0.47/0.44/0.40/0.45 |
| acetyl | 33/13/24 | 26/16/28 | 24/16/30 | 30/15/25 | 9 | 4 | 5 | 0.38/0.33/0.29/0.36 |
| carboxymethyl | 36/13/21 | 33/12/25 | 31/13/26 | 36/12/22 | 6 | 2 | 5 | 0.45/0.40/0.37/0.44 |
| phosphate | 19/15/36 | 15/18/37 | 13/18/39 | 21/14/35 | 4 | 4 | 7 | 0.19/0.16/0.14/0.21 |
| succinyl | 29/16/25 | 24/16/30 | 20/18/32 | 28/17/25 | 8 | 6 | 7 | 0.36/0.31/0.27/0.35 |
| octenylsuccinyl | 26/14/30 | 20/17/33 | 18/18/34 | 24/18/28 | 7 | 2 | 7 | 0.32/0.28/0.25/0.32 |

## Literature comparison (substituted_starch_hydrolysis_literature.md, §7, §8, §11)

Our columns: all_glucoselike set, 3.0 A (category, groups); 2.8 and 3.2 A categories; 45 state; share of units in which the oxygen is held. Agreement rule (stated, mechanical): literature "tolerated" agrees with F; "not tolerated" agrees with N; s or a mixed literature verdict is "partial"; otherwise "no".

| cell | substituent | literature | our fit 3.0 A | 2.8 / 3.2 / polar-aware | glc units only (3.0) | 45 state | held | agree? (3.0) | agree? (polar-aware) |
|---|---|---|---|---|---|---|---|---|---|
| GH13 -1 O3 | phosphate | not productive (PPA, B. subtilis, A. oryzae; Takeda 1983/1986/1987) | N 0.00 (0/29; ent 0.00/70) | N / N / N | N 0.00 (0/18; ent 0.00/28) | INV | 0.99 | yes | yes |
| GH13 -1 O6* | phosphate | not productive (PPA, B. subtilis, A. oryzae; Takeda 1983/1986/1987) | N 0.00 (0/29; ent 0.00/70) | N / N / N | N 0.00 (0/18; ent 0.00/28) | INV | 0.96 | yes | yes |
| GH13 -1 O2 | hydroxyethyl | not tolerated at -1 (PPA; Yoshida 1984 II; position mixed O2/O3/O6) | s 0.00 (0/29; ent 0.01/70) | s / s / s | s 0.00 (0/18; ent 0.04/28) | INV | 1.00 | partial | partial |
| GH13 -1 O3 | hydroxyethyl | not tolerated at -1 (PPA; Yoshida 1984 II; position mixed) | s 0.00 (0/29; ent 0.01/70) | s / N / s | s 0.00 (0/18; ent 0.04/28) | INV | 0.99 | partial | partial |
| GH13 -1 O6* | hydroxyethyl | not tolerated at -1 (PPA; Yoshida 1984 II; position mixed) | s 0.28 (8/29; ent 0.29/70) | s / s / s | s 0.28 (5/18; ent 0.21/28) | INV | 0.96 | partial | partial |
| GH13 +1 O6* | phosphate | tolerated (PPA, A. oryzae); not productive (B. subtilis) | F 0.53 (9/17; ent 0.28/43) | F / s / F | s 0.33 (1/3; ent 0.33/3) | BUR | 0.33 | partial | partial |
| GH13 +1 O3 | phosphate | not productive (PPA, B. subtilis, A. oryzae) | N 0.00 (0/17; ent 0.00/43) | N / N / N | N 0.00 (0/3; ent 0.00/3) | INV | 1.00 | yes | yes |
| GH13 -2 O3 | phosphate | tolerated (PPA, B. subtilis); not productive (A. oryzae) | s 0.09 (2/22; ent 0.22/54) | s / s / s | s 0.10 (2/21; ent 0.23/53) | VAR | 0.75 | partial | partial |
| GH13 -2 O6* | phosphate | tolerated (B. subtilis); not productive (PPA, A. oryzae) | s 0.09 (2/22; ent 0.07/54) | s / s / s | s 0.10 (2/21; ent 0.07/53) | VAR | 0.70 | partial | partial |
| GH13 +2 O6* | phosphate | tolerated (B. subtilis) | F 0.67 (10/15; ent 0.83/41) | F / F / F | F 0.67 (10/15; ent 0.83/41) | EXP | 0.06 | yes | yes |
| GH13 +2 O3 | phosphate | tolerated (A. oryzae); not productive (PPA, B. subtilis) | N 0.00 (0/15; ent 0.00/41) | N / N / s | N 0.00 (0/15; ent 0.00/41) | VAR | 0.92 | partial | partial |
| GH13 -3 O6* | phosphate | tolerated (PPA, A. oryzae) | F 0.53 (9/17; ent 0.34/38) | F / F / F | F 0.50 (4/8; ent 0.35/20) | EXP | 0.20 | yes | yes |
| GH13 -3 O3 | phosphate | tolerated (A. oryzae; liquefying amylase, Kamasaka 1997) | s 0.47 (8/17; ent 0.42/38) | s / s / F | F 0.62 (5/8; ent 0.50/20) | EXP | 0.12 | partial | yes |
| GH13 +3 O3 | phosphate | tolerated (PPA, B. subtilis) | s 0.20 (2/10; ent 0.15/13) | s / s / s | s 0.20 (2/10; ent 0.15/13) | MC | 0.43 | partial | partial |
| GH13 +1 O2 | methyl | terminal residues of B. subtilis limit dextrins >90 % unsubstituted (van der Burgt 2000b; position mixed) | s 0.35 (6/17; ent 0.28/43) | F / s / s | F 0.67 (2/3; ent 0.67/3) | VAR | 0.86 | partial | partial |
| GH13 +1 O3 | methyl | idem (position mixed) | s 0.06 (1/17; ent 0.09/43) | s / s / s | s 0.33 (1/3; ent 0.33/3) | INV | 1.00 | partial | partial |
| GH13 +1 O6* | methyl | idem (position mixed) | F 1.00 (17/17; ent 0.98/43) | F / F / F | F 1.00 (3/3; ent 1.00/3) | BUR | 0.33 | no | no |
| GH14 -1 O6* | phosphate | not productive (sweet-potato beta-amylase; Takeda 1987) | N 0.00 (0/2; ent 0.00/16) | N / N / N | N 0.00 (0/2; ent 0.00/13) | INV | 1.00 | yes | yes |
| GH14 -2 O6* | phosphate | not productive (sweet-potato beta-amylase; Takeda 1987) | N 0.00 (0/2; ent 0.00/15) | N / N / N | N 0.00 (0/2; ent 0.00/15) | INV | 1.00 | yes | yes |
| GH14 +1 O6* | phosphate | tolerated (sweet-potato beta-amylase; Takeda 1987) | N 0.00 (0/2; ent 0.00/3) | N / N / N | N 0.00 (0/2; ent 0.00/3) | INV | 0.67 | no | no |
| GH14 +1 O3 | phosphate | not productive (Takeda 1987) | N 0.00 (0/2; ent 0.00/3) | N / N / N | N 0.00 (0/2; ent 0.00/3) | INV | 1.00 | yes | yes |
| GH14 -1 O3 | phosphate | not productive (Takeda 1987) | s 0.00 (0/2; ent 0.07/15) | s / s / s | N 0.00 (0/2; ent 0.00/13) | MC | 0.33 | partial | partial |
| GH14 -2 O3 | phosphate | not productive (Takeda 1987) | N 0.00 (0/2; ent 0.00/15) | N / N / N | N 0.00 (0/2; ent 0.00/15) | INV | 1.00 | yes | yes |
| GH14 -2 O2 | methyl | not tolerated at -2/-1 (B. cereus beta-amylase; van der Burgt 2000b; position mixed) | s 0.00 (0/2; ent 0.13/15) | s / s / s | s 0.00 (0/2; ent 0.13/15) | VAR | 0.89 | partial | partial |
| GH14 -1 O2 | methyl | idem (position mixed) | F 1.00 (2/2; ent 0.94/16) | F / F / F | F 1.00 (2/2; ent 0.92/13) | MC | 1.00 | no | no |
| GH14 +1 O2 | methyl | not tolerated at +1 (van der Burgt 2000b; position mixed) | N 0.00 (0/2; ent 0.00/3) | N / N / N | N 0.00 (0/2; ent 0.00/3) | MC | 1.00 | yes | yes |
| GH15 -1 O6* | phosphate | not productive (A. niger GA; Abe 1982, Takeda 1983/1987) | N 0.00 (0/2; ent 0.00/4) | N / N / N | — | INV | 1.00 | yes | yes |
| GH15 -1 O3 | phosphate | not productive (A. niger GA) | N 0.00 (0/2; ent 0.00/4) | N / N / N | — | INV | 1.00 | yes | yes |
| GH15 +1 O6* | phosphate | tolerated (A. niger GA; Takeda 1983/1987; Kamasaka 1997) | N 0.00 (0/2; ent 0.00/4) | N / N / N | — | — | — | no | no |
| GH15 +1 O3 | phosphate | not productive (A. niger GA) | N 0.00 (0/2; ent 0.00/4) | N / N / N | — | INV | 1.00 | yes | yes |
| GH15 +2 O3 | phosphate | tolerated (implied: one Glc left on non-reducing side) | F 0.50 (1/2; ent 0.50/4) | F / F / F | F 0.50 (1/2; ent 0.50/4) | EXP | 0.00 | yes | yes |
| GH15 -1 O2 | methyl | O-methyl blocks (Lemieux 1996, methyl isomaltoside) | N 0.00 (0/2; ent 0.00/4) | N / N / N | — | INV | 1.00 | yes | yes |
| GH15 -1 O3 | methyl | O-methyl blocks (Lemieux 1996) | s 0.00 (0/2; ent 0.25/4) | F / N / s | — | INV | 1.00 | partial | partial |
| GH15 -1 O6* | methyl | O-methyl blocks (Lemieux 1996) | F 1.00 (2/2; ent 0.75/4) | F / F / F | — | INV | 1.00 | no | no |
| GH15 +1 O2 | methyl | O-methyl tolerated, kinetics unaffected (Lemieux 1996; isomaltoside +1) | F 1.00 (2/2; ent 1.00/4) | F / F / F | — | INV | 1.00 | yes | yes |
| GH15 +1 O3 | methyl | O-methyl costs 1.9 kcal/mol (~25-fold kcat/KM) (Lemieux 1996) | N 0.00 (0/2; ent 0.00/4) | N / N / N | — | INV | 1.00 | partial | partial |
