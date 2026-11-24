# 09_mmgbsa — MM/GBSA (calculation 9)

Started 2026-11-23. Plan: `docs/03_supporting_calculations.md`. Database opened read-only. Runs on GUANE <cluster-account>.

## Status

- MD-based part: pipeline tested on the 1 ns test of S01_6Z8L (`md_test_S01_1ns/`); waits for the 100 ns runs of
  calculation 8.
- Series part: GAA iminosugars, jobs 33442 (first trial, 5NN5) and 33445 (all four, final protocol) on <cluster-account>.
- ChEMBL tables appeared in the database at 15:31 (`chembl_wait.log`); `resultados/calc/07_docking_chembl/` did not
  exist at any check, so no docking poses are used.

## Tools and versions

AmberTools on GUANE <cluster-account>, env `~/enz_md/soft/amber` (micromamba, conda-forge): installed as AmberTools 24
(used for the tleap builds of calculation 8), upgraded to AmberTools 26 by the solver when RDKit was added; MMPBSA.py,
sander, cpptraj, antechamber of AmberTools 26 run the MM/GBSA. AmberTools 26 has no `reduce`; for the series the
reduce output (AmberTools 24, laptop) is supplied as `series/<PDB>_<COMP>/reduce_out.pdb`.
gmx_MMPBSA was not used: MMPBSA.py reads the Amber topologies the systems were built with.

## MD-based part (`scripts/calc/c09_mmgbsa_md.sh`)

Input: the solute trajectory of calculation 8 (protein + glucan + crystal ions, every 10 ps) and the topologies of
`c08_solute.py` (complex = solute, receptor = solute without the glucan, ligand = glucan; mbondi2 radii).
1. cpptraj: frames 11, 21, … (every 100 ps from 0.1 ns, 1000 frames for 100 ns), autoimage.
2. MMPBSA.py, single trajectory: GB igb 5 (OBC II), saltcon 0.150 M, no cut-off, SASA non-polar term (gbsa 2);
   ΔG on all 1000 frames; per-residue decomposition (idecomp 1, 1-4 terms in "internal") on every 10th frame
   (every 1 ns), residues within 6 Å of the glucan in the first frame. No entropy term.

Test (1 ns of S01_6Z8L, 10 frames): ΔG = −22.6 ± 3.1 kcal/mol (SD; SEM 1.0); vdW −21.1, EEL −42.2, EGB +43.8,
ESURF −3.1. Residues with the largest share: Trp284 −4.6, Tyr276 −3.9, Glu272 −3.2 kcal/mol. 14.3 min on one legacy core for
10 frames with decomposition, which set the split above.

## Series part (`scripts/calc/c09_series.py`, `c09_series_summary.py`)

Congeneric series = the inhibitors that have both a deposited complex in the database and a pChEMBL against the same
enzyme: human lysosomal α-glucosidase GAA (UniProt P10253, ChEMBL CHEMBL2608), four iminosugars:

| PDB | Inhibitor | pChEMBL used (median) | Activities used |
|---|---|---|---|
| 5NN5 | 1-deoxynojirimycin (NOJ) | 6.82 | Ki 6.82, IC50 6.00, IC50 5.82, Ki 7.23, Ki 7.23 |
| 5NN6 | miglitol (MIG) | 6.46 | IC50 6.46 |
| 9GTC | 2-pentyl-DNJ analogue (A1IO7) | 7.72 | Ki 7.72 |
| 9GTW | 2-heptyl pyrrolidine (A1IPC) | 5.68 | IC50 5.68 |

Activities kept: target CHEMBL2608, assay description naming the human lysosomal enzyme (or GAA), cell assays
(Caco-2) dropped. Other ChEMBL molecules with structures act on enzymes with no second congeneric partner.
Method: chain A, inhibitor with its ring N protonated (+1), GAFF2/AM1-BCC; ff14SB; GB minimisation with protein
heavy atoms > 6 Å from the inhibitor restrained; single-structure MMPBSA.py (igb 5, 0.15 M) with decomposition.
No MD for the series (no free node for four ~13 000-atom complexes; docking poses absent).

## Failures and their cause

- `print_res="within 6"` is not accepted by MMPBSA.py (gmx_MMPBSA syntax): residue lists computed with cpptraj.
- Series first trial (job 33442): GB minimisation without cut-off, 2000 steps, ~9 s per step → replaced by
  cut-off 16 Å, 1000 steps (job 33445, waits on 33442).
