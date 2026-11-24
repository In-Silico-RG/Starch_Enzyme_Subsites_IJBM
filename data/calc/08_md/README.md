# 08_md — short MD of representative enzyme–glucan complexes (calculation 8)

Started 2026-11-23. Plan and rules: `docs/03_supporting_calculations.md`. Database opened read-only.

## Status

- (a) selection: done, 12 systems (`selection/selection.tsv`, `selection/candidates.tsv`).
- (b) parametrisation: done, all 12 built on GUANE (<cluster-account>, legacy node), energy check passed.
- (c) one system end to end: done, S01_6Z8L, minimisation to 1 ns (job 33424, <cluster-account>), analysed.
- (d) production 100 ns × 12: RUNNING since 2026-11-23 15:44 -05 (jobs 33435, 33437, 33440; see table).
- (e) analysis: script ready and tested on the 1 ns run; to run on each 100 ns trajectory.

## (a) Selection

Script `scripts/calc/c08_select.py` → `selection/candidates.tsv` (161 oligosaccharides, one row each).
A candidate passes when the branched oligosaccharide (one `label_asym` of kind `branched` in table `instance`) is
glucose only (GLC/BGC), 2–7 residues, every residue RSCC ≥ 0.80, the entry is X-ray at ≤ 2.2 Å, the title,
entity description, `_entity.pdbx_mutation` and `_struct_ref_seq_dif` (engineered mutation) show no mutation,
it is not cyclic (glycosidic bonds = residues) and has no covalent bond to anything outside it.

| Filter | Oligosaccharides removed |
|---|---|
| mutant (title/description/mutation records) | 70 |
| cyclic (cyclodextrins) | 15 |
| covalent bond to another molecule (e.g. 9FZ2) | 4 |
| pass all | 83 (42 entries; AA13 2, GH13 29, GH14 5, GH49 1, GH57 1, GT35 4) |

No wild-type glucose-only glucan passes in GH15, GH31, GH63, GH97, GH4, GH126, GH133, GH176: those families are
not represented. The mutant flag is conservative: 1WO2 is flagged only because of a pyroglutamate (Q16PCA)
record; it is the same enzyme type as S01 and was not needed.

Chosen (one per family that has a candidate, plus GH13 enzyme types; CAZy subfamily numbers are not in the
database, so GH13 types are by EC number and protein):

| System | PDB | Family | Enzyme | Res. (Å) | Glucan | Why |
|---|---|---|---|---|---|---|
| S01_6Z8L | 6Z8L | GH13 | human pancreatic α-amylase | 1.40 | maltose, linkage anti-ψ (φ 90.8, ψ 73.7) | the required anti-ψ entry; 1MFU, 1MFV, 3WEO fail "glucose only" (AGL or AC1 in the same oligosaccharide) |
| S02_2FHF | 2FHF | GH13 | K. pneumoniae pullulanase | 1.65 | maltotetraose | debranching enzyme; barley limit dextrinase entries have chain gaps |
| S03_6TP1 | 6TP1 | GH13 | B. licheniformis α-amylase | 1.94 | maltotetraose | bacterial liquefying amylase, Ca–Na–Ca site |
| S04_2GVY | 2GVY | GH13 | A. niger neutral α-amylase | 1.80 | maltose | fungal (Taka-type) amylase |
| S05_1QHO | 1QHO | GH13 | G. stearothermophilus maltogenic amylase | 1.70 | maltose | five-domain CGTase-like enzyme |
| S06_7DCG | 7DCG | GH13 | W. cibaria α-glucosidase | 1.53 | maltotriose | α-glucosidase type |
| S07_6A0K | 6A0K | GH13 | A. globiformis cyclic maltosyl-maltose hydrolase | 1.94 | trisaccharide α(1→4) + α(1→6) | only wild-type GH13 glucan with an α(1→6) linkage |
| S08_1J18 | 1J18 | GH14 | B. cereus β-amylase | 2.00 | maltotetraose (one linkage φ 175.9, ψ 72.6) | longest wild-type GH14 glucan |
| S09_9KYQ | 9KYQ | GH57 | A. aeolicus amylopullulanase | 1.60 | maltotriose | only non-cyclic wild-type GH57 glucan |
| S10_2Z8G | 2Z8G | GH49 | A. brasiliensis isopullulanase | 1.70 | isopanose | only wild-type GH49 glucan |
| S11_5T7N | 5T7N | AA13 | A. oryzae AA13 LPMO | 1.60 | branched tetrasaccharide | longer of the two AA13 glucans |
| S12_1L5W | 1L5W | GT35 | E. coli maltodextrin phosphorylase | 1.80 | maltotetraose | best GT35 glucan; needs PLP and phosphate parameters |

## (b) Parametrisation

Force fields: ff14SB (protein) + GLYCAM06j-1 (glucan) + TIP3P, Joung–Cheatham Na⁺/Cl⁻, Li–Merz 12-6 Ca²⁺/Zn²⁺;
the combination of the project's GLYCAM work (Proyecto_C). 0.15 M NaCl (n = round(0.15 × N_water / 55.5))
after neutralisation; truncated octahedron, 10 Å of water around the solute. Deposited Ca²⁺ (and Na⁺, Cl⁻, Zn²⁺)
within 3.5 Å of the kept chain are kept; crystal waters within 3.5 Å kept unless closer than 2.2 Å to a kept heavy
atom or kept water. Other ligands, other sugars, buffers and N-glycans are removed.

Scripts: `scripts/calc/c08_build.py` (structure, protonation, GLYCAM names, tleap), `c08_gmx.py` (GROMACS
topology), `c08_solute.py` (solute and MM/GBSA topologies), `c08_check_energy.py` (energy check),
`c08_gt35_params.py` (GT35 cofactors). Tools: AmberTools 24 tleap/ParmEd 4.3 (laptop env `sc_modeling`, and
AmberTools 24 on GUANE <cluster-account> for the final builds), PROPKA 3.5.1, reduce (AmberTools 24), RDKit, antechamber/sqm.

Rules in the build (each stated once):
- protonation at pH 7 by PROPKA 3 on protein + glucan + ions: Asp/Glu pKa > 7 → ASH/GLH, except carboxylates with
  an O within 3.0 Å of a metal ion (PROPKA ignores the ions); His pKa > 7 → HIP; neutral His tautomer from reduce;
  His N within 2.6 Å of a metal left unprotonated; disulfides SG–SG < 2.5 Å.
- residues without ff14SB parameters: pyroglutamate PCA (6Z8L, N-terminal) → ACE cap from its CA, C, O;
  N-methyl-His HIC (5T7N, N-terminal) → His without the methyl.
- glucan → GLYCAM residues (0GA/4GA/6GA/UGA, B for β), reducing-end O1 in an ROH residue, deposited glycosidic
  bonds added explicitly.
- GT35 (S12): Lys645–PLP Schiff base as residue LLP (net −2: phosphate 2−, imine N protonated, O3′ phenolate);
  backbone + CB/CG/CD keep ff14SB LYS types and charges, the rest GAFF2 with AM1-BCC charges shifted by −0.00525 e
  per atom to total −2; the 6 bonds/angles/torsions across the ff14SB/GAFF2 boundary taken from GAFF2 analogues
  (`params_gt35/frcmod.llp_boundary`); phosphate 998 (2.6 Å from the glucan and from the PLP phosphate) as HPO4²⁻,
  GAFF2/AM1-BCC. Parameters in `params_gt35/`.
- 1-4 scaling: ff14SB 1/1.2 and 1/2, GLYCAM none. GROMACS takes one pair of factors, so the glucan (its own
  moleculetype) has every 1-4 pair written as a type-2 pair with factor 1.0 and unscaled LJ (`c08_gmx.py`).

Energy check (`systems/*/energy_check.txt`): single point of system.gro, OpenMM 8.2 on the prmtop against GROMACS
2024.4 on the converted top (same PME/LJ settings). Bonds, angles, dihedrals and the 1-4 Coulomb and LJ sums agree
within 0.4 kJ/mol in all 12 systems (e.g. S01: 1-4 Coulomb 85 300.36 vs 85 300.40; S12: 147 006.07 vs 147 006.00).
A wrong glucan 1-4 factor would show as about −416 kJ/mol in S01 (126 pairs; 1-4 Coulomb 2 053.5, LJ 147.3 kJ/mol).
The totals differ by 4–55 kJ/mol (S02_2FHF: 438 of 1.6 × 10⁶), from the PME and dispersion-correction
implementations, not from the topology.

| System | Atoms | Solute atoms | Protein residues | Crystal ions kept | Crystal waters | Solute charge |
|---|---|---|---|---|---|---|
| S01_6Z8L | 56 533 | 7 708 | 496 | Ca906, Cl907 | 340 | +1 |
| S02_2FHF | 144 128 | 15 903 | 1 052 | Ca2401, 2402, 2404, 2405 (Ca2403 > 3.5 Å from chain A) | 938 | −37 |
| S03_6TP1 | 64 325 | 7 653 | 481 | Ca501, Ca502, Na503, Na504 | 344 | −1 |
| S04_2GVY | 61 434 | 7 215 | 476 | Ca801 | 316 | −21 |
| S05_1QHO | 59 381 | 10 432 | 686 | Ca701, 702, 703 | 839 | −1 |
| S06_7DCG | 77 541 | 9 475 | 588 | Ca1306 | 482 | −26 |
| S07_6A0K | 69 369 | 6 845 | 450 | none deposited | 187 | −27 |
| S08_1J18 | 69 745 | 8 194 | 516 | Ca1500 | 351 | +1 |
| S09_9KYQ | 47 731 | 8 092 | 476 | none deposited | 338 | −14 |
| S10_2Z8G | 60 438 | 8 310 | 549 | none deposited | 540 | −10 |
| S11_5T7N | 28 215 | 3 510 | 233 | Zn406–409 | 264 | −10 |
| S12_1L5W | 76 639 | 12 840 | 796 (+ LLP, HP4) | none deposited | 489 | −6.001 (AM1-BCC rounding of HP4, 0.001 e) |

## (c) Test of one system: S01_6Z8L, job 33424 (<cluster-account>, exadell, 64 threads)

em1 (heavy atoms restrained 1000 kJ mol⁻¹ nm⁻²; converged Fmax < 500 in 3 252 steps), em2 (7 steps), NVT 100 ps
(restrained), NPT 500 ps (1000), NPT 500 ps (100), production 1 ns unrestrained, 300 K V-rescale, 1 bar C-rescale,
2 fs, h-bond constraints, PME 1.0 nm, LJ 1.0 nm potential-shift + dispersion correction. Speed 92.6 ns/day.
Output 2 970 684 bytes for 101 solute frames (7 708 atoms) → 3.8 bytes/atom/frame, used to budget the quota.
Multi-rank benchmarks failed (one task per job; mpirun had one slot) and were dropped.
1 ns analysis (`analysis/S01_6Z8L/` from `systems/S01_6Z8L/test_run/test.xtc`, to be overwritten by the 100 ns
analysis): glucan RMSD 2.1 Å mean; anti-ψ linkage stayed anti-ψ in all 101 frames (φ 86 ± 9°, ψ 64 ± 14°);
5 of 6 deposited contacts kept ≥ 50 %; Ca²⁺ kept its 5 protein ligands in every frame. The maltose of 6Z8L sits at
a surface site (Glu272, Tyr276, Trp284), not in the catalytic cleft. Its reducing ring is deposited distorted
(θ = 82.7°, census value) and moved to θ ≈ 130–175° within 1 ns; the anomeric configuration was checked (α).

## (d) Production

Partition `partition.tsv` written and checked (disjoint, complete) by `scripts/calc/c08_partition.py`; job scripts
`scripts/calc/c08_md/job_*.sh` (from `make_jobs.py`), runner `run_system.sh` / `run_system_gpu.sh`.
100 ns, 1 replica, 300 K; solute (protein + glucan + crystal ions) written every 10 ps (`compressed-x-grps =
Protein_Lig`). Waters are not written: GROMACS writes one fixed atom group, so "nearby water" cannot be a
group; the full final frame (md.gro) is kept.

| Job | Account | Node | Systems (threads) |
|---|---|---|---|
| 33435 | <cluster-account> | smexa (amd, 124 cores) | S12 (28), S08 (25), S03 (23), S05 (21), S09 (17), S11 (10) |
| 33437 | <cluster-account> | exadell (amd, 126 cores) | S06 (30), S07 (27), S04 (24), S10 (23), S01 (22) |
| 33440 | <cluster-account> | felix (Titan X, 24 cores) | S02 |

<cluster-account> was not used: 3.0 GB free in its HOME at submission. Failures: 33436 and 33438 (felix) stopped at
em1, GROMACS 2024.4 cannot put PME or bonded terms on the GPU with steepest descent; `run_system_gpu.sh` puts them
on the CPU for em1/em2. The amd jobs first waited on memory (the other user's jobs hold 370 of 385 GB per node);
their request was lowered to 14 000 MB with `scontrol update` while pending.
NVT speeds: 21–30 ns/day per system on amd, 43 ns/day for S02 on the GPU → about 4 days for 100 ns.

## (e) Analysis

`scripts/calc/c08_analyse.py SYSTEM` → `analysis/<SYSTEM>/`: ligand RMSD (protein CA superposed on the deposited
CA), persistence of each deposited glucan–residue contact (≤ 4.0 Å, deposited pairs from the asymmetric unit),
φ/ψ(/ω) of every glycosidic linkage with the census anti-ψ criterion |ψ + 117.8| > 90, Cremer–Pople puckering,
ion coordination (protein + glucan O/N within 3.0 Å; water is not in the trajectory) and ion displacement.

## Failures and their cause

- First laptop build of S01: ParmEd refused the GROMACS export (mixed 1-4 scaling) → `c08_gmx.py`.
- Chloride of 6Z8L lost at a 3.2 Å ion cut-off → 3.5 Å.
- 5T7N first build: total energy 1.6 × 10⁹ kJ/mol from crystal waters overlapping the glucan and His90 → waters
  closer than 2.2 Å dropped; all 12 rebuilt with the rule.
- Builds moved from the laptop to GUANE at AFC's order (laptop kept free) on 2026-11-23.
