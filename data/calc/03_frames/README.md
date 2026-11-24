# 03_frames — family frames, conserved waters and ions, ligand centres in the frame

Calculation 3 of `docs/03_supporting_calculations.md`. Run on the laptop, 2026-11-23.

## What was done

For each family with ≥ 2 entries (10 families; GH4, GH126, GH133, GH176 have one entry and are not done):

1. **Reference structure.** Part 3 reference where one exists (first row of `part3_site_geometry` for the family), with its two
   catalytic residues (nucleophile/base, acid) in author numbering of that structure. GH57: 1K1X (Part 3 family reference;
   not one of the 637 entries; mmCIF downloaded to `ref_cif/1K1X.cif`). Families without a Part 3 reference (AA13, GH49, GH63,
   GH97, GT35): the best-resolution entry holding a sugar/'other' ligand; its "site residues" are the protein residues within
   4.0 Å of those ligands (database table `contact`), on the chain with most of them.
2. **Core.** Every protein chain (≥ 30 amino acids with Cα) of every entry is aligned by sequence to the reference chain
   (Biopython 1.88 `PairwiseAligner`, global, BLOSUM62, gap open −10, extend −0.5, end gaps free). Core = aligned pairs
   whose reference Cα is within 12.0 Å of any atom of the reference catalytic (site) residues.
3. **Superposition.** Kabsch on the core Cα; two starting fits (core pairs; all aligned pairs trimmed at 3 Å), each refined
   by refitting on the pairs within 2.0 Å until stable; the fit with more retained pairs is kept. Accepted when ≥ 20 pairs
   and ≥ 40 % of the core pairs lie within 2.0 Å. **Fallback** when the sequence alignment fails this test: structural
   matching of Cα in the reference core (mutual nearest Cα within 3.0 Å, iterated), started from the same two fits; the
   method is stated in the `status` column.
4. **Into the frame.** Waters (HOH, from the mmCIF), ions (database class `ion`) and the centres of all sugar/'other'
   instances (database `instance.x,y,z`) are given to the nearest accepted chain (nearest protein heavy atom) and moved with
   that chain's rotation and translation.
5. **Clusters.** Waters and ions separately: greedy density clustering — seed = the point whose 1.0 Å neighbourhood holds
   the most distinct entries; cluster = unassigned points within 1.0 Å of the seed; repeated while the seed reaches 50 % of
   the family's superposed entries. Clusters present in ≥ 50 % of the superposed entries are reported.
6. **Catalytic residues of every chain.** The reference catalytic residues carried to each accepted chain by the alignment;
   if that residue is not of the reference type or its functional atoms (carboxylate O etc.) lie > 3.0 Å from the reference
   ones after superposition, the nearest residue of the reference type within 3.0 Å is taken (column `method`).

Tools: Python 3.13 (`~/miniconda3/bin/python`), gemmi 0.7.5, numpy, scipy 1.16.2, Biopython 1.88.
UniProt numbers from the database table `sifts_residue` (SIFTS; 619 of 645 entries mapped at the time
of this run). **Correction (2026-11-23, same day):** the 8 GH63 entries identified below as crescentin
filaments (not the YgjK α-glucosidase) were all SIFTS-mapped and were removed from `sifts_residue` when
excluded from the database; the current, corrected figure is 611 of 637 entries mapped
(`resultados/db/sifts_summary.json`).

Command: `~/miniconda3/bin/python scripts/calc/c03_frames.py` (all families; one process, about 7 min).

## Files

| File | One row is |
|---|---|
| `superposition.tsv` | one chain of one entry: status and method, sequence identity over aligned pairs, core pairs, retained pairs, RMSD (first fit, all core pairs, retained pairs), rotation R11..R33 and translation t1..t3 (x_frame = R x + t) |
| `catalytic_residues_mapped.tsv` | one reference catalytic (or site) residue in one chain: residue found, author number, UniProt number, type check, functional-atom distance to the reference after superposition, method |
| `conserved_clusters.tsv` | one water or ion cluster in ≥ 50 % of the superposed entries: entries, members, rms and max radius, centre in the frame, distance to each catalytic residue (functional atoms), nearest reference residue (with UniProt number) |
| `conserved_cluster_members.tsv` | one water or ion in a reported cluster: PDB id, chain, residue number, assigned protein chain, nearest protein residue (UniProt number), frame coordinates |
| `ligand_centres_in_frame.tsv` | one sugar/'other' instance: PDB id, instance id, chain, residue number, component, RSCC, assigned chain, nearest residue (UniProt number), frame coordinates, distance to each catalytic residue, nearest reference residue |
| `families.json` | per family: reference, rule, residues defining the core, counts |
| `ref_cif/1K1X.cif` | GH57 reference coordinates (RCSB) |
| `run.log` | output of the run |

## Counts

| Family | Reference (chain) | Residues defining the core | Core Cα | Entries | Entries superposed | Chains accepted (sequence / structural fallback) | Chains rejected | RMSD retained, median (max) Å | Waters | Ions | Ligand centres | Water clusters ≥ 50 % | Ion clusters ≥ 50 % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AA13 | 5LSV (A) | 13 ligand-contact residues | 115 | 8 | 8 | 9 / 0 | 0 | 0.48 (0.57) | 2495 | 61 | 6 | 212 | 4 |
| GH13 | 1AMY (A) | Asp179 nuc., Glu204 acid | 76 | 412 | 411 | 530 / 69 | 18 | 0.79 (1.52) | 191540 | 1077 | 1705 | 5 | 1 |
| GH14 | 1BYA (A) | Glu380 base, Glu186 acid | 99 | 50 | 50 | 74 / 3 | 0 | 0.43 (1.36) | 24046 | 38 | 216 | 144 | 0 |
| GH15 | 1AGM (A) | Glu400 base, Glu179 acid | 95 | 18 | 18 | 20 / 4 | 0 | 0.47 (1.50) | 10033 | 3 | 125 | 151 | 0 |
| GH31 | 2QMJ (A) | Asp443 nuc., Asp542 acid | 91 | 73 | 73 | 149 / 0 | 5 | 0.58 (1.12) | 50144 | 124 | 219 | 55 | 0 |
| GH49 | 2Z8G (B) | 17 ligand-contact residues | 149 | 4 | 4 | 10 / 0 | 0 | 0.26 (0.40) | 3656 | 0 | 10 | 664 | 0 |
| GH57 | 1K1X (A) | Glu123 nuc., Asp214 acid | 95 | 25 | 25 | 67 / 0 | 0 | 1.07 (1.12) | 12580 | 0 | 295 | 198 | 0 |
| GH63 | 3W7T (B) | 29 ligand-contact residues | 402 | 17 | 9 | 17 / 0 | 120 | 0.38 (0.90) | 9974 | 31 | 40 | 533 | 1 |
| GH97 | 9BS5 (B) | 20 ligand-contact residues | 228 | 7 | 7 | 14 / 0 | 0 | 0.24 (0.27) | 10584 | 16 | 16 | 807 | 1 |
| GT35 | 1L5W (A) | 20 ligand-contact residues | 215 | 27 | 27 | 56 / 6 | 0 | 0.80 (1.43) | 13496 | 0 | 156 | 32 | 0 |

- 2 788 sugar/'other' instance centres placed in a family frame (2 623 sugar, 165 other) of the 2 818 in the database; the
  rest are in single-entry families or in entries not superposed.
- Catalytic residues (Part 3 families): 1 832 residue-in-chain rows; 1 566 confirmed (reference type, functional atoms
  within 3.0 Å of the reference after superposition), 177 of them by the spatial step.
- GH13 conserved sites: one calcium site in 51.6 % of the 411 entries (245 Ca, 14 Mn, 8 Na, 4 Sr, 1 K; nearest
  reference residue Gly183, 10.0 Å from Asp179); five water sites, two of them at the catalytic residues (W3, 3.1 Å from
  the Asp179 carboxylate, 61 % of entries; W4, 2.4 Å from the Glu204 carboxylate, 57 %).
- The number of water clusters depends on how uniform the family is: GH49 (4 entries of one protein) gives 664, GH13
  (many proteins) gives 5. The 50 % threshold of small families is 2–5 entries.

## Failures and their causes

- **GH63: 8 of 17 entries not superposed** (8AFE, 8AFH, 8AFL, 8AFM, 8AHL, 8AIA, 8AIX, 8AJB; 120 chains). These are cryo-EM
  structures of crescentin filaments (`_struct.title` of 8AFE: "Cryo-EM structure of crescentin filaments …"), not the
  GH63 α-glucosidase YgjK under which the database lists them (entry.family = GH63, UniProt P42592). Sequence identity to
  YgjK ≈ 0.28 over aligned pairs, no core fit. This is an error of the entry list (Part 2), not of the superposition.
- **GH13: 4OKD not superposed** (Chlamydomonas isoamylase ISA1 + maltoheptaose): neither the sequence alignment nor the
  structural matching retains 40 % of the core pairs within 2 Å of barley α-amylase (19/38, 17/35).
- **Other rejected chains.** Partner proteins (entity description in the mmCIF): GH13 1BVN T (tendamistat), 1CLV I
  (α-amylase inhibitor), 1TMQ B (ragi bifunctional inhibitor), 1KXQ E–H (VHH antibody fragment), 4X0N B (helianthamide),
  4CVW C–D (limit dextrinase inhibitor). Not checked one by one: 1KXV C–D, 9XVP/9XVV B (in the ISA2–ISA1–ISA1
  heterotrimer entries), 1BPL A (entity 1 of a two-entity α-1,4-glucan-4-glucanohydrolase entry; 19/36 pairs), 9US6 B
  (second copy of the maltohexaose-producing amylase; chain A accepted; 19/42 pairs). GH31: 5 chains with 0–18 core pairs
  (7P2Z AAA, 8CB1 G, 8CB6 F, 8D43 B, 8EMR B).
- Structural fallback chains (82) have fewer retained pairs; GH15 1LF6/1LF9 (bacterial glucoamylase) keep 27–33 of
  46–54 core pairs (RMSD 1.3–1.5 Å) and their catalytic residues are not found (gap or wrong type): their frame is coarse.
- Single linkage at 1.0 Å (the first clustering tried) chained the GH13 waters into one cluster of 139 214 waters with a
  45.6 Å radius (GH14: 3.9 Å, GH31: 5.8 Å); replaced by the greedy density clustering above (all clusters now ≤ 1.5 Å max
  radius).
- Insertion codes: residues are matched by number and insertion code throughout; SIFTS-unmapped entries (26) have empty
  UniProt columns.
