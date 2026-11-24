# Figure 1 candidates (protein overview; replaces the methyl β-maltoside figure) — 2026-11-24

Files: `fig1_protein_A.{png,pdf}` (7.00 × 8.40 in), `fig1_protein_B.{png,pdf}` (7.06 × 8.34 in), 300 dpi.
Style: matched to the approved Figure 2 (`scripts/figures/fig2_v3_*.py`): ChimeraX, soft lighting and black
silhouettes on white; faint grey protein; sugar gold #F4B942 (O #E03C31, N #2E6DB4) with the subsite numbers
written on the rings (8 pt bold, white halo); catalytic residues magenta #B0308F; Liberation Sans, 7 pt text,
7.5 pt titles, 9 pt bold panel letters. The scissors glyph uses DejaVu Sans because Liberation Sans lacks it.

Scripts (all in `scripts/figures/`): `fig1_protein_common.py` (paths, palette); `fig1_protein_frame.py`
(ring centroids in the family frame); `fig1_protein_struct.py` (entries, domains, catalytic residues, the
view search with PyMOL montages); `fig1_protein_cx.py` (ChimeraX renders); `fig1_protein_chem.py` (RDKit
depictions); `fig1_protein_layouts.py` (assembly). Run order: frame → cx → layouts. Working files: `work/`.

## What each candidate shows

**A (proteins first).** a, GH13 6JQB; b, GH31 3WEO: cartoon with a faint tint per domain under a 60 %
transparent grey surface. c, GH13 frame; d, GH31 frame. e, chemistry (maltohexaose over −3…+3 with the
scissile bond; glucose unit with R at O2, O3, O6 and the four substituents in one row).

**B (question first).** a, chemistry (substituents as a block at the right). b, 6JQB as an opaque surface
(the cleft reads as a groove holding the chain); c, GH13 frame. d, 3WEO opaque surface (−1 sits in a pocket);
e, GH31 frame.

The catalytic close-up first planned is left out: Figure 2 already shows both active sites.

## Draft caption (A lettering)

**Figure 1.** Starch-degrading glycoside hydrolases, their subsites and the substituents of functionalised
starch.
(a) GH13 maltotetraose-forming amylase of *Pseudomonas saccharophila* STB07 (PDB 6JQB, 1.10 Å). An
acarbose-derived heptasaccharide is bound from subsite −4 to +3. Domains are A ((β/α)8 barrel), B and C. Units at
−3 to +3 are gold and numbered; the −4 unit is grey. The catalytic Asp193, Glu219 and Asp294 are magenta.
(b) GH31 sugar beet α-glucosidase (PDB 3WEO, 1.45 Å). Acarviosyl-maltohexaose is bound from −1 to +7. Domains
are the N-terminal β-sandwich (N), the (β/α)8 barrel (A), and the proximal (C) and distal (C′) C-terminal
domains. Units −1 to +3 are gold; +4 to +7 are grey. The catalytic Asp469 and Asp568 are magenta.
(c, d) Ring centroids of every density-supported sugar unit (RSCC ≥ 0.80) at subsites −3 to +3. The units come
from 79 GH13 entries (308 units) and 18 GH31 entries (82 units). Each is superposed on its family reference
(1AMY, 2QMJ) and shown in the orientation of panel a or b. The number under each subsite label is the number of
units. Grey: the reference protein cut ±5 Å around the median depth of the centroids. Magenta +: catalytic
carboxylates of the reference.
(e) Maltohexaose drawn over subsites −3 to +3, non-reducing end at left. The bar and scissors mark the bond
hydrolysed, C1(−1)–O4(+1). O2, O3 and O6 are named on the −2 unit. Below: a glucose unit with R at O2, O3 and O6,
and the four common starch substituents (hydroxypropyl, acetyl, octenyl succinyl, carboxymethyl), each drawn on
the ring oxygen it is attached to.

## Entries chosen and why

| Family | Entry | Protein, organism | Resolution | Ligand, subsites | RSCC of the sugar units | Mutation | Why |
|---|---|---|---|---|---|---|---|
| GH13 | 6JQB chain A | maltotetraose-forming amylase, *P. saccharophila* STB07 (UniProt P22963) | 1.101 Å | acarbose-derived heptasaccharide 7SA (one residue, 7 rings), −4..+3 | 0.971 (one value for the whole 7SA residue) | none: `_entity.pdbx_mutation` '?', no 'engineered mutation' in `_struct_ref_seq_dif`, RCSB `rcsb_mutation_count` 0 | Of the ten GH13 entries in units_used.tsv with units at all six subsites −3..+3, the only one at ≤ 1.2 Å with every unit on one chain and RSCC ≥ 0.97. The others: 1G94 1.74 Å, min RSCC 0.810; 1W9X 2.1 Å; 1QHO 1.7 Å; 3BC9 1.35 Å, units on three chains, min 0.807; 2GJP 1.9 Å; 26WD 1.9 Å; 1WPC 1.9 Å, min 0.803; 6SAU 1.35 Å, two chains; 6SAO 1.2 Å, min 0.934. Primary citation: Zhang et al., Int. J. Biol. Macromol. 154 (2020) 1303, the target journal |
| GH31 | 3WEO chain A (sugar chain D) | α-glucosidase, *Beta vulgaris* | 1.45 Å | acarviosyl-maltohexaose (AC1 + 6 GLC), −1..+7 | −1/+1 (AC1) 0.987; +2 0.991; +3 0.979; +4 0.983; +5 0.922; +6 0.926; +7 0.805 | none (as above; `rcsb_mutation_count` 0) | Best resolution and highest RSCC of the nine GH31 entries with units at −1..+3. The others: 3W37 1.7 Å; 3WEL 1.84 Å; 2X2I 2.6 Å; 3PHA 2.17 Å, a W169Y mutant; 3TOP 2.88 Å; 3WEM, 3WEN 2.59 Å; 5NN8 2.45 Å. Also the GH31 entry of Figure 2. Primary citation: Tagami et al., J. Biol. Chem. 290, 1796 |

Sources: resolution, protein and organism from database table `entry`. Subsites and RSCC from
`resultados/analysis/33_subsites/units.tsv` (and `44_rule_test/units_used.tsv` for −3..+3). Mutations from the
mmCIF files in `resultados/census/cif` and the RCSB REST endpoint `/core/polymer_entity/<id>/1`.

## Domains (source)

- **6JQB**
  - CATH-B, from the PDBe API `/mappings/cath_b/6jqb` (read 2026-11-24):
    - 6jqbA01, 3.20.20.80, residues 1–360;
    - 6jqbA02, 2.60.40.1180, residues 362–417.
  - PDBe returns no data from `/mappings/cath` or `/mappings/scop` for 6JQB, and the database table
    `sifts_domain` holds only Pfam entries for it.
  - Domain B is not a CATH domain. It was taken as the β3–α3 excursion of the barrel (the usual GH13 definition),
    located on the PyMOL `dss` secondary structure:
    - β3 = 108–113; it ends just before His117, the region I His mapped to 1AMY H92;
    - α3 = 168–185, the helix before β4 (189–192), which ends at the nucleophile Asp193;
    - so **B = 114–167**, **A = 1–113 and 168–361**, **C = 362–418**. Residue 418 is the last one modelled; the
      starch-binding domain of the 530-residue construct is not modelled.
  - The domain boundaries in Zhang et al. (2020) were not checked; the paper could not be reached from here.
    [PENDIENTE]
- **3WEO**: CATH domains from SIFTS (database table `sifts_domain`):
  - 3weoA01, 2.60.40.1760, 38–313 (N);
  - 3weoA02, 3.20.20.80, 314–655 (A);
  - 3weoA03, 2.60.40.1180, 656–756 (C);
  - 3weoA04, 2.60.40.1180, 757–909 (C′).
  - GH31 subdomain B (an insertion in the barrel) is not shown separately.
- Domain letters are placed at the median Cα of each domain on the page, pushed 12 Å clear of the sugar.

## Catalytic residues (source)

- `resultados/calc/03_frames/catalytic_residues_mapped.tsv`:
  - 6JQB: nucleophile Asp193 (1AMY D179), acid/base Glu219 (E204);
  - 3WEO: nucleophile Asp469 (2QMJ D443), acid/base Asp568 (D542).
- `resultados/analysis/44_rule_test/posmap.tsv`: 6JQB Asp294 sits at reference position D289, its Cα 0.15 Å
  from the reference Cα after superposition.

## Counts used (panels c and d; `resultados/analysis/44_rule_test/units_used.tsv`, recounted at run time)

| Family | Subsite | Units | Entries |
|---|---|---|---|
| GH13 | −3 | 43 | 38 |
| GH13 | −2 | 64 | 55 |
| GH13 | −1 | 88 | 74 |
| GH13 | +1 | 51 | 44 |
| GH13 | +2 | 47 | 41 |
| GH13 | +3 | 15 | 14 |
| GH13 | total | 308 | 79 |
| GH31 | −1 | 37 | 18 |
| GH31 | +1 | 16 | 9 |
| GH31 | +2 | 16 | 9 |
| GH31 | +3 | 13 | 9 |
| GH31 | total | 82 | 18 |

These agree with `resultados/analysis/33_subsites/subsite_summary.tsv` (n_units and n_entries at −3..+3).

## Method notes

- **Centroid.** The mean of the six ring atoms of a unit. Rings are found by `ring_units()` of
  `scripts/analysis/a33_subsites.py`, the function that built the subsite table. Each centroid is moved into the
  family frame with the stored superposition of the unit's frame chain. Check: the reference atom moved the same
  way reproduces `x_frame, y_frame, z_frame` of units_used.tsv to within 0.0008 Å (largest deviation over 390
  units; `work/frame_check.txt`). Tables: `work/centroids_GH13.tsv`, `work/centroids_GH31.tsv`.
- **Orientation of panels c and d.** Each centroid goes from the frame to the entry with the entry's own
  superposition, then through the render basis of panel a or b. From `superposition.tsv`:
  - 6JQB→1AMY: 62 of 68 core Cα retained, RMSD 0.623 Å;
  - 3WEO→2QMJ: 83 of 86 retained, RMSD 0.362 Å.
- **Labels in c and d.** Each label is joined by a leader to the median of its cluster. For GH31 +2 and +3 the
  units sit at two positions, so the median falls between them (frame rms radius 1.55 and 2.92 Å in
  subsite_summary.tsv).
- **Views of the whole proteins.** A search over 4000 directions scored three things: the on-page length of the
  −3..+3 chain, no backbone in front of the sugar, and few protein atoms in front of the sugar. The directions
  were then chosen from the montages `work/views_*.png` (deposited frame, pointing toward the viewer):
  - 6JQB: [0.612, −0.389, −0.689];
  - 3WEO: [−0.615, 0.585, 0.529].
- **ChimeraX 1.11.1, offscreen.**
  - Settings: `lighting soft multiShadow 256`, `material dull`, `graphics silhouettes true width 3`, supersample 2.
  - The orthographic camera is set from Python, because `camera ortho` fails offscreen. The label overlay uses
    the same projection, written to `work/whole_*.json`.
  - Sticks: sugar 0.34 Å (outer units 0.20 Å), catalytic side chains 0.28 Å.
- **Chemistry.** RDKit 2025.03.5: rdDepictor with CoordGen, then MolDraw2DCairo with Liberation Sans at 7 pt.
  - Maltohexaose is assembled with `Chem.molzip` from one α-D-glucopyranose template. The template's CIP labels
    (C1 S, C2 R, C3 S, C4 S, C5 R) are checked at run time; the formula is C36H62O31.
  - Octenyl succinyl is drawn as the half-ester given in the request, −C(O)CH2CH(CH2CH=CHC5H11)COOH, with the
    (E)-alkene. OSA also gives the other regio-isomer, with the ester at the carbonyl next to the octenyl branch.

## Problems left

- The 6JQB domain B boundaries come from the secondary structure, not from the primary paper. [PENDIENTE]
- In panel b of A, the 3WEO chain lies under the A/C domains seen through the surface. The −1 to +3 labels sit
  close together but do not overlap.
- In the 6JQB surface there is a through-gap next to the −3/−4 units. It is real geometry: a magenta-background
  test shows background through it. In candidate B it appears as a dark hole.
- B is 7.06 in wide at the tight bounding box, within the IJBM 190 mm double-column width.
