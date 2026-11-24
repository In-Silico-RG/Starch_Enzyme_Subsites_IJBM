# Figure 1 candidates (protein overview, replaces the methyl β-maltoside figure) — 2026-11-24

Files: `fig1_protein_A.{png,pdf}`, `fig1_protein_B.{png,pdf}` (7.0 in wide, 300 dpi). Working renders and label
anchors: `work/`. Scripts: `scripts/figures/fig1_protein_common.py`, `fig1_protein_struct.py` (PyMOL renders),
`fig1_protein_frame.py` (centroids in the family frame), `fig1_protein_chem.py` (RDKit depictions),
`fig1_protein_layouts.py` (assembly). Run order: frame → struct → layouts (chem is called by layouts).

## What each candidate shows

**Candidate A (proteins first).** a, GH13 6JQB whole protein; b, GH31 3WEO whole protein (cartoon coloured by
domain under a transparent molecular surface); c, GH13 frame; d, GH31 frame; e, catalytic centre of 6JQB;
f, chemistry (maltohexaose spanning −3..+3 with the scissile bond; glucose unit with R at O2/O3/O6 and the four
substituents).

**Candidate B (question first).** a, chemistry (same content as A-f); b–d, GH13: 6JQB whole protein as an opaque
molecular surface coloured by domain (the cleft reads as a groove holding the chain), GH13 frame in the same
orientation, catalytic centre; e–g, the same for GH31 3WEO.

## Draft caption (candidate A lettering)

**Figure 1.** Starch-degrading glycoside hydrolases, their subsites and the substituents of functionalised
starch. (a) GH13 maltotetraose-forming amylase of *Pseudomonas saccharophila* STB07 (PDB 6JQB, 1.10 Å) with the
acarbose-derived heptasaccharide bound from subsite −4 to +3; domains A ((β/α)8 barrel), B and C; units at
subsites −3 to +3 coloured by subsite (the −4 unit grey); catalytic Asp193 (nucleophile), Glu219 (acid/base) and
Asp294 (transition-state stabiliser) in dark grey. (b) GH31 sugar beet α-glucosidase (PDB 3WEO, 1.45 Å) with
acarviosyl-maltohexaose bound from −1 to +7; N-terminal β-sandwich, (β/α)8 barrel (A), proximal and distal
C-terminal domains; units −1 to +3 coloured, +4 to +7 grey; catalytic Asp469 (nucleophile) and Asp568
(acid/base). (c, d) Ring centroids of every density-supported sugar unit (RSCC ≥ 0.80) assigned to subsites −3
to +3 across the GH13 (308 units, 79 entries) and GH31 (82 units, 18 entries) structures, superposed on the
family reference (1AMY, 2QMJ) and shown in the orientation of panels a and b; grey: section of the reference
protein within ±5 Å of the median depth of the centroids; +: catalytic carboxylates of the reference. Number of
units per subsite in the legends. (e) Catalytic centre of 6JQB: units at −1 and +1 and the three catalytic
residues. (f) Maltohexaose drawn over subsites −3 to +3 (non-reducing end left); bar and scissors: the
glycosidic bond C1(−1)–O4(+1) hydrolysed; O2, O3, O6 named on the −2 unit; a glucose unit with R at O2, O3 and
O6 and the four common starch substituents (hydroxypropyl, acetyl, octenyl succinyl, carboxymethyl), each drawn
on the ring oxygen it esterifies or etherifies. Subsite colours are the same in every panel.

## Entries chosen and why

| Family | Entry | Protein, organism | Resolution | Ligand, subsites | RSCC of the sugar units | Mutation | Why |
|---|---|---|---|---|---|---|---|
| GH13 | 6JQB chain A | maltotetraose-forming amylase (MTA), *P. saccharophila* STB07 (UniProt P22963) | 1.101 Å | acarbose-derived heptasaccharide 7SA (one residue, 7 rings), −4..+3 | 0.971 (one value for the whole 7SA residue; RCSB/PDBe per-residue RSCC) | none: `_entity.pdbx_mutation` '?', no 'engineered mutation' in `_struct_ref_seq_dif`, RCSB `rcsb_mutation_count` 0 | the only GH13 entry of units_used.tsv with all six subsites −3..+3 on one chain at better than 1.2 Å with RSCC ≥ 0.97 (the other six-subsite entries: 1G94 1.74 Å, min RSCC 0.810; 1W9X 2.1 Å; 1QHO 1.7 Å; 3BC9 1.35 Å, units on three chains, min 0.807; 2GJP 1.9 Å; 26WD 1.9 Å; 1WPC 1.9 Å, min 0.803; 6SAU 1.35 Å, two chains; 6SAO 1.2 Å, min 0.934). Primary citation: Zhang et al., Int. J. Biol. Macromol. 154 (2020) 1303 — the target journal |
| GH31 | 3WEO chain A (sugar chain D) | α-glucosidase (SBG), *Beta vulgaris* | 1.45 Å | acarviosyl-maltohexaose (AC1 + 6 GLC), −1..+7 | −1/+1 (AC1) 0.987; +2 0.991; +3 0.979; +4 0.983; +5 0.922; +6 0.926; +7 0.805 | none: as above, `rcsb_mutation_count` 0 | best resolution and highest RSCC of the nine GH31 entries with units at −1..+3 (3W37 1.7 Å, 3WEL 1.84 Å, 2X2I 2.6 Å, 3PHA 2.17 Å and W169Y mutant, 3TOP 2.88 Å, 3WEM/3WEN 2.59 Å, 5NN8 2.45 Å); already the GH31 entry of Figure 2. Primary citation: Tagami et al., J. Biol. Chem. 290, 1796 |

Sources: resolution, protein, organism — database table `entry`; RSCC and subsites — `resultados/analysis/33_subsites/units.tsv`
(and `44_rule_test/units_used.tsv` for −3..+3); mutations — mmCIF in `resultados/census/cif` and RCSB REST
`/core/polymer_entity/<id>/1`.

## Domains (source)

- **6JQB**: CATH-B via PDBe API `/mappings/cath_b/6jqb` (read 2026-11-24): 6jqbA01 3.20.20.80, residues 1–360;
  6jqbA02 2.60.40.1180, residues 362–417. CATH and SCOP have no mapping for 6JQB (PDBe `/mappings/cath`,
  `/mappings/scop` empty); the database `sifts_domain` holds only Pfam for it. Domain B is not a CATH domain; it
  was taken as the β3–α3 excursion of the barrel (the usual GH13 definition), located on the PyMOL `dss`
  secondary structure of 6JQB: β3 = 108–113 (ends just before His117, the region I His mapped to 1AMY H92),
  α3 = 168–185 (the helix before β4 189–192, which ends at the nucleophile Asp193) → **B = 114–167**,
  **A = 1–113 + 168–361**, **C = 362–418** (418 = last modelled residue; the starch-binding domain of the
  530-residue construct is not modelled). The literature boundaries of Zhang et al. (2020) were not checked
  (paper not accessible from here): [PENDIENTE] compare with their domain definition.
- **3WEO**: CATH from SIFTS (database table `sifts_domain`): 3weoA01 2.60.40.1760, 38–313 (N-terminal
  β-sandwich); 3weoA02 3.20.20.80, 314–655 ((β/α)8 barrel, domain A); 3weoA03 2.60.40.1180, 656–756 (proximal
  C); 3weoA04 2.60.40.1180, 757–909 (distal C). GH31 subdomain B (an insertion in the barrel) is not separated.

## Catalytic residues (source)

`resultados/calc/03_frames/catalytic_residues_mapped.tsv`: 6JQB nucleophile Asp193 (reference 1AMY D179),
acid/base Glu219 (E204); 3WEO nucleophile Asp469 (2QMJ D443), acid/base Asp568 (D542).
`resultados/analysis/44_rule_test/posmap.tsv`: 6JQB Asp294 at reference position D289 (Cα 0.15 Å from the
reference Cα after superposition).

## Counts used (panel c/d legends; `resultados/analysis/44_rule_test/units_used.tsv`, recounted at run time)

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

These agree with `resultados/analysis/33_subsites/subsite_summary.tsv` (n_units, n_entries at −3..+3).

## Method notes

- Centroid = mean of the six ring atoms of the unit (ring found by `ring_units()` of
  `scripts/analysis/a33_subsites.py`, the function that made the subsite table), moved into the family frame
  with the stored superposition of the unit's frame chain. Check: the reference atom moved the same way
  reproduces `x_frame, y_frame, z_frame` of units_used.tsv to 0.0008 Å (max over 390 units;
  `work/frame_check.txt`). Table: `work/centroids_GH13.tsv`, `work/centroids_GH31.tsv`.
- Panels c/d use the orientation of panels a/b: frame → entry with the entry's own superposition
  (6JQB→1AMY: 62 of 68 core Cα retained, RMSD 0.623 Å; 3WEO→2QMJ: 83 of 86, 0.362 Å; superposition.tsv),
  then the render basis.
- Views of the whole proteins: search over 4000 directions (on-page length of the −3..+3 chain, sugar not
  covered by backbone in front, few protein atoms in front of the sugar); chosen 6JQB [0.612, −0.389, −0.689],
  3WEO [−0.615, 0.585, 0.529] (deposited frame, toward the viewer) from the montages `work/views_*.png`.
  Close-ups: search over 20 000 directions (both rings tilted 0.6–0.97 to the page, the −1→+1 step not
  foreshortened, label anchors apart, no residue in front of a sugar atom); chosen 6JQB [0.8358, −0.0189, 0.5487],
  3WEO [−0.0125, −0.153, 0.9882].
- Rendering: PyMOL (orthoscopic, ray-traced, ray_trace_mode 1 outlines, white background); sticks 0.30 Å in the
  whole-protein panels, 0.15/0.13 Å (sugar/residues) in the close-ups. Chemistry: RDKit 2025.03.5,
  rdDepictor with CoordGen, MolDraw2DCairo, DejaVu Sans at 7 pt; maltohexaose assembled with `Chem.molzip` from
  one α-D-glucopyranose template (CIP C1 S, C2 R, C3 S, C4 S, C5 R checked at run time; formula C36H62O31).
- Octenyl succinyl is drawn as the half-ester written in the request, −C(O)CH2CH(CH2CH=CHC5H11)COOH, (E)-alkene;
  the other regio-isomer (ester at the carbonyl next to the octenyl branch) also forms from OSA.
