# Figure 2 (v2) — draft captions for candidates A, B, C

Scripts (run in this order, from the project root, with `~/miniconda3/bin/python`):
`scripts/figures/fig2_v2_data.py` → `scripts/figures/fig2_v2_layouts.py A B C`
(`scripts/figures/fig2_v2_struct.py` holds the PyMOL rendering, the view search and the label placement).

Where the numbers come from:
- Residue types per protein (panel c): `resultados/figures/fig2_v2_work/per_protein_counts.tsv`
  (per protein and position: `per_protein.tsv`; every row: `contact_rows.tsv`). Built by `fig2_v2_data.py` with
  the same method as `scripts/analysis/a40_contact_types.py`. Its rows for the Table 1 positions are identical
  to `resultados/analysis/40_contact_types/contact_types.tsv` (977 rows; the script checks this each time it runs).
  Lys182 (+2) is new: `contact_types.tsv` does not hold it.
- Lys182 added and Gly183 dropped: `resultados/analysis/33_subsites/subsite_residues_reference_numbering.tsv`
  (Lys182 at +2: 37 of 41 entries, 0.902; Gly183 at +1: 34 of 44, 0.773).
- Protein = UniProt accession of the chain that carries the residue (`sifts_residue` in
  `resultados/db/enzimas_experimental.sqlite`; `entry.uniprot` when SIFTS has no mapping; the PDB code when
  neither exists: 1W9X, 24ZQ and 7P01 have no accession, so each counts as one protein).
- Engineered mutants left out of the counts: 1QI5, 5TD4, 9UNS, 9XG6 (Asp→Asn at Asp289), 7WJC, 7WJF (D542A).
  The other entries in the counts that carry an engineered mutation according to RCSB (`pdbx_mutation`; list in
  `resultados/figures/fig2_v2_work/rcsb_mutations.tsv`: 1NM9, 1Z32, 3BLK, 3BLP, 3DHP, 8CQF, 3PHA, 1UA7) are
  mutated at residues that are not figure positions, so they stay in.
- Entries shown: resolution from `entry.resolution_A`; RSCC from `resultados/analysis/33_subsites/units.tsv`;
  mutation status from the RCSB Data API (`pdbx_mutation` empty and `rcsb_mutation_count` = 0 for 3BC9, 1G94,
  3WEO); distances and angles computed by `fig2_v2_struct.py` from the deposited mmCIF files
  (`resultados/census/cif/`).

---

## Candidate A (`fig2_v2_A.png/.pdf`; 7.0 × 6.6 in)

**Figure 2.** Each class of contact at the catalytic subsite, shown on a deposited structure, and the residue
types found at every conserved position. (a) GH13, PDB 3BC9 (α-amylase AmyB of *Halothermothrix orenii*,
UniProt B8CZ54, 1.35 Å, bound to acarbose); (b) GH31, PDB 3WEO (sugar beet α-glucosidase, UniProt L0N7E5,
1.45 Å, bound to acarviosyl-maltohexaose). Both are wild-type proteins, and every sugar unit drawn fits the
electron density with RSCC ≥ 0.93 (3BC9: −1 0.935, +1 0.971, +2 0.958; 3WEO: −1 and +1 0.987). Residues are
named by their position in the family reference structure (GH13: 1AMY; GH31: 2QMJ). Where the entry carries a
different residue, the entry's residue is given first and the reference position in parentheses. Side-chain
carbons are coloured by contact class: blue, H (a side-chain N or O lies within 3.5 Å of a sugar O or N, read as
a hydrogen bond); orange, Ar (apolar contact made by an aromatic side chain); pink, Al (apolar contact made by an
aliphatic side chain). Sugar carbons are grey, O red, N blue; the protein chain is drawn as a grey cartoon.
−1, +1, +2 mark the sugar subsites (cleavage between −1 and +1). Dashed lines join each side-chain N or O to its
nearest sugar O or N when that distance is ≤ 3.5 Å. Sugar atoms are labelled where they take part in such a
contact; they are named after the equivalent glucose atom, and N is the glycosidic nitrogen of the inhibitor,
which takes the place of the glycosidic O4. In 3BC9, His92 (3.83 Å) and in 3WEO the nucleophile Asp443 (3.73 Å)
are just beyond 3.5 Å from the sugar, so no line is drawn to them. Insets: the aromatic ring stacked on the −1
ring (Tyr51 in 3BC9, Phe575 in 3WEO), seen from the side; the dotted line joins the two ring centroids, and the
inset gives the centroid–centroid distance and the angle between the ring planes. (c) Residue types at each
conserved position, as the fraction of proteins (one vote per UniProt accession; engineered mutants
excluded), with *n* the number of proteins. Rows are grouped by the class of contact made in most entries;
bar segments are the residue types, majority first, in tints of the class colour; types too narrow to label are
listed to the right of *n* with their protein count. Positions: those of Table 1, without Gly183 (+1; contacted
in 34 of 44 entries, below the 80% threshold) and with Lys182 (+2; contacted in 37 of 41 entries).

Data in panel c (proteins; native entries in brackets):
GH13 H — Asp179 (−1): Asp 42 [69]; Asp289 (−1): Asp 42 [70]; His288 (−1): His 41, Gln 1 [69];
Arg177 (−1): Arg 42 [69]; His92 (−1): His 42 [70]; Glu204 (+1): Glu 25 [41]; Lys182 (+2): Lys 16, Arg 3,
Gln 1 [36]. GH13 Ar — Tyr51 (−1): Tyr 41 [67]. GH13 Al — Phe180 (+1): Ala 14, Phe 4, Thr 4, Tyr 1, Val 1 [40];
Trp206 (+2): Trp 11, Ile 5, Phe 2, Tyr 1, Leu 1 [36].
GH31 H — Asp327, His600, Asp443, Asp542 (−1): 6 of 6 each [16]. GH31 Ar — Trp441: Trp 6 [16]; Trp539:
Trp 5, Phe 1 [16]; Phe575: Phe 4, His 1 [14]. GH31 Al — Ile328: Ile 3, Val 1, Leu 1, Met 1 [15]; Ile364:
Ile 4, Asn 2 [14].

---

## Candidate B (`fig2_v2_B.png/.pdf`; 6.8 × 7.0 in)

Same entries, views and data as A. The structure panels are stacked in the left column, so they print larger
(about 3.7 in wide instead of 3.4 in), and panel c is a single column on the right. The narrower bars use
one-letter codes where the three-letter code does not fit.

**Figure 2.** Caption as for candidate A, with "(c)" placed in the right-hand column. Add to the panel c
sentence: "Residue types are given by three-letter codes, or one-letter codes where the segment is narrow."

---

## Candidate C (`fig2_v2_C.png/.pdf`; 7.0 × 6.8 in)

GH13 is shown on PDB 1G94 (wild-type α-amylase of *Pseudoalteromonas haloplanktis*, UniProt P29957, 1.74 Å,
RSCC −1 and +1 0.966, +2 0.947, bound to an acarbose-derived heptasaccharide). In this entry every H position
draws a line to the sugar, His92 included (3.06 Å to O6), and the variable position Trp206 is Ile, so the
panel itself shows a variable aliphatic position. Panel c is arranged by contact class across both families
(H positions of GH13 and GH31 on the left; Ar and Al on the right), which is the form of the claim.

**Figure 2.** Each class of contact at the catalytic subsite, shown on a deposited structure, and the residue
types found at every conserved position. (a) GH13, PDB 1G94 (α-amylase of *Pseudoalteromonas haloplanktis*,
UniProt P29957, 1.74 Å, with an acarbose-derived heptasaccharide); (b) GH31, PDB 3WEO (sugar beet
α-glucosidase, UniProt L0N7E5, 1.45 Å, with acarviosyl-maltohexaose). Both are wild-type proteins, and every
sugar unit drawn fits the electron density with RSCC ≥ 0.94 (1G94: −1 and +1 0.966, +2 0.947; 3WEO: −1 and
+1 0.987). Residues are named by their position in the family reference structure (GH13: 1AMY; GH31: 2QMJ).
Where the entry carries a different residue, the entry's residue is given first and the reference position in
parentheses. Side-chain carbons are coloured by contact class: blue, H (a side-chain N or O lies within 3.5 Å of
a sugar O or N, read as a hydrogen bond); orange, Ar (apolar contact made by an aromatic side chain); pink, Al
(apolar contact made by an aliphatic side chain). Sugar carbons are grey, O red, N blue; the protein chain is
drawn as a grey cartoon. −1, +1, +2 mark the sugar subsites (cleavage between −1 and +1). Dashed lines join each
side-chain N or O to its nearest sugar O or N when that distance is ≤ 3.5 Å. Sugar atoms are labelled where they
take part in such a contact; they are named after the equivalent glucose atom, and N is the glycosidic nitrogen
of the inhibitor, which takes the place of the glycosidic O4. In 3WEO the nucleophile Asp443 is 3.73 Å from the
nearest sugar O, so no line is drawn to it. Insets: the aromatic ring stacked on the −1 ring (Tyr51 in 1G94,
Phe575 in 3WEO), seen from the side; the dotted line joins the two ring centroids, and the inset gives the
centroid–centroid distance and the angle between the ring planes (1G94: 4.9 Å, 21°; 3WEO: 5.2 Å, 34°). (c)
Residue types at each conserved position, as the fraction of proteins (one vote per UniProt accession;
engineered mutants excluded), with *n* the number of proteins; left, positions whose side chain
hydrogen-bonds the sugar in most entries (H); right, positions that make an apolar contact (Ar, aromatic; Al,
aliphatic). Bar segments are the residue types, majority first, in tints of the class colour; types too narrow
to label are listed to the right of *n* with their protein count. Positions: those of Table 1, without Gly183
(+1; contacted in 34 of 44 entries, below the 80% threshold) and with Lys182 (+2; contacted in 37 of 41
entries).

(The data in panel c are the same as for candidate A.)

---

## Points the author should know before choosing

1. **Lys182 is not invariant.** Counted one vote per protein, it is Lys in 16 of 20 proteins, Arg in 3
   (1RPK barley AMY1 P00693; 4E2O *Geobacillus thermoleovorans* A0ABF7PFV7; 6JQB *Roseateles saccharophilus*
   P22963) and Gln in 1 (1UH3 *Thermoactinomyces vulgaris* TVA I, Q60053). Its side chain hydrogen-bonds the
   +2 unit in 19 of the 20 proteins. None of the four entries carries an engineered mutation (RCSB).
2. **Trp441 and Trp539 (GH31) do not stack.** In 3WEO their rings are at 76° and 77° to the −1 ring
   (centroid distances 6.0 and 6.9 Å): edge-on contacts. Only Phe575 lies roughly parallel (34°, 5.2 Å). The
   captions therefore define Ar as "apolar contact made by an aromatic side chain", which is the definition of
   Table 1, and not as "stacking". The insets show the two real stacks (Tyr51 in GH13, Phe575 in GH31).
3. **Counted per protein, the majority type at Trp206 changes.** Per entry (Table 1) it is Ile 20, Trp 13; per
   protein it is Trp 11, Ile 5, because the 20 Ile entries come mostly from human amylases deposited many
   times. Phe180 per protein: Ala 14, Phe 4, Thr 4, Tyr 1, Val 1 (per entry Ala 31).
4. **His92 in 3BC9** (candidates A, B) is 3.83 Å from O6, so no line is drawn to an H residue; candidate C's
   1G94 avoids this.
5. Protein = UniProt accession from SIFTS. Three entries have no accession and count as one protein each
   (1W9X, 24ZQ, 7P01). 1W9X (*Sutcliffiella halmapala* α-amylase) may be the same protein as 2GJP, which SIFTS
   maps to P19571; if so, the +1/+2 counts would each drop by one. Not checked.
