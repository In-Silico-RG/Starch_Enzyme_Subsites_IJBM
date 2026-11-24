# Two checks of the subsite −1 grip finding (2026-11-23, late)

For the short draft `manuscript_short/grip_finding.tex`.

## Check 1: residue types one ring away were miscounted

**Problem.** The draft (and `manuscript_JCIM/atlas_jcim.tex`, Section 3.3 and the Discussion)
said that at the GH13 position Phe180, at subsite +1, the residue is "Ala in 33, Thr in 11,
Val in 9, Leu in 8, Phe in 5". That sums to 66, while the same text gives 44 entries at +1.

**Data.** `manuscript_JCIM/figures/subsite_residue_type_check.tsv` (one row per subsite, entry,
reference position; columns family, subsite, pdb, ref, res, agree).

**Method.** Rows filtered to GH13 and to one reference position, then split by subsite; one
vote per entry.

**What is read off.**

| Position | Subsite | Entries | Residue types |
|---|---|---|---|
| F180 | +1 | 43 | Ala 31, Phe 4, Thr 4, Val 3, Tyr 1 |
| F180 | −1 | 28 | Val 8, Leu 8, Thr 7, Ala 3, Phe 1, Met 1 |
| W206 | +2 | 37 | Ile 20, Trp 13, Phe 2, Tyr 1, Leu 1 |
| W206 | +1 | 20 | Trp 8, Ile 7, Phe 4, Tyr 1 |
| G183 | +1 | 34 | His 33, Tyr 1 |

The old breakdown pooled F180 contacts from subsites +1 and −1 (and a −1~N1J row).

**Reading (mine).** One ring out, each position has one dominant type (54–97 % of entries),
often different from the reference enzyme (barley AMY2), not "freely varying". The contrast
with −1 (reference type in every entry at the catalytic positions, Tyr51, His92, Arg177;
70/71 at His288) still holds, narrower than first written.

**Not done.** One vote per entry; lineage weighting pending (`grip_independence_check.md`).
The 80-page draft keeps the pooled numbers: AFC said not to touch it.

## Check 2: the grip is not an acarbose artefact

**Problem.** Acarbose-type inhibitors (transition-state mimics) supply many subsite −1 units;
the conserved contacts could reflect how inhibitors bind rather than how substrate binds.

**Data.** `resultados/analysis/33_subsites/units.tsv`, `contacts_mapped.tsv`.

**Method.** GH13 units at −1 with in_tables = True and RSCC ≥ 0.80 (the filter of
`scripts/analysis/a33_subsites.py`, line 354): 88 units, 74 entries. An entry is an
"inhibitor/pseudo-sugar complex" if any labelled unit in it is not GLC/BGC, or it is 1QHO
(valienamine deposited as GLC); otherwise "glucose-only". Contacts (asymmetric unit only)
assigned to units by (pdb, instance_id, ring). The one assumption: the entry-level label
captures whether the −1 unit belongs to an inhibitor-derived chain.

**What is read off.**

| Group | Entries | D289 | H92 | D179 | R177 | Y51 | H288 |
|---|---|---|---|---|---|---|---|
| inhibitor/pseudo-sugar | 50 | 50 | 50 | 50 | 49 | 47 | 49 |
| glucose-only | 24 | 24 | 24 | 23 | 22 | 24 | 22 |

50 + 24 = 74, the table total. Glucose-only entries: 1J0I, 1QHP, 1QI5, 1WO2, 2FHB, 2FHC,
2FHF, 2GVY, 3A4A, 3BCD, 3L2L, 3WDI, 3WDJ, 3WY2, 5A2B, 5A2C, 5DO8, 5TD4, 6A0K, 6A0L, 6JFJ,
6JFX, 7P07, 9XG6.

**Check.** A first pass without the RSCC filter gave 117 units in 95 entries; corrected.

**Not done.** Glucose-only entries are not separated by protein; some may be mutants.

## Check 3: GH31 residue identity at −1, and contact type at every conserved position

**Problem.** The draft calls the GH31 grip "chemically fixed"; and a proposed reading was that
hydrogen-bonding residues are invariant while packing residues vary.

**Data.** `subsite_residue_type_check.tsv`; the deposited mmCIF files in `resultados/census/cif`;
units and asymmetric-unit contacts from `33_subsites`.

**Method.** `scripts/analysis/a40_contact_types.py`: for each density-supported unit at the
subsite and each reference position, the mapped residue and the unit's contacting ligand atoms
are read from model 1, first altloc. Side-chain polar = a side-chain N/O within 3.5 Å of a ligand
N/O; backbone polar = only a main-chain N/O; apolar otherwise. One vote per entry (strongest kind).
The one assumption: a 3.5 Å N/O–N/O distance stands for a hydrogen bond (no angle test).

**What is read off** (`resultados/analysis/40_contact_types/summary.tsv`):

| Fam | Sub | Pos | n | SC polar | BB polar | apolar | Types |
|---|---|---|---|---|---|---|---|
| GH13 | −1 | D179 | 73 | 69 | 0 | 4 | Asp 73 |
| GH13 | −1 | D289 | 74 | 73 | 0 | 1 | Asp 70, Asn 4 |
| GH13 | −1 | H288 | 71 | 71 | 0 | 0 | His 70, Gln 1 |
| GH13 | −1 | H92 | 74 | 59 | 0 | 15 | His 74 |
| GH13 | −1 | R177 | 71 | 69 | 0 | 2 | Arg 71 |
| GH13 | −1 | Y51 | 71 | 2 | 4 | 65 | Tyr 71 |
| GH13 | +1 | E204 | 44 | 44 | 0 | 0 | Glu 44 |
| GH13 | +1 | G183 | 34 | 34 | 0 | 0 | His 33, Tyr 1 |
| GH13 | +1 | F180 | 43 | 0 | 0 | 43 | Ala 31, Phe 4, Thr 4, Val 3, Tyr 1 |
| GH13 | +2 | W206 | 37 | 0 | 6 | 31 | Ile 20, Trp 13, Phe 2, Tyr 1, Leu 1 |
| GH31 | −1 | D327 | 18 | 18 | 0 | 0 | Asp 18 |
| GH31 | −1 | H600 | 18 | 18 | 0 | 0 | His 18 |
| GH31 | −1 | D443 | 18 | 11 | 0 | 7 | Asp 18 |
| GH31 | −1 | D542 | 18 | 16 | 0 | 2 | Asp 16, Ala 2 (7WJC, 7WJF) |
| GH31 | −1 | W441 | 18 | 0 | 0 | 18 | Trp 18 |
| GH31 | −1 | W539 | 18 | 0 | 0 | 18 | Trp 15, Phe 3 (7WJB/C/F) |
| GH31 | −1 | F575 | 16 | 0 | 0 | 16 | Phe 13, His 3 (7WJB/C/F) |
| GH31 | −1 | I328 | 17 | 0 | 0 | 17 | Ile 7, Leu 4, Val 3, Met 3 |
| GH31 | −1 | I364 | 16 | 5 | 0 | 11 | Ile 9, Asn 7 |

GH31 glucose-only entries at −1: 5 (5KZW, 5KZX, 7WJB, 7WJC, 7WJF), all nine positions 5/5;
inhibitor entries: 13, positions 11–13/13.

**Reading (mine).** "Polar invariant, apolar variable" fails (Tyr51, Trp441 are apolar and
invariant). Three classes fit every position: side-chain hydrogen-bonding residues invariant;
aromatic stacking residues invariant or aromatic; aliphatic packing residues interchangeable.
The same rule holds at +1. "Tight at the cut, loose around it" may follow from where each kind
of residue sits, rather than from distance to the cut. The GH31 "chemically fixed" wording in
the draft is wrong for I328 and I364.

**Not done.** No angle test; stacking inferred from residue type, not ring geometry; one vote
per entry; the 4 Asn at D289 not checked for mutants; GH31 glucose-only set is 2 proteins with
2 catalytic mutants.

## Check 4: the contact-type rule in GH14 and GH15 (folds documented as unrelated to GH13)

**Problem.** GH13 and GH31 may share a remote ancestor (`gh13_gh31_clan_check.md`), so they
cannot test whether the three-class rule (hydrogen-bonding residues invariant, aromatic stackers
aromatic, aliphatic packers interchangeable) arose independently.

**Data.** `scripts/analysis/a40_contact_types.py GH14 GH15`: every position contacted in at least
half the entries at subsites −2, −1, +1; output `40_contact_types/all_summary.tsv`,
`all_contact_types.tsv`. Protein identities from UniProt (primary accession of `entry.uniprot`,
which holds several accessions separated by ";" for GH14).

**Lineages.** GH14: barley P16098 and soybean P10538, 0.66 identity (over the shorter sequence);
*B. cereus* P36924, 0.29 and 0.32 to them: two lineages. GH15: *S. fibuligera* P08017 and
*A. resinae* Q03045, 0.34; Q12537 returned an empty record from UniProt (unresolved).

**What is read off.** GH14 (16 entries): side-chain hydrogen-bonders Asp53, His93, Arg420,
Glu380, Lys295, Glu186 invariant; aromatic Trp55 invariant; aliphatic Leu419, Leu20, Met17,
Ala184, Ile89 invariant; Thr342 Thr/Ser/Ala. GH15 (4 entries): hydrogen-bonders Asp55, Arg54,
Glu179, Glu180, Arg305 invariant; aromatic Tyr48, Trp52, Trp317, Trp417, Trp120, Tyr311
invariant; aliphatic Leu177, Leu415, Ala39 invariant.

**Reading (mine).** The (α/α)6 fold of GH15 builds its −1 site from the same three kinds of
residue as GH13 and GH31. The first two classes are invariant there too. The third class
(aliphatic interchangeable) cannot be tested with two lineages per family: no variation is
expected to show with so few. The full rule rests on GH13 and GH31.

**Not done.** GH15 third protein unresolved; GH14 +1 has 3 entries only (in the table, not read).

## Check 5 (2026-11-24): control on sugars bound away from the catalytic cleft

**Problem.** Residues that hydrogen-bond a sugar may be conserved anywhere sugars bind; then the
invariance at subsite −1 would not point to catalysis. (Quirke & Crich 2020 used lectins and
non-catalytic sites as the equivalent control.)

**Data.** `scripts/analysis/a41_control_noncatalytic.py` (MIN_ENTRIES=4, MIN_PROT=2); output
`resultados/analysis/41_control_noncatalytic/control_positions.tsv`, `control_entries.tsv`.

**Method.** GH13 density-supported sugar units with no subsite, > 12 Å from the catalytic
midpoint (365), keeping only those touching no position that occurs in the GH13 subsite tables
(296; a first pass without this filter included catalytic residues, from sugars in the cleft of
another chain or in shifted poses). Contact kind as in a40. Positions contacted in ≥ 4 entries
from ≥ 2 proteins: 32. One residue type per protein (most common in its entries). Catalytic
hydrogen-bonding positions (GH13 D179, D289, H288, R177, H92 at −1; E204, G183 at +1) subsampled
200 times to each control position's protein count.

**What is read off.**

| Set | Positions | All proteins same residue type |
|---|---|---|
| control, side-chain hydrogen bond | 15 (2–5 proteins each) | 7 |
| control, aromatic | 7 | 1 |
| control, aliphatic | 7 | 3 |
| catalytic hydrogen-bonders, subsampled to the control counts | 7 | 94 % of draws |

Catalytic hydrogen-bonders, one vote per protein: D179 46/46, D289 43/46, H288 43/44, R177 44/44,
H92 46/46, E204 28/28, G183 17/18.

**Reading (mine).** A side chain that hydrogen-bonds a sugar is fixed when it does so at the
catalytic subsite and only about half the time at other sugar sites of the same proteins, at
equal numbers of proteins. The invariance is not a property of sugar binding as such.

**Not done.** Control positions rest on 2–5 proteins each; GH31 has 16 far units (too few).
Protein counts here use the first UniProt accession or the PDB code when empty (46 at D179,
against 42 in `grip_independence_check.md`).
