# 01_propka — pKa of the titratable groups of every entry

Calculation 1 of `docs/03_supporting_calculations.md`. Run 2026-11-23.

## What was done

- **Input per entry** (`scripts/calc/c01_propka.py`): model 1 of the cached mmCIF (`resultados/census/cif/<PDB>.cif`),
  hydrogens and waters removed, first alternate conformation kept, protein **and all hetero groups as deposited**
  (sugars, ions, buffers, inhibitors) written to PDB format with gemmi 0.7.5. Chains with names longer than one
  character (7 entries: 6YQ9, 6YQA, 6YQB, 6YQC, 7P2Z, 7P32, 9JU0) are renamed to single letters for PDB format; the
  mapping is in `runs/<PDB>/chainmap.tsv` and the collected tables use the author chain names.
- **PROPKA 3.5.1**, default options: `python -m propka --quiet <PDB>.pdb`. PROPKA reads the hetero groups: their atoms
  enter the desolvation and hydrogen-bond terms, and titratable ligand atoms get their own pKa (e.g. the acarbose amine
  AC1 N4A, PLP N1, acetate). Test on 1QHO: Glu256 (acid/base) 8.69 with ligands, 6.75 with HETATM removed; Asp228
  (nucleophile) 5.71 with ligands, 8.23 without.
- **Collection** (`scripts/calc/c01_collect.py`): the SUMMARY section of every `.pka` file; UniProt numbers from the
  database table `sifts_residue`; contacts from the database (`contact` joined to `instance`, ligand_class in
  ('sugar','other'), partner_kind = 'protein', 4.0 Å); catalytic residues from `../03_frames/catalytic_residues_mapped.tsv`.

**Where it ran.** 525 entries on the laptop (6, then 5 parallel processes); the remaining 120 (the largest structures:
7U39, 7U3B, 7XOI, …) on GUANE, account <cluster-account>, job 33429 (legacy node guane04, 22 parallel processes, 20 min),
with the same script, PROPKA 3.5.1 and gemmi 0.7.5 in a Python 3.9 venv (`~/enzimas_calc123/`), after AFC asked that the
laptop be kept free. 1QHO run on both: identical pKa values.

## Files

| File | One row is |
|---|---|
| `propka_all.tsv` | one titratable group of one entry: PDB id, family, author chain, residue name, author residue number, UniProt accession and number (SIFTS), group type, pKa, model pKa, shift, catalytic role, reference residue and superposition distance, minimum distance to a sugar/'other' ligand and the list of those ligands within 4 Å |
| `propka_near_ligand.tsv` | the same columns, restricted to protein residues within 4.0 Å of a sugar or 'other' instance |
| `runs/<PDB>/` | PROPKA input `<PDB>.pdb`, output `<PDB>.pka`, `propka.log`, `chainmap.tsv` |
| `entries.txt`, `guane_entries.txt` | all 645 entries; the 120 run on GUANE |
| `batch.log`, `guane_batch_33429.log` | per-entry status lines (laptop; GUANE job 33429) |

**Catalytic flag** (`catalytic_role`): `nuc_or_base` or `acid` when the Part 3 reference residue of the family, carried
to the chain by calculation 3, has the reference residue type and its functional atoms lie within 3.0 Å of the
reference ones after superposition; a trailing `?` when the residue found fails one of these tests. Families without a
Part 3 reference (AA13, GH49, GH63, GH97, GT35) and single-entry families carry no flag.

## Counts

- 645 of 645 entries have a PROPKA result; 197 516 titratable groups (Asp 43 398, Glu 37 653, Tyr 31 122, Lys 30 773,
  Arg 30 745, His 15 084, Cys 6 290, N-termini 1 179, C-termini 508, ligand atoms 764).
- **Correction (2026-11-23).** This calculation ran on all 645 entries in scope at the time. 8 of them (8AFE, 8AFH,
  8AFL, 8AFM, 8AHL, 8AIA, 8AIX, 8AJB) were found the same day to be mislabelled as GH63 (cryo-EM crescentin filaments,
  not the YgjK α-glucosidase) and were excluded from the experimental database, which now holds 637 entries
  (`docs/02_experimental_database.md`). Their PROPKA results are still in `propka_all.tsv`/`runs/` but are out of
  scope for any count that cites the 637-entry database; none of them appear in `propka_near_ligand.tsv` or the
  catalytic-flag counts below (they carried no sugar/'other' ligand).
- 186 471 groups carry a UniProt number (SIFTS maps 619 entries).
- `propka_near_ligand.tsv`: 6 347 protein groups within 4.0 Å of a sugar/'other' instance, in 388 entries.
- Catalytic flags: 1 614 rows — nuc_or_base 789 (+13 with `?`), acid 769 (+43 with `?`).
- Mean pKa of the flagged catalytic residues (all chains, including `?`): GH13 nucleophile 6.57 (n = 524), acid/base
  6.27 (n = 541); GH14 base 4.90, acid 4.74 (n = 73); GH15 base 6.12, acid 3.03 (n = 20); GH31 nucleophile 4.81
  (n = 131), acid/base 5.86 (n = 132); GH57 nucleophile 7.69 (n = 54), acid/base 5.78 (n = 46).

## Failures and limits

- No entry failed. Five runs on the laptop were stopped part-way when the laptop was handed back (7U39, 7U3B, 7WJ9,
  7WLG, 7XOI); their partial folders were deleted and they were run on GUANE.
- PROPKA keys residues by name, number and chain only: insertion codes are lost. One entry has insertion codes in a
  protein chain, 1KXT chain B (the VHH antibody fragment, Kabat numbers 82A–C, 100A–G); its groups at those numbers
  may be merged. The enzyme chains are unaffected.
- Ligand groups are named by component, atom and chain; the residue number is taken from the input PDB and is a
  `;`-joined list when the same atom name occurs in several copies of the component in that chain.
- PROPKA treats each structure as deposited: one model, first alternate conformation, no hydrogens, no
  conformational relaxation. How PROPKA 3.5.1 weighs the metal ions (Ca²⁺ at most GH13 sites) was not examined.
