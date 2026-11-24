# PROPKA review: the catalytic pKa order in retaining families (checked 2026-11-23)

**Problem.** The textbook mechanism of a retaining glycosidase (already cited in the manuscript,
Davies1995, McCarter1994) expects the nucleophile to stay deprotonated (a low pKa, so it can attack)
and the general acid/base to have an elevated pKa (so it can first donate a proton, then pick one
up again from the incoming water). Do the PROPKA pKa values computed for this database's reference
structures actually show that pattern?

**Data.** `resultados/calc/01_propka/propka_near_ligand.tsv` (6 347 rows, one per titratable
residue within 4 Å of a sugar/other ligand), restricted to rows flagged `catalytic_role`, cross-checked
against the official residue assignment in `resultados/part3/part3_step1_catalytic_residues.tsv`.

**Method.** For the three retaining families with a confirmed catalytic pair (GH13, GH31, GH57),
the residue type behind each role was confirmed first (not assumed): GH13 nucleophile=Asp,
acid/base=Glu; GH31 nucleophile=Asp, acid/base=Asp (a different Asp); GH57 nucleophile=Glu,
acid/base=Asp — all three match `part3_step1_catalytic_residues.tsv` exactly. The PROPKA pKa of
every such residue, across every reference structure that has one, was then compared, role against
role, within each family.

**Result.**

| Family | Nucleophile (n) | median pK$_a$ | Acid/base (n) | median pK$_a$ | Textbook order? |
|---|---|---|---|---|---|
| GH13 | Asp (203) | **7.46** | Glu (196) | 5.30 | **reversed** |
| GH31 | Asp (81) | 5.16 | Asp (88) | 5.78 | matches (weakly) |
| GH57 | Glu (35) | **7.41** | Asp (27) | 5.39 | **reversed** |

In two of the three retaining families with a confirmed pair, the computed nucleophile pK$_a$ is
*higher*, not lower, than the acid/base's — the opposite of the simple textbook picture.

**Reading (mine, kept apart from the numbers).** This could mean several different things, and I
have not settled which:
1. PROPKA predicts a single, static, resting-state pK$_a$ from one crystal geometry. The textbook
   rule is about the *reaction coordinate* — the nucleophile must be able to reach a low-pKa,
   reactive state at the moment of attack, and the acid/base must cycle between two states during
   the two chemical steps. A static resting-state number does not have to honour that order; it is
   not necessarily wrong for it to differ.
2. PROPKA is known in the literature to struggle with buried carboxylates near a bound ligand,
   where the empirical pKa model's assumptions are most stretched; this could be a real limitation
   of the tool on exactly these residues, not a fact about the enzymes.
3. The reference structures compared here are not all the same kind (some bound, some not); no
   check was done for whether a ligand nearby systematically shifts one role's pKa more than the
   other's.

**Checks.** The residue-type-to-role mapping was independently confirmed against
`part3_step1_catalytic_residues.tsv` for all three families before comparing pKa (not assumed from
column names). GH14/GH15 (inverting families) were excluded from this comparison: an inverting
mechanism has no covalent nucleophile, so the textbook rule being tested here does not apply to
them.

**Limits.** This is a comparison of PROPKA numbers, not a validated mechanistic claim — I have not
checked whether this reversal (or something like it) is already reported and explained in the
glycosidase literature, nor run a second pKa method for comparison. Separately found while checking
this table: 1 041 of 6 347 rows have pKa > 14, almost all (1 017) for Arg/Tyr/Lys (PROPKA correctly
reporting these as effectively non-titratable near neutral pH in this environment, not an error);
24 Cys rows report PROPKA's standard 99.99 flag for a disulfide-bonded cysteine (also not an error).
Neither of those affects the catalytic-residue comparison above, since none of the flagged rows are
Asp/Glu catalytic residues.

**AFC's call needed:** whether this is worth a line in the manuscript (as an open, honestly-flagged
question) or should wait for a literature check on whether resting-state pKa reversals of this kind
are already known for GH13/GH57-type retaining mechanisms.
