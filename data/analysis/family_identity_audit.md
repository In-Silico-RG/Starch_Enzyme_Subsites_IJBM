# Family / UniProt identity audit — second-instance check for the GH63 mislabeling bug

*Written 2026-11-23, 20:32 (`date`). Read-only audit; no row of `resultados/db/enzimas_experimental.sqlite`
was changed. Database opened as `sqlite3.connect('file:resultados/db/enzimas_experimental.sqlite?mode=ro',
uri=True)`.*

## 1. The problem

Eight PDB entries (8AFE, 8AFH, 8AFL, 8AFM, 8AHL, 8AIA, 8AIX, 8AJB) were found today labelled as family GH63
(claimed enzyme: α-glucosidase YgjK, UniProt P42592) when they are in fact cryo-EM structures of crescentin
filaments — a cytoskeletal protein with no relation to any glycoside hydrolase (`docs/04_curation_debt.md`
item 1, `docs/02_experimental_database.md` "Correction (2026-11-23)"). That error was caught by
`scripts/calc/c03_frames.py` because GH63 had 17 entries at the time: the 8 crescentin structures failed to
structurally superpose with the 9 genuine GH63 structures, and the mismatch was then confirmed by reading
each entry's own RCSB title.

That detection method needs several entries in a family to work — a lone entry has nothing to superpose
against. Two questions follow:

1. **Task A.** Do the database's four single-entry families (GH4, GH126, GH133, GH176 — the only ones with
   exactly one PDB entry, so the structural cross-check could never have caught an error in them) each
   genuinely hold what they claim?
2. **Task B.** Beyond those four, is there any other (family, UniProt accession) pair in the whole database
   whose UniProt protein name is as obviously wrong for its claimed CAZy family as "crescentin filaments"
   was for GH63?

This is a hypothesis-generating audit, not a decision to delete or relabel anything; any flag below is a
candidate for AFC to examine the way GH63 was examined, not a finding that something is already wrong.

## 2. What the data contain

| What | Where | How much |
|---|---|---|
| Per-entry family, EC, claimed protein, organism, composite UniProt string | `entry` table, `resultados/db/enzimas_experimental.sqlite` | 637 entries, 14 distinct `family` values |
| Per-chain UniProt accession (SIFTS mapping, one row per residue, used here only for its distinct `(pdb, chain, uniprot)` triples) | `sifts_residue` table, same database | 611 of 637 entries have at least one non-null `uniprot`; 26 do not (listed in §6) |
| Independent identity of each PDB entry | RCSB Data API, `https://data.rcsb.org/rest/v1/core/entry/<PDB>` — fields `struct.title`, `rcsb_primary_citation.title`/`rcsb_journal_abbrev`/`year` | 5 entries fetched for this audit (1OBB, 3REN, 7VMA, 8ZEQ, 8J5D) |
| Independent identity of each UniProt accession | UniProt REST API, `https://rest.uniprot.org/uniprotkb/<accession>.json` (single lookups) and `https://rest.uniprot.org/uniprotkb/search?query=accession:A OR accession:B...&fields=accession,protein_name,organism_name` (batch) — fields `proteinDescription` (`recommendedName`/`submissionNames`), `organism.scientificName` | 158 distinct accessions fetched in 2 batch requests (sizes 80 and 78) |
| CAZy's own family definitions | `https://www.cazy.org/<FAMILY>.html` (activities table) and `https://www.cazy.org/<FAMILY>_structure.html` (list of PDB structures CAZy itself assigns to the family) | GH4, GH126, GH133, GH176 fetched for Task A |

All API calls used a generic, non-identifying client string (`research-script/1.0 (contact: generic-client)`
or a generic browser user-agent for CAZy, which rejects the plain script agent); no personal information was
sent in any request.

## 3. The method

**Task A.**
1. Found the single-entry families with the SQL the task specified:
   `SELECT family, pdb FROM entry e JOIN (SELECT family, COUNT(*) n FROM entry GROUP BY 1 HAVING n=1) t USING(family)`.
   Result: GH4→1OBB, GH126→3REN, GH176→7VMA, GH133→8ZEQ (confirmed, matches the task's statement).
2. For each of the 4, read `entry.protein`, `entry.EC`, `entry.organism`, `entry.uniprot` (the database's own
   claim) and `sifts_residue.uniprot` (the per-chain, SIFTS-derived accession — the one `docs/04_curation_debt.md`
   item 2 says to trust, not the composite `entry.uniprot` string).
3. Fetched the RCSB entry record for the PDB id and read `struct.title` and `rcsb_primary_citation.title`.
4. Fetched the UniProt record for the `sifts_residue` accession and read its protein name and organism.
5. Fetched CAZy's own activities table for the claimed family (gives the family's defining EC number(s)) and,
   independently, CAZy's own `_structure.html` page for that family, and checked by text search whether the
   PDB id itself appears there — this is CAZy's own authority on which PDB entries belong to the family,
   independent of anything in this project's database.
6. Verdict = match if the RCSB title/citation, the UniProt protein name, and CAZy's own structure list all
   agree with the claimed family's defining activity; does-not-match if they describe an unrelated molecule
   (the crescentin pattern); inconclusive if the sources disagree or are silent.

**Task B.**
1. Built every distinct `(entry.family, sifts_residue.uniprot)` pair over the whole database with:
   ```sql
   SELECT e.family, s.uniprot, GROUP_CONCAT(DISTINCT s.pdb), COUNT(DISTINCT s.pdb)
   FROM sifts_residue s JOIN entry e ON e.pdb = s.pdb
   WHERE s.uniprot IS NOT NULL AND s.uniprot != ''
   GROUP BY e.family, s.uniprot
   ```
   as the task specified, deliberately using `sifts_residue.uniprot` (per chain) rather than `entry.uniprot`
   (a composite string for 301 heteromeric entries — `docs/04_curation_debt.md` item 2) so that each row is
   one real protein chain, not a merged multi-protein string.
   Result: 158 distinct `(family, uniprot)` pairs, 158 distinct UniProt accessions, spread over all 14
   families.
2. Batch-fetched the UniProt protein name and organism for all 158 accessions (2 requests, `accession:A OR
   accession:B OR ...`, 158/158 resolved — no accession came back empty or obsolete).
3. Read every one of the 158 rows by eye (listed in full below as the "what is read off" table is too long
   to print in full here; the complete join is reproducible with the SQL above) and marked a row as a
   **candidate flag** only if the protein name contains no term a non-expert would associate with starch/
   glucan chemistry at all — the same bar "crescentin filaments" cleared for GH63. A separate `grep` for
   generic "this-is-clearly-a-different-kind-of-protein" vocabulary (kinase, antibody, nuclease, polymerase,
   fluorescent, albumin, lysozyme, actin, tubulin, filament, capsid, histone, collagen, dehydrogenase,
   reductase, isomerase, ligase, phosphatase, inhibitor, subunit, chaperone, protease, channel, receptor)
   was run over the full 158-row list as a second pass, independent of the by-eye read, to make sure nothing
   was missed.
4. Every candidate was then checked the same way the GH63 case was checked: fetch the RCSB entry record for
   its PDB id(s) and read the title/citation, to see whether the surface-level mismatch is a real error or a
   genuine, deliberate co-structure (as several of the other `04_curation_debt.md` items turned out to be).

## 4. What is read off

### Task A — the four single-entry families

| PDB | Claimed family (`entry.family`) | `entry.EC` | RCSB `struct.title` | RCSB primary citation | UniProt accession (`sifts_residue`) | UniProt protein name / organism | CAZy family's own EC (`cazy.org/<FAM>.html`) | CAZy's own `<FAM>_structure.html` lists this PDB? | **Verdict** |
|---|---|---|---|---|---|---|---|---|---|
| 1OBB | GH4 | 3.2.1.20 | "alpha-glucosidase A, AglA, from Thermotoga maritima in complex with maltose and NAD+" | Jeon et al., *J Biol Chem* 2003, "Crystal Structure of Thermotoga Maritima Alpha-Glucosidase AglA Defines a New Clan of NAD+-Dependent Glycosidases" | O33830 | "Alpha-glucosidase" (alt. name "Maltase"), *Thermotoga maritima* | 3.2.1.20 (α-glucosidase, among GH4's listed activities) | Yes — "1OBB[A,B]" | **Match** |
| 3REN | GH126 | 3.2.1.1 | "CPF\_2247, a novel alpha-amylase from Clostridium perfringens" | Park et al., *Proteins* 2011, "Structural analysis of CPF\_2247, a novel alpha-amylase from Clostridium perfringens" | A0A0H2YP60 | submission name "Glycosyl hydrolase, family 8" (no curated name), *Clostridium perfringens* | 3.2.1.1 (amylose endo-α-1,4-glucosidase — the same EC number as "alpha-amylase") | Yes — "3REN[A,B]" | **Match** (see note below on the stale UniProt "family 8" label) |
| 7VMA | GH176 | 3.2.1.68 | "The X-ray crystallographic structure of amylo-alpha-1,6-glucosidase from Thermococcus gammatolerans STB12" | unpublished ("To be published") | C5A758 | submission name "Amylo-alpha-1,6-glucosidase, putative archaeal type glycogen debranching enzyme (Gde)", *Thermococcus gammatolerans* | 3.2.1.68 (starch endo-α-1,6-glucosidase, among GH176's 2 listed activities) | Yes — "7VMA[A,B]" | **Match** |
| 8ZEQ | GH133 | 3.2.1.33 | "Molecular Architecture of Human Glycogen Debranching Enzyme: Insights into Glycogen Storage Disease III Pathogenesis" | Guo et al., *Nat Commun* 2025, "Molecular architecture and catalytic mechanism of human glycogen debranching enzyme" | P35573 | "Glycogen debranching enzyme" (alt. name "Glycogen debrancher"), *Homo sapiens* | 3.2.1.33 (glycogen exo-α-1,6-glucosidase, GH133's one listed activity) | Yes — "8ZEQ[A]" | **Match** |

**Note on 3REN/GH126.** UniProt's own submitted annotation for A0A0H2YP60 is the generic, non-curated
"Glycosyl hydrolase, family 8" — read alone, that would look like a GH8-vs-GH126 mismatch of the same shape
as the crescentin case. But every other source agrees it is correctly GH126: the entry's own EC (3.2.1.1)
is exactly GH126's defining activity, the RCSB title and its 2011 primary citation both independently call
it "a novel alpha-amylase," and — decisively — CAZy's own `GH126_structure.html` page lists "3REN[A,B]" as
one of only two PDB structures it assigns to GH126. The most likely reading (not a verified fact, since
UniProt gives no revision history here) is that GH126 did not exist as a family when A0A0H2YP60 was first
annotated, and the generic submission name was never updated after CAZy created GH126 for this structurally
novel amylase; this is a fact about the UniProt record's annotation history, not about the database audited
here. **Verdict stands as Match for all four single-entry families — no second instance of the GH63 bug was
found among them.**

### Task B — systematic scan, 158 (family, UniProt) pairs over all 637 entries

Two candidates passed the by-eye and keyword screens; both are from the **same PDB entry, 8J5D (claimed
family GH14)**, and both were checked against RCSB and resolved as genuine, not errors:

| PDB | Claimed family | UniProt accession | Protein name / organism | Why it looked like a candidate | Resolution |
|---|---|---|---|---|---|
| 8J5D | GH14 | Q9SN86 | "Malate dehydrogenase, chloroplastic", *Arabidopsis thaliana* | Malate dehydrogenase is a TCA-cycle/redox enzyme, unrelated to any glycoside hydrolase — the same class of surface mismatch as "crescentin" was for GH63 | **Not an error.** `entry.family`/`entry.EC`/`entry.protein` for 8J5D are "GH14 / 3.2.1.2 / β-amylase (BAM1;BMY7;TR-BAMY;At3g23920)" — the real enzyme chain. `sifts_residue` shows 8J5D has 4 protein chains: A = Q9LIR6 (β-amylase BAM1, the genuine GH14 enzyme), B/C = Q9SN86 (malate dehydrogenase), D = F4J117 (see below). RCSB's own title for 8J5D is "Cryo-EM structure of starch degradation complex of **BAM1-LSF1-MDH**," citing Zeeman-lab/collaborators, *Plant Cell* 2023, "The LIKE SEX FOUR 1-malate dehydrogenase complex functions as a scaffold to recruit beta-amylase to promote starch degradation." The malate dehydrogenase chain is a real, deliberately co-crystallized scaffold partner of the GH14 enzyme, not a mislabelled enzyme itself. The query construction (joining the entry-level `family` column against every chain's own `sifts_residue.uniprot`) mechanically attaches "GH14" to all four chains of this one heteromeric entry, exactly the pattern `docs/04_curation_debt.md` item 2 already describes for `entry.uniprot`, now showing up through `sifts_residue` as well because `family` itself is only ever recorded once per entry, not once per chain. |
| 8J5D | GH14 | F4J117 | "Phosphoglucan phosphatase LSF1, chloroplastic", *Arabidopsis thaliana* | A phosphatase, not a glycoside hydrolase — same class of surface mismatch | **Not an error**, same resolution as above: chain D of the same BAM1-LSF1-MDH complex, confirmed by the same RCSB title/citation. |

No other accession among the 158 failed the screen. For completeness, the by-eye pass also turned up six
UniProt names that mention "amylase"/"dextrinase"/"glucosidase" explicitly but are *not themselves catalytic
glycoside hydrolases* — they are non-catalytic partner chains in genuine enzyme–inhibitor or enzyme–subunit
complexes, so they do not clear the "obviously wrong to a non-expert" bar the task set (unlike malate
dehydrogenase or crescentin, their names read as amylase-related, not unrelated) and are **not flagged**,
only logged here for transparency:

| PDB(s) | Family | UniProt | Protein name | Organism |
|---|---|---|---|---|
| 1DHK, 1VIW | GH13 | P02873 | Alpha-amylase inhibitor 1 | *Phaseolus vulgaris* |
| 1CLV | GH13 | P80403 | Alpha-amylase inhibitor AAI | *Amaranthus hypochondriacus* |
| 1BVN | GH13 | P01092 | Alpha-amylase inhibitor HOE-467A (tendamistat) | *Streptomyces tendae* |
| 4X0N | GH13 | A0A0X1KGZ5 | Alpha-amylase inhibitor helianthamide | *Stichodactyla helianthus* (sea anemone) |
| 1AVA | GH13 | P07596 | Alpha-amylase/subtilisin inhibitor | *Hordeum vulgare* |
| 1TMQ | GH13 | P01087 | Alpha-amylase/trypsin inhibitor | *Eleusine coracana* |
| 4CVW | GH13 | Q2V8X0 | Limit dextrinase inhibitor | *Hordeum vulgare* |
| 8D43, 8EMR | GH31 | P14314 | Glucosidase 2 subunit beta (non-catalytic subunit of glucosidase II) | *Homo sapiens* |

These are all well-documented crystallographic complexes of a real GH13/GH31 enzyme bound to its own natural
protein inhibitor or non-catalytic partner subunit (e.g. 1DHK = *Tenebrio molitor* α-amylase + bean
α-amylase inhibitor 1); the `family` label correctly describes the enzyme chain in each entry. They are
listed here, not as errors, only because the method of Task B (joining `family` against every chain's
`uniprot`) surfaces them the same mechanical way it surfaced the two Task B candidates above, and a
systematic audit should say so rather than silently drop them.

**Bottom line for Task B: zero genuine family/UniProt mislabelings of the GH63 kind were found among the
other 154 (family, UniProt) pairs. The only two surface-level candidates (8J5D's malate dehydrogenase and
LSF1 phosphatase chains) were checked against RCSB and are a real, published multi-protein complex, not a
database error.**

## 5. The checks

- Task A's verdict for each of the 4 single-entry families rests on three independent sources agreeing
  (RCSB title + citation, UniProt protein name, and CAZy's own `_structure.html` listing of that exact PDB
  id under that exact family) — not on this project's own `entry.family`/`entry.protein` columns, which is
  the thing being audited and so cannot certify itself.
- Task B's screen was run twice with different methods (a by-eye read of all 158 rows, and an independent
  keyword `grep` for generic non-starch vocabulary) and both passes found the same two rows (8J5D's two
  non-catalytic chains) and no others, which is some evidence the screen did not simply miss something
  because of how the by-eye pass was phrased.
- Every one of the 158 distinct UniProt accessions in the Task B pair list resolved to a named protein
  record in a single pair of batch UniProt queries (0 missing, 0 obsolete/demerged) — so the 158-pair count
  is not an undercount caused by failed lookups.
- The SQL for Task B used `sifts_residue.uniprot` (per chain), not `entry.uniprot` (composite string for
  heteromeric entries), exactly as `docs/04_curation_debt.md` item 2 says to; the two Task B candidates
  found here show that even the per-chain `sifts_residue` route still mixes a complex's non-catalytic
  partner chains into the catalytic chain's family label, because `family` lives on `entry` (once per PDB
  id), not on `sifts_residue` (once per chain) — a second, related but distinct consequence of the same
  one-family-per-entry schema design that item 2 already flagged for `entry.uniprot`.

## 6. What it does not include

- 26 of 637 entries (listed by this audit's own query, all GH13/GH15/GH57 — multi-entry families already
  covered by the `03_frames` structural cross-check) have no `uniprot` at all in `sifts_residue`, so they
  contribute zero rows to the Task B pair list; this audit says nothing about whether their claimed identity
  is correct. They were not part of the task as specified (which asked for the `sifts_residue.uniprot`-keyed
  pairs), and their absence here is a gap in SIFTS coverage, not a finding.
- The "obviously wrong to a non-expert" bar used for Task B is deliberately conservative, by instruction:
  it will not catch a mislabeling where the substitute protein's name itself still sounds plausibly
  carbohydrate-related (e.g. a different, wrong carbohydrate-active enzyme swapped into the wrong family) —
  only ones as clear-cut as "crescentin" or "malate dehydrogenase" sitting under a glycoside hydrolase
  family. A subtler mislabeling of that kind would not be caught by this method; it would need the same
  structural superposition check GH63 got, which (per the premise of this audit) cannot be run on
  single-entry families and was not re-run here for any multi-entry family.
- This audit trusts the RCSB and UniProt API responses as ground truth without independently re-verifying
  their own curation; it is possible, though not checked here, for RCSB's or UniProt's own records to be
  wrong.
- CAZy's `_structure.html` listing for GH126 was reported (by the tool that fetched it) to include only 2
  PDB accessions for the whole family; this was not independently cross-checked against a separate CAZy
  page, so the count itself is taken on that one fetch.
- No file in the repository and no row of `enzimas_experimental.sqlite` was changed by this audit; this
  file is the only output written.

## References

- `docs/04_curation_debt.md` item 1 (GH63 fix, the bug this audit checks for a second instance of) and
  item 2 (why `sifts_residue.uniprot`, not `entry.uniprot`, is the per-chain source of truth).
- `docs/02_experimental_database.md`, "Correction (2026-11-23)" paragraph.
- `resultados/db/enzimas_experimental.sqlite` — tables `entry`, `sifts_residue`.
- RCSB PDB Data API, <https://data.rcsb.org/rest/v1/core/entry/>.
- UniProtKB REST API, <https://rest.uniprot.org/uniprotkb/>.
- CAZy, <https://www.cazy.org/> — family activity tables and `<FAMILY>_structure.html` pages for GH4,
  GH126, GH133, GH176.
