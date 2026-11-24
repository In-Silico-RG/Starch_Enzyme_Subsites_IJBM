# Independence check of the subsite −1 grip (GH13, GH31)

Date: 2026-11-23. Read-only analysis; no manuscript file and no database was modified.
Scripts and intermediate tables: `resultados/analysis/grip_independence_check/` (step1.py … step5.py, `entries_minus1.tsv`,
`contact_matrix.tsv`, `ident.json`, `step4.out`, `step5.out`). Python: `~/miniconda3/bin/python` (Biopython 1.88, pandas, gemmi).

## 1. The problem

`manuscript_short/grip_finding.tex` (lines 108–137, Table `tab:grip`) states that at subsite −1 (the glucose ring whose
C1–O4 glycosidic bond is cleaved) six reference positions are contacted in 71–74 of 74 GH13 entries "from 39 distinct
proteins", and nine positions in 16–18 of 18 GH31 entries "from 6 distinct proteins". The argument is that a position
contacted in almost every entry is a fixed feature of the family. That reading needs the entries to be many independent
proteins. A protein deposited many times, or a group of near-identical orthologues, counts once in evolution but many
times in a per-entry fraction. The per-entry count cannot tell these apart, so the entries have to be grouped by
protein and by sequence similarity and the count repeated with one vote per group.

## 2. What the data contain

| file | content used | count |
|---|---|---|
| `resultados/analysis/33_subsites/units.tsv` | one row per sugar ring; columns `family, pdb, instance_id, ring, subsite, rscc, in_tables, frame_chain` | 2926 units |
| `resultados/analysis/33_subsites/contacts_mapped.tsv` | one row per ligand atom–residue contact (4.0 Å, from DB table `contact`); `rings, asu, ref_position, same_chain_as_frame` | 15 411 contacts |
| `resultados/analysis/33_subsites/subsite_summary.tsv` | per family/subsite: `n_entries`, `n_proteins` | GH13 −1: 74 entries, 39 "proteins"; GH31 −1: 18, 6 |
| `resultados/db/enzimas_experimental.sqlite` (read-only) | table `entry` (`uniprot, organism, protein`); table `sifts_residue` (`pdb, chain, uniprot`) | — |
| `resultados/census/cif/<PDB>.cif` | deposited chain sequence of the frame chain (entity full sequence, read with gemmi) | 92 entries |
| UniProt REST, `https://rest.uniprot.org/uniprotkb/<ACC>.json` and `.fasta` (fetched 2026-11-23) | sequence, recommended name, organism lineage | 47 accessions, 2 of them inactive |

**Reproduction of the entry counts.** `scripts/analysis/a33_subsites.py` counts an entry at subsite −1 when it has at least
one unit with `subsite == '-1'`, `rscc >= 0.80` and `in_tables == True` (sugar class, or 'other' class bonded into a set
that holds a sugar). A position is counted for that entry when a contact row of that unit has `asu == True`, the unit's ring
index in `rings`, a non-empty `ref_position`, and `same_chain_as_frame == True`. Rebuilt this way (`step3.py`,
output `contact_matrix.tsv`): GH13 74 entries, GH31 18 entries; per-position counts GH13 D289 74, H92 74, D179 73, R177 71,
Y51 71, H288 71; GH31 D327, H600, D443, D542, W539, W441 18 each, I328 17, F575 16, I364 16. These equal
`subsite_residues_reference_numbering.tsv` and the draft. No difference.

## 3. The method

1. **One protein per entry** (`step1.py`, `step2.py`). For each of the 92 entries the UniProt accession is taken from the
   SIFTS mapping of the frame chain (the chain whose superposition placed the −1 unit, `units.tsv: frame_chain`), because
   `entry.uniprot` is a `;`-joined list of every accession in the entry (e.g. `P04746 ; P04746 ; P19961 ; Q53F26` for
   human pancreatic amylase) or is empty. When SIFTS gives no single accession the `entry.uniprot` value is used if it is a
   single accession. When the accession is inactive in UniProt (A0A3P8MUS3 for 6GXV: demerged; A7M087 for 8DL2: deleted)
   or absent (1W9X, 24ZQ, 7P01, 7P07), the protein is represented by the deposited chain sequence in the mmCIF file and
   keyed `PDB:<id>`.
2. **Sequences.** UniProt full-length sequences for the accessions; deposited chain sequences for the six `PDB:` keys.
3. **Pairwise identity.** Biopython `PairwiseAligner`, global, BLOSUM62, gap open −10, gap extend −0.5, end gaps free
   (score 0). Two identities per pair:
   - $I_{\rm al} = N_{\rm id} / N_{\rm al}$, identical positions over aligned (non-gap) columns, the definition requested;
   - $I_{\rm sh} = N_{\rm id} / \min(L_1, L_2)$, identical positions over the length of the shorter sequence.
4. **Clusters.** Single linkage: two proteins share a cluster if a chain of pairs at or above the threshold joins them.
   Thresholds 95 %, 90 %, 50 %, 30 % (40 % added as a check).
5. **Cluster-weighted counting.** For each cluster $c$ and position $q$, $f_{c,q}$ = fraction of the cluster's entries in
   which $q$ contacts the −1 unit. Reported: number of clusters with $f_{c,q} > 0$ ("any"), with $f_{c,q} = 1$ ("all"), and the
   cluster-weighted fraction $\bar f_q = \frac{1}{C}\sum_c f_{c,q}$, to be compared with the script's conserved threshold 0.80.
6. **Taxonomy.** UniProt `organism.lineage` for each accession; for the six `PDB:` keys the lineage of the organism named
   in `entry.organism`, taken from the UniProt taxonomy service.

**The one assumption.** Proteins in different sequence clusters count as independent observations of the site; proteins
in one cluster count as one. Sequence identity stands in for evolutionary distance.

## 4. What is read off

### 4.1 Protein counts

| quantity | GH13 | GH31 | source |
|---|---|---|---|
| entries at −1 (RSCC ≥ 0.80) | 74 | 18 | `units.tsv`, reproduced |
| distinct `entry.uniprot` strings (what the draft counted) | 39 | 6 | `subsite_summary.tsv: n_proteins` |
| protein keys (accession of frame chain, or `PDB:` chain) | 45 | 6 | `entries_minus1.tsv: protein_key` |
| distinct proteins after merging identical sequences (≥ 99.8 %) | **42** | **6** | `ident.json` |
| clusters at 95 % ($I_{\rm al}$) | 41 | 6 | `step4.out` |
| clusters at 90 % ($I_{\rm al}$) | 40 | 6 | `step4.out` |
| clusters at 50 % ($I_{\rm al}$ / $I_{\rm sh}$) | **29** / 31 | **6** / 6 | `step4.out` |
| clusters at 40 % ($I_{\rm al}$ / $I_{\rm sh}$) | 18 / 20 | 4 / 5 | `step5.out` |
| clusters at 30 % ($I_{\rm al}$ / $I_{\rm sh}$) | 1 / **5** | 1 / **4** | `step4.out`, `step5.out` |

**Why 39 is wrong for GH13.** `a33_subsites.py` computes `n_proteins` as the number of distinct `entry.uniprot` strings.
Eight entries have an empty string and were counted as one protein: 1W9X (*Sutcliffiella halmapala*), 24ZQ (*Alkalimonas*),
26WD (*Apis mellifera*, A0A7M6W876), 2GJP (P19571, also present through 1WPC), 7P01 and 7P07 (*Blastobotrys
adeninivorans*, identical), 9UNS and 9XG6 (*Klebsiella grimontii*, A0ABD7ANA7 / A0A7H4P971, 99.8 % identical). These eight
are five proteins not counted elsewhere. In the other direction two strings (`A0A1S9DH83 ; …` for 6YQB and
`A2QL05 ; A2QL05 ; P0C1B3` for 2GVY) are the same Taka-amylase sequence (A0A1S9DH83 vs P0C1B3: 100 % over 499 columns).
38 non-empty strings − 1 duplicate + 5 = **42 distinct proteins**. GH31: 6 is right; SIFTS gives one accession per frame
chain and the six are distinct.

### 4.2 How the entries are spread over proteins

GH13 (`entries_minus1.tsv`): 1 protein with 11 entries (human pancreatic α-amylase P04746), 1 with 7 (human salivary
α-amylase P0DUB6), 2 with 3, 9 with 2, 32 with 1 (counts over the 45 keys). The two human amylases, 96.9 % identical,
supply 18 of 74 entries (24 %).
GH31: 5 (sugar-beet α-glucosidase L0N7E5), 4 (human lysosomal α-glucosidase P10253), 4 (*Gracilariopsis* α-1,4-glucan
lyase Q9STC1), 3 (*Lactococcus* A2RM80), 1, 1.

Protein table, GH13 (grp95 = cluster at 95 % $I_{\rm al}$; cl50 = cluster at 50 % $I_{\rm al}$; cl30s = cluster at 30 % $I_{\rm sh}$).
Note: 2GVY is listed in `entry.organism` as *Aspergillus niger*; its frame chain maps in SIFTS to *A. oryzae* P0C1B3.

| # | protein key | organism | protein (UniProt name) | entries | PDB entries | grp95 | cl50 | cl30s | sequence source |
|---|---|---|---|---|---|---|---|---|---|
| 1 | A0A0C5GWS2 | *Paenibacillus barengoltzii* | Pulullanase | 2 | 6JFJ 6JFX | 1 | 1 | 1 | UniProt A0A0C5GWS2 |
| 2 | A0A0H2UNG0 | *Streptococcus pneumoniae serotype 4* | Pullulanase A | 1 | 2YA2 | 2 | 2 | 1 | UniProt A0A0H2UNG0 |
| 3 | P0C1B3 | *Aspergillus oryzae* | Alpha-amylase A type-1/2 | 2 | 2GVY 7TAA | 3 | 3 | 2 | UniProt P0C1B3 |
| 4 | A0A1S9DH83 | *Aspergillus oryzae* | Alpha-amylase | 1 | 6YQB | 3 | 3 | 2 | UniProt A0A1S9DH83 |
| 5 | A0A5S8WFA7 | *Thamnidium elegans* | alpha-amylase | 1 | 6SAO | 4 | 4 | 2 | UniProt A0A5S8WFA7 |
| 6 | A0A5S8WFA8 | *Cordyceps farinosa* | Alpha-amylase | 1 | 6SAU | 5 | 5 | 2 | UniProt A0A5S8WFA8 |
| 7 | A0A7H4P971 | *Klebsiella grimontii* | Maltodextrin glucosidase | 1 | 9XG6 | 6 | 6 | 2 | UniProt A0A7H4P971 |
| 8 | A0ABD7ANA7 | *Klebsiella grimontii* | Maltodextrin glucosidase | 1 | 9UNS | 6 | 6 | 2 | UniProt A0ABD7ANA7 |
| 9 | P04746 | *Homo sapiens* | Pancreatic alpha-amylase | 11 | 1XCW 1XCX 1XD0 3BAJ 3IJ8 3IJ9 3OLD 3OLE 3OLG 3OLI 5TD4 | 21 | 7 | 3 | UniProt P04746 |
| 10 | P0DUB6 | *Homo sapiens* | Alpha-amylase 1A | 7 | 1MFU 1MFV 1NM9 1Z32 3BLK 3BLP 3DHP | 21 | 7 | 3 | UniProt P0DUB6 |
| 11 | P29957 | *Pseudoalteromonas haloplanktis* | Alpha-amylase | 3 | 1G94 1G9H 8CQF | 28 | 7 | 3 | UniProt P29957 |
| 12 | P00690 | *Sus scrofa* | Pancreatic alpha-amylase | 2 | 1WO2 3L2L | 17 | 7 | 3 | UniProt P00690 |
| 13 | A0A7M6W876 | *Apis mellifera* | Alpha-amylase | 1 | 26WD | 7 | 7 | 3 | UniProt A0A7M6W876 |
| 14 | P08144 | *Drosophila melanogaster* | Alpha-amylase A | 1 | 8ORP | 23 | 7 | 3 | UniProt P08144 |
| 15 | I1VWH9 | *Anoxybacillus ayderensis* | alpha-amylase | 2 | 5A2B 5A2C | 13 | 8 | 2 | UniProt I1VWH9 |
| 16 | A0ABF7PFV7 | *Geobacillus thermoleovorans CCB_US3_UF5* | alpha-amylase | 1 | 4E2O | 8 | 8 | 2 | UniProt A0ABF7PFV7 |
| 17 | B8CZ54 | *Halothermothrix orenii* | Alpha amylase | 2 | 3BC9 3BCD | 9 | 9 | 2 | UniProt B8CZ54 |
| 18 | PDB:24ZQ | *Alkalimonas sp. NCh-2* | α-amylase (AmyY) | 1 | 24ZQ | 32 | 9 | 2 | chain A of 24ZQ mmCIF (no UniProt mapping) |
| 19 | D2YYD7 | *Arthrobacter globiformis* | Alpha-glucosidase | 1 | 8YIE | 10 | 10 | 2 | UniProt D2YYD7 |
| 20 | D2YYE1 | *Arthrobacter globiformis* | Cyclic maltosyl-maltose hydrolase | 2 | 6A0K 6A0L | 11 | 11 | 2 | UniProt D2YYE1 |
| 21 | H3K096 | *Halomonas sp. H11* | Alpha-glucosidase | 1 | 3WY2 | 12 | 12 | 2 | UniProt H3K096 |
| 22 | P19571 | *Bacillus sp.* | Glucan 1,4-alpha-maltohexaosidase | 2 | 1WPC 2GJP | 26 | 13 | 2 | UniProt P19571 |
| 23 | I3P686 | *Bacillus licheniformis* | Amylase | 1 | 6TOZ | 14 | 13 | 2 | UniProt I3P686 |
| 24 | PDB:1W9X | *Sutcliffiella halmapala* | α-amylase | 1 | 1W9X | 31 | 13 | 2 | chain A of 1W9X mmCIF (no UniProt mapping) |
| 25 | PDB:6GXV | *Alicyclobacillus sp. 18711* | α-amylase (AliC) | 1 | 6GXV | 33 | 13 | 2 | chain A of 6GXV mmCIF (UniProt A0A3P8MUS3 inactive) |
| 26 | K9L0H1 | *Anoxybacillus sp. LM18-11* | Type I pullulanase | 2 | 3WDI 3WDJ | 15 | 14 | 1 | UniProt K9L0H1 |
| 27 | O06994 | *Bacillus subtilis* | Oligo-1,6-glucosidase 1 | 1 | 5WCZ | 16 | 15 | 2 | UniProt O06994 |
| 28 | Q8YAE6 | *Listeria monocytogenes serovar 1/2a* | oligo-1,6-glucosidase | 1 | 5DO8 | 40 | 15 | 2 | UniProt Q8YAE6 |
| 29 | P00691 | *Bacillus subtilis* | Alpha-amylase | 1 | 1UA7 | 18 | 16 | 4 | UniProt P00691 |
| 30 | P00693 | *Hordeum vulgare* | Alpha-amylase type A isozyme | 1 | 1RPK | 19 | 17 | 2 | UniProt P00693 |
| 31 | P04063 | *Hordeum vulgare* | Alpha-amylase type B isozyme | 1 | 1BG9 | 20 | 17 | 2 | UniProt P04063 |
| 32 | P07206 | *Klebsiella pneumoniae* | Pullulanase | 3 | 2FHB 2FHC 2FHF | 22 | 18 | 1 | UniProt P07206 |
| 33 | P13507 | *Stutzerimonas stutzeri* | Glucan 1,4-alpha-maltotetraohydrolase | 1 | 1QI5 | 24 | 19 | 2 | UniProt P13507 |
| 34 | P22963 | *Roseateles saccharophilus* | Glucan 1,4-alpha-maltotetraohydrolase | 1 | 6JQB | 27 | 19 | 2 | UniProt P22963 |
| 35 | P19531 | *Geobacillus stearothermophilus* | Maltogenic alpha-amylase | 2 | 1QHO 1QHP | 25 | 20 | 2 | UniProt P19531 |
| 36 | P38940 | *Geobacillus stearothermophilus* | Neopullulanase | 1 | 1J0I | 29 | 21 | 2 | UniProt P38940 |
| 37 | P53051 | *Saccharomyces cerevisiae* | Oligo-1,6-glucosidase IMA1 | 1 | 3A4A | 30 | 22 | 2 | UniProt P53051 |
| 38 | PDB:7P01 | *Blastobotrys adeninivorans CBS 8244* | maltase (BaAG2) | 1 | 7P01 | 34 | 23 | 2 | chain A of 7P01 mmCIF (no UniProt mapping) |
| 39 | PDB:7P07 | *Blastobotrys adeninivorans CBS 8244* | maltase (BaAG2) | 1 | 7P07 | 34 | 23 | 2 | chain A of 7P07 mmCIF (no UniProt mapping) |
| 40 | PDB:8DL2 | *Bacteroides ovatus ATCC 8483* | α-amylase / neopullulanase (BoSusG; BoGH13A; BoGH13ASus; Bovatus_03803; BACOVA_03514) | 1 | 8DL2 | 35 | 24 | 5 | chain A of 8DL2 mmCIF (UniProt A7M087 inactive) |
| 41 | Q08751 | *Thermoactinomyces vulgaris* | Neopullulanase 2 | 1 | 3A6O | 36 | 25 | 2 | UniProt Q08751 |
| 42 | Q60053 | *Thermoactinomyces vulgaris* | Neopullulanase 1 | 1 | 1UH3 | 37 | 26 | 2 | UniProt Q60053 |
| 43 | Q7LYT7 | *Pyrococcus woesei* | Alpha amylase | 1 | 1MXD | 38 | 27 | 2 | UniProt Q7LYT7 |
| 44 | Q8A1G3 | *Bacteroides thetaiotaomicron* | Alpha-amylase SusG | 1 | 3K8M | 39 | 28 | 2 | UniProt Q8A1G3 |
| 45 | Q9FYY0 | *Hordeum vulgare* | Limit dextrinase | 1 | 4J3V | 41 | 29 | 1 | UniProt Q9FYY0 |


Protein table, GH31:

| # | protein key | organism | protein (UniProt name) | entries | PDB entries | grp95 | cl50 | cl30s | sequence source |
|---|---|---|---|---|---|---|---|---|---|
| 1 | A2RM80 | *Lactococcus lactis subsp. cremoris* | Alpha-xylosidase | 3 | 7WJB 7WJC 7WJF | 1 | 1 | 1 | UniProt A2RM80 |
| 2 | A5ZY13 | *Blautia obeum ATCC 29174* | Glycosyl hydrolase, family 31 | 1 | 3PHA | 2 | 2 | 2 | UniProt A5ZY13 |
| 3 | L0N7E5 | *Beta vulgaris* | Alpha-glucosidase | 5 | 3W37 3WEL 3WEM 3WEN 3WEO | 3 | 3 | 3 | UniProt L0N7E5 |
| 4 | O43451 | *Homo sapiens* | Maltase-glucoamylase | 1 | 3TOP | 4 | 4 | 3 | UniProt O43451 |
| 5 | P10253 | *Homo sapiens* | Lysosomal alpha-glucosidase | 4 | 5KZW 5KZX 5NN5 5NN8 | 5 | 5 | 3 | UniProt P10253 |
| 6 | Q9STC1 | *Gracilariopsis lemaneiformis* | Alpha-1,4-glucan lyase isozyme 1 | 4 | 2X2I 2X2J 4AMW 4AMX | 6 | 6 | 4 | UniProt Q9STC1 |

### 4.3 Pairwise identity: the close pairs

Every GH13 pair with $I_{\rm al} \ge 0.50$ (from `ident.json`; $I_{\rm al}$, aligned columns):

| pair | $I_{\rm al}$ | columns |
|---|---|---|
| A0A1S9DH83 – P0C1B3 (Taka-amylase, two accessions) | 1.000 | 499 |
| PDB:7P01 – PDB:7P07 (*B. adeninivorans*) | 1.000 | 585 |
| A0A7H4P971 – A0ABD7ANA7 (*K. grimontii*) | 0.998 | 605 |
| P04746 – P0DUB6 (human pancreatic – salivary) | 0.969 | 511 |
| P13507 – P22963 (maltotetraose-forming amylases) | 0.909 | 548 |
| P00690 – P04746 (pig – human pancreatic) | 0.863 | 511 |
| P19571 – PDB:1W9X | 0.859 | 481 |
| P00690 – P0DUB6 | 0.857 | 511 |
| P00693 – P04063 (barley AMY1 – AMY2) | 0.742 | 426 |
| PDB:1W9X – I3P686 | 0.715 | 478 |
| A0ABF7PFV7 – I1VWH9 | 0.700 | 504 |
| I3P686 – P19571 / PDB:6GXV / 1W9X–6GXV / P19571–6GXV | 0.637–0.692 | 478–483 |
| A0A7M6W876 – P08144 (honeybee – *Drosophila*) | 0.615 | 488 |
| animal-type amylases (human ×2, pig, *Drosophila*, honeybee) among themselves | 0.512–0.566 | 482–489 |
| P00690 – P29957 (pig – *Pseudoalteromonas*) | 0.511 | 466 |
| O06994 – Q8YAE6 (*B. subtilis* – *Listeria* oligo-1,6-glucosidases) | 0.555 | 546 |
| B8CZ54 – PDB:24ZQ | 0.503 | 477 |

All other GH13 pairs: $I_{\rm al}$ 0.235–0.493, median 0.308. GH31: all 15 pairs 0.281–0.453 (largest: MGAM O43451 – GAA
P10253, 0.453 over 912 columns).

The 50 % clusters of GH13 (29): one cluster holds 6 proteins and **25 of 74 entries** — the animal-type α-amylases
(human pancreatic, human salivary, pig pancreatic, *Drosophila*, honeybee) together with *Pseudoalteromonas haloplanktis*
AHA. The next largest hold 5 entries (four *Bacillus*-type liquefying/maltohexaose amylases), then 3, 3, 3, 3; 9 clusters
of 2 entries and 14 of 1 (`step5.out`).

The 30 % clusters of GH13 by $I_{\rm sh}$ (5): (1) the pullulanases and limit dextrinase — 5 proteins, 9 entries;
(2) 32 proteins, 38 entries (all remaining α-amylases, neopullulanases, oligo-1,6-glucosidases, maltases, the fungal and
plant amylases); (3) the animal-type cluster above — 6 proteins, 25 entries; (4) *B. subtilis* α-amylase P00691 (1UA7);
(5) *B. ovatus* SusG-like 8DL2. GH31 (4): (1) A2RM80; (2) A5ZY13; (3) sugar-beet α-glucosidase + human MGAM + human GAA —
3 proteins, 10 entries; (4) *Gracilariopsis* lyase.

### 4.4 Taxonomic spread (`step5.py`, UniProt lineage; counts over protein keys / entries)

| family | domain / kingdom-level group / phylum | proteins | entries |
|---|---|---|---|
| GH13 | Bacteria / Bacillota | 17 | 23 |
| GH13 | Bacteria / Pseudomonadota | 8 | 12 |
| GH13 | Bacteria / Actinomycetota | 2 | 3 |
| GH13 | Bacteria / Bacteroidota | 2 | 2 |
| GH13 | Archaea / Methanobacteriota (UniProt lineage of *Pyrococcus woesei*, class Thermococci) | 1 | 1 |
| GH13 | Eukaryota / Fungi / Ascomycota | 6 | 7 |
| GH13 | Eukaryota / Fungi / Mucoromycota | 1 | 1 |
| GH13 | Eukaryota / Metazoa / Chordata (human, pig) | 3 | 20 |
| GH13 | Eukaryota / Metazoa / Arthropoda | 2 | 2 |
| GH13 | Eukaryota / Viridiplantae / Streptophyta (barley) | 3 | 3 |
| GH31 | Bacteria / Bacillota | 2 | 4 |
| GH31 | Eukaryota / Metazoa / Chordata (human) | 2 | 5 |
| GH31 | Eukaryota / Viridiplantae / Streptophyta (sugar beet) | 1 | 5 |
| GH31 | Eukaryota / Rhodophyta (red alga) | 1 | 4 |

GH13: all three domains; bacteria 4 phyla (40 entries), archaea 1 (1 entry), eukaryotes 3 kingdoms — fungi, animals,
plants — 5 phyla (33 entries). GH31: bacteria (1 phylum) and eukaryotes (animals, land plant, red alga); no archaea, no fungi.

### 4.5 Key check: votes per cluster instead of per entry

GH13, six positions (reference 1AMY numbering). Per entry from `contact_matrix.tsv`; per cluster from `step4.out` /
`step5.out`.

| position | per entry (of 74) | 50 % $I_{\rm al}$: any / all (of 29) | $\bar f_q$ (50 %) | 30 % $I_{\rm sh}$: any / all (of 5) | $\bar f_q$ (30 %) |
|---|---|---|---|---|---|
| Asp289 | 74 | 29 / 29 | 1.000 | 5 / 5 | 1.000 |
| His92 | 74 | 29 / 29 | 1.000 | 5 / 5 | 1.000 |
| Asp179 | 73 | 29 / 28 | 0.983 | 5 / 4 | 0.995 |
| Arg177 | 71 | 28 / 27 | 0.964 | 5 / 3 | 0.981 |
| Tyr51 | 71 | 29 / 27 | 0.980 | 5 / 3 | 0.979 |
| His288 | 71 | 28 / 27 | 0.964 | 5 / 3 | 0.981 |

Entries without the contact: Asp179 — 1QHP; Arg177 and His288 — 3L2L (pig), 9UNS and 9XG6 (*K. grimontii* maltodextrin
glucosidase, both entries of that 50 % cluster); Tyr51 — 1XD0, 3BAJ (human pancreatic), 6JQB.

GH31, nine positions (reference 2QMJ numbering):

| position | per entry (of 18) | 50 % $I_{\rm al}$: any / all (of 6) | $\bar f_q$ (50 %) | 30 % $I_{\rm sh}$: any / all (of 4) | $\bar f_q$ (30 %) |
|---|---|---|---|---|---|
| Asp327 | 18 | 6 / 6 | 1.000 | 4 / 4 | 1.000 |
| His600 | 18 | 6 / 6 | 1.000 | 4 / 4 | 1.000 |
| Asp443 | 18 | 6 / 6 | 1.000 | 4 / 4 | 1.000 |
| Asp542 | 18 | 6 / 6 | 1.000 | 4 / 4 | 1.000 |
| Trp539 | 18 | 6 / 6 | 1.000 | 4 / 4 | 1.000 |
| Trp441 | 18 | 6 / 6 | 1.000 | 4 / 4 | 1.000 |
| Ile328 | 17 | 6 / 5 | 0.958 | 4 / 3 | 0.938 |
| Phe575 | 16 | 5 / 4 | **0.800** | 4 / 3 | 0.950 |
| Ile364 | 16 | 6 / 5 | 0.933 | 4 / 3 | 0.950 |

Entries without the contact: Ile328 — 4AMW; Phe575 — 3TOP (human MGAM, the only entry of its protein) and 3WEN;
Ile364 — 3WEM, 3WEN.

**Facts.** With one vote per 50 % cluster, every GH13 position keeps $\bar f_q \ge 0.964$ and is contacted in at least
28 of 29 clusters; with one vote per 30 % ($I_{\rm sh}$) cluster, every position is contacted in all 5 clusters. Every
GH31 position keeps $\bar f_q \ge 0.80$ at 50 % and $\ge 0.938$ at 30 %; Phe575 at 50 % sits exactly on the 0.80 threshold.

**Interpretation.** The grip does not rest on the redundancy of the set: weighting by cluster changes the GH13 fractions by
at most 0.02 and none falls below 0.80. The number of independent groups behind the GH13 claim is 29 at 50 % identity,
not 39 or 42 proteins, and about a third of the entries come from one animal-type amylase cluster; the claim should be
worded in clusters (or proteins with that caveat). For GH31 the evidence is six proteins that are all below 46 % identical,
from four lineages at 30 %; Phe575 is the one position whose conserved status depends on how the votes are counted
(0.800 per protein, 0.889 per entry).

## 5. The checks

| check | predicted | measured | source |
|---|---|---|---|
| reproduce entries at −1 | 74 / 18 | 74 / 18 | `contact_matrix.tsv` |
| reproduce per-position entry counts | draft values | identical for all 15 positions | `step3.py` output |
| draft "39" recomputed as in the script | 39 | 39 (distinct `entry.uniprot` strings, '' counted once) | `step1.py` |
| UniProt clustering vs deposited chain sequences (one sequence per entry, same aligner) | same cluster number | GH13 41 / 40 / 29 at 95 / 90 / 50 %; GH31 6 / 6 / 6 — same as UniProt | `step4.out` ("entry-level chain-seq clusters") |
| background identity of unrelated sequences (12 pairs, one partner shuffled) | — | $I_{\rm al}$ mean 0.293 (0.249–0.331); $I_{\rm sh}$ mean 0.182 (0.084–0.249) | `step5.out` |

The background check decides the 30 % threshold. $I_{\rm al}$ of a sequence aligned to a shuffled sequence is 0.25–0.33,
so a 30 % cut on $I_{\rm al}$ sits inside the noise and single linkage then joins everything (1 cluster in each family).
$I_{\rm sh}$ has a background of at most 0.25, so the 30 % clusters are read from $I_{\rm sh}$ (5 for GH13, 4 for GH31).
At 50 % the two definitions agree (29 vs 31 for GH13; 6 vs 6 for GH31).

## 6. What it does not include

- Sequence identity is a proxy for evolutionary independence. All GH13 members are homologous by family definition, as
  are all GH31 members; no two clusters are phylogenetically independent in the strict sense. The count says how many
  divergent lineages carry the contact, not how many independent origins it had. No phylogeny was built.
- Single linkage can chain clusters together; the 30 % $I_{\rm sh}$ cluster 2 of GH13 (32 proteins) is such a chain.
- Full-length UniProt sequences include signal peptides and extra domains (pullulanases, MGAM 2753 aa); the deposited
  chain sequences gave the same cluster counts at 95, 90 and 50 %.
- Inactive UniProt accessions (A0A3P8MUS3, A7M087) and unmapped entries (1W9X, 24ZQ, 7P01, 7P07) were represented by
  their deposited chain sequences; their organism comes from `entry.organism`.
- The contact definition, RSCC filter, reference-position mapping (2.0 Å Cα) and subsite assignment are those of
  `a33_subsites.py`, taken as given; they were reproduced, not re-validated.
- Only subsite −1 and only GH13 and GH31 were checked. The draft sentence "In GH13 (412 structures, 39 distinct proteins)"
  reads as if 39 referred to all 412 structures; the 39 came from the 74 entries at −1. The protein count over all 412
  GH13 structures was not computed here.
- Taxonomy uses the UniProt lineage strings current on 2026-11-23 (e.g. *Pyrococcus* placed under Methanobacteriota).
