# Structure–activity knowledge gap (exploratory check, 2026-11-23)

**Problem.** Of everything that binds these enzymes (deposited structures) and everything that has been
measured against them (ChEMBL), how much is known by both routes, and how much by only one?

**Data.** `resultados/analysis/31_chemspace/species.tsv`: 112 deposited species, 940 ChEMBL molecules,
clustered by chemical similarity (Butina, Tanimoto distance 0.4) into 323 clusters.

**Result.**

| Cluster type | Count |
|---|---|
| Deposited only (never measured in ChEMBL) | 52 clusters |
| Measured only (never co-crystallised) | 257 clusters |
| Both | 14 clusters |

Of the 112 deposited species, only 37 (33%) fall in a cluster with a ChEMBL-measured molecule. Of the 940
ChEMBL molecules, only 123 (13%) fall in a cluster with a deposited ligand.

By family (deposited species with no chemically similar ChEMBL-measured molecule):

| Family | Deposited species | With no ChEMBL neighbour |
|---|---|---|
| GH13 | 56 | 38 |
| GH31 | 27 | 12 |
| GH14 | 9 | 8 |
| GT35 | 4 | 4 |

**Reading (mine).** Most of what is in the crystals has never had its potency measured, and most of what
has been measured has never been seen bound in a structure. GH13 and GH31 — the families with the human
drug targets — have the largest absolute gaps.

**Limits.** Cluster membership depends on the 0.4 Tanimoto cutoff; a different cutoff would move some
borderline cases between "gap" and "shared". Being in the same cluster as a measured molecule does not mean
the exact deposited species itself was measured.
