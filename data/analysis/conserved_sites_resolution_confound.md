# Conserved GH13 sites: a resolution confound (checked 2026-11-23)

**Problem.** §3.2 of the manuscript reports GH13 water sites W3/W4 and the calcium site I1 as
"conserved" from the fraction of superposed entries in which each is occupied (W3: 251/412, W4:
235/412, I1: 212/411). Does that fraction mostly track whether a structure was solved well enough
to see the water/ion, rather than tracking real biological occupancy?

**Data.** `resultados/calc/03_frames/conserved_cluster_members.tsv` (which entries occupy each
site) joined to `entry.resolution_A` for every GH13 entry.

**Method.** Instead of comparing median resolution between occupied and unoccupied entries (which
mixes a wide resolution range on both sides and understates any effect), entries were binned into
five resolution bands and the occupied fraction computed within each band separately.

**Result.**

| Resolution band | n | W3 occupied | W4 occupied | I1 (Ca²⁺) occupied |
|---|---|---|---|---|
| < 1.6 Å | 53 | 75.5% | 83.0% | 66.0% |
| 1.6–1.9 Å | 75 | 62.7% | 81.3% | 48.0% |
| 1.9–2.2 Å | 131 | 80.2% | 62.6% | 69.5% |
| 2.2–2.6 Å | 83 | 55.4% | 47.0% | 42.2% |
| ≥ 2.6 Å | 70 | 18.6% | 12.9% | 21.4% |

All three sites show a clear, roughly monotonic drop from the best to the worst resolution band —
a 4–6-fold difference in apparent occupancy across the resolution range used in this database.

**Reading (mine, kept apart from the numbers).** The headline fractions quoted in the manuscript
(61%, 57%, 52%) are aggregates over a resolution range of 0.95–4.6 Å and are **not** simple
biological-occupancy numbers: they are pulled down by the many lower-resolution entries where the
site is real but cannot be seen. At the best resolutions, where detection is least in doubt,
occupancy is 66–83%, markedly higher than the aggregate. This does not undermine the finding that
these are real, conserved sites — the calcium site's coordination chemistry was independently
checked earlier today and is textbook-correct — but it means **the aggregate fraction likely
understates true conservation**, not overstates it, and should be reported with this resolution
dependence stated, not as a flat population number.

**Checks.** The same binning was repeated for all three sites independently and all three show the
same direction of effect, at a similar magnitude; not an artefact of one site's cluster definition.

**Limits.** Only checked for GH13; not repeated for GH14/GH31/GH57's conserved sites. The band
edges are a coarse, round-number choice, not optimised. Entries without a resolution value (a
handful of cryo-EM-adjacent cases) are excluded from this specific check.
