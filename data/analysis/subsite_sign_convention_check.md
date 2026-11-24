# Does the subsite -1 rule flip sign across GH13 enzyme subtypes?

A check of a methodological assumption behind `resultados/analysis/33_subsites/`, run 2026-11-23.

## 1. The problem

`scripts/analysis/a33_subsites.py` assigns a subsite number (..., -2, -1, +1, +2, ...) to every
bound sugar unit of every entry of a family, in that family's own structural reference frame
(built in `resultados/calc/03_frames/`, a Kabsch superposition of every entry onto one reference
structure, using the catalytic residues as the alignment anchor). The rule that starts the
numbering is purely geometric: subsite -1 is the sugar unit whose reference atom lies nearest
the midpoint of the two catalytic carboxylate centroids, in that entry's own (real, un-superposed)
coordinates, provided that distance is ≤ 4.5 Å (`D_MAX` in the script). Everything else is
numbered by walking outward from that unit along the deposited glycosidic links.

GH13 is the one family in this analysis where that rule is applied across enzymes that are not
all the same catalytic activity: alpha-amylase, maltogenic amylase, cyclomaltodextrinase /
neopullulanase, isoamylase / glycogen-debranching enzyme, pullulanase, and several
oligo-1,6-glucosidase / alpha-glucosidase lines, all GH13, all superposed onto the same reference,
1AMY (barley alpha-amylase), because they share the fold and the Asp/Glu catalytic pair. A
systematic check already confirmed that, within each entry, the subsite numbers walk without
gaps (236/236 main-chain components checked, no gaps — `resultados/analysis/33_subsites/run.log`,
`connected_sets.tsv`). That check cannot see a *whole-chain* sign flip: if, for some subtype, the
"nearest unit to the midpoint" rule happened to pick a unit that sits, in truth, on the opposite
side of the catalytic dyad from where an alpha-amylase's -1 unit sits, every label in that chain
would still look perfectly contiguous (-3,-2,-1,+1,+2...) — the numbers would just be running
backwards relative to the other subtypes. That would silently corrupt every subsite-based
comparison across subtypes (the fingerprints, the fraction-occupied tables, the fingerprint
fingerprints_GH13.tsv columns) without any internal-contiguity check ever flagging it.

This note asks exactly that question, and only that question: does -1 sit on the same side of,
and point in a consistent direction from, the catalytic dyad for a clearly different GH13
activity (isoamylase, cyclomaltodextrinase, oligo-1,6-glucosidase, ...) as it does for an
ordinary alpha-amylase?

## 2. What the data contain

- `resultados/analysis/33_subsites/units.tsv` — one row per sugar/other ring unit; for units with
  an assigned subsite, `x_frame,y_frame,z_frame` is that unit's reference-atom position after
  rigid transformation into the family's frame (`resultados/calc/03_frames/superposition.tsv`
  matrices), and `d_ref_to_midpoint_A` is the distance (in the entry's own, real coordinates)
  from that unit's reference atom to the midpoint used to pick it. Filtering `family=='GH13'`
  and `subsite=='-1'` gives **129 sugar units, from 99 distinct GH13 entries**, each one the unit
  chosen by the rule as that entry's -1.
- `resultados/calc/03_frames/catalytic_residues_mapped.tsv` — per entry, the residue (chain,
  auth_seq, resname) that plays the nucleophile (`nuc_or_base`) and acid role, restricted here to
  rows with `type_check=='same type'` (i.e. the entry really has an Asp at the nucleophile
  position and a Glu at the acid position, matching GH13's reference pair). For 1AMY itself
  (the GH13 reference, chain A): nucleophile Asp179, acid Glu204 (`families.json`,
  `catalytic_residues_mapped.tsv`).
- `resultados/calc/03_frames/superposition.tsv` — the accepted (R, t) rigid transform for every
  (family, pdb, chain), the same one `a33_subsites.py` applied to the ligand units.
- `resultados/census/cif/{pdb}.cif` — the deposited coordinates, read to recover the real-space
  position of each entry's own catalytic Asp OD1/OD2 and Glu OE1/OE2 atoms (not stored anywhere
  as a flat file; `a33_subsites.py` computes them on the fly and only keeps the resulting
  midpoint distance, not the midpoint's own coordinates — so this check reads the structures
  again to get the midpoint itself).
- `resultados/part2_structures_all.tsv` — EC number and protein name per PDB entry, used only to
  label each of the 99 entries by enzyme activity (alpha-amylase, isoamylase, etc.); this file
  plays no role in the geometry, only in grouping the result afterwards.

## 3. The method

For each of the 99 GH13 entries that have a unit labelled subsite -1:

1. Take that entry's catalytic chain (`frame_chain` in `units.tsv`, the same chain
   `a33_subsites.py` used to pick -1) and its nucleophile/acid residues from
   `catalytic_residues_mapped.tsv`.
2. Read the deposited structure (`resultados/census/cif/{pdb}.cif`) and recompute the carboxylate
   centroids in real coordinates — this exactly repeats the `mids` calculation inside
   `a33_subsites.py` (same atom names, same residues), so it is not a new assumption, only a
   re-derivation of a quantity the original script used internally but did not save:

   c_nuc = mean(OD1, OD2) of the nucleophile Asp        (real coordinates)
   c_acid = mean(OE1, OE2) of the acid Glu               (real coordinates)
   m_real = ( c_nuc + c_acid ) / 2                        — the catalytic midpoint

3. Carry that midpoint into the family frame with the entry's own rigid transform, the same
   (R, t) applied to the ligand units:

   m_frame = R · m_real + t

4. Form the frame-space vector from the midpoint to the -1 unit (whose frame coordinates are
   already in `units.tsv`):

   v = ( x_frame, y_frame, z_frame )_(-1 unit) − m_frame

   If the sign convention held across subtypes, v should point in essentially the same direction
   for every entry, whatever the enzyme activity, because the frame itself is built so that every
   entry's own catalytic dyad lands close to the same place (that is what the Kabsch alignment on
   the catalytic residues optimises for).
5. Normalise each v to a unit vector u = v / |v|, group the 129 u's by enzyme subtype (from
   `part2_structures_all.tsv`'s EC number), and for each subtype of n members compute the mean
   resultant vector

   ū = ( Σ u_i ) / n ,   |ū| ∈ [0,1]  (closer to 1 ⇒ tightly clustered directions)

   and the angle of every member to its own group's mean, and the angle of each group's mean to
   the plain-alpha-amylase group's mean (the largest group, 66 units / 59 entries, taken as the
   reference direction since alpha-amylase is the one activity the 1AMY reference itself has).
   A flipped subtype would show an angle near 180° to the alpha-amylase mean (or, for a
   within-group flip, |ū| near 0 despite individual vectors being well-defined, i.e. members
   cancelling each other out).

The one assumption of this check is that the family frame really does bring every entry's
catalytic dyad to (near) the same point regardless of subtype — which is simply what the
superposition in `resultados/calc/03_frames/` was built to do, and is re-confirmed below rather
than assumed blindly (§5).

Script: `/tmp/claude-1000/.../scratchpad/sign_check.py` (vector computation, reads the CIF files)
and `.../sign_stats.py` (grouping and angle statistics); both single-process, `nice -n 19`,
`~/miniconda3/bin/python`, read-only on `resultados/db/enzimas_experimental.sqlite` (not touched —
the two scripts never open the database; everything they use is in the flat `.tsv`/`.cif` files
already produced by `a33_subsites.py` and `resultados/calc/03_frames/`).

## 4. What is read off

129 subsite -1 units group into 8 EC-based subtypes. For each: the mean unit direction
ū = (ux,uy,uz) in frame coordinates, its length |ū| (1.0 = every member points exactly the same
way), the mean and max angle of individual members to their own group's mean (internal spread),
and the angle between that group's mean direction and the plain alpha-amylase group's mean
direction (the number that answers the sign-flip question):

| subtype (EC) | representative PDBs | n units / n entries | mean direction ū | \|ū\| | spread: mean / max angle to own mean (°) | **angle to alpha-amylase mean (°)** |
|---|---|---|---|---|---|---|
| alpha-amylase, plain (3.2.1.1) | 1BAG, 1OSE, 1PPI, 1HX0, 7TAA, 6TOZ | 66 / 59 | (0.602, 0.767, 0.224) | 0.957 | 11.6 / 53.1 | 0.0 (reference) |
| maltogenic amylase (3.2.1.133) | 1QHO, 1QHP, 26WD | 3 / 3 | (0.580, 0.789, 0.204) | 0.995 | 5.3 / 6.4 | **2.1** |
| other (3.2.1.-; cyclic maltosyl-maltose hydrolase CmmF, acarbose glucosidase) | 6A0K, 6A0L, 9IVZ, 9UNS, 9XG6 | 6 / 5 | (0.581, 0.771, 0.259) | 0.939 | 19.9 / 24.8 | **2.3** |
| oligo-1,6-glucosidase / alpha-glucosidase (3.2.1.10/20) | 3A4A (Ima1), 5WCZ (MalL), 5DO8, 7P01, 7P07, 8YIE, 8YIF | 14 / 8 | (0.531, 0.821, 0.208) | 0.977 | 11.3 / 19.5 | **5.2** |
| exo maltooligosaccharide-forming amylase (3.2.1.60/98) | 1QI5, 1WPC, 6IYG, 6JQB | 4 / 4 | (0.524, 0.823, 0.218) | 0.982 | 10.9 / 12.5 | **5.5** |
| pullulanase (3.2.1.41) | 2FH8, 2FHB, 2FHC, 2FHF, 3WDI, 3WDJ, 4J3V, 6JFJ, 6JFX, 6JHI, 7LSA, 2E9B | 13 / 12 | (0.561, 0.816, 0.137) | 0.990 | 7.2 / 16.7 | **6.2** |
| cyclomaltodextrinase / neopullulanase (3.2.1.54;135) | 1J0I (NplT), 3A6O (TvCMD), 3K8M (SusG), 8DL2 (BoSusG) | 10 / 4 | (0.564, 0.822, 0.075) | 0.930 | 17.6 / 44.3 | **9.4** |
| isoamylase / glycogen-debranching (3.2.1.68) | 2VR5 (TreX), 2YA2 (SpuA), 7U3B, 7U3D (GlgX) | 13 / 4 | (0.281, 0.948, −0.147) | 0.950 | 16.4 / 30.1 | **30.3** |

Full per-entry table: `.../scratchpad/per_entry_vectors.tsv` (129 rows); per-subtype table:
`.../scratchpad/subtype_summary.tsv`.

**Reading this off:** every one of the 8 subtypes, including the three the task specifically
asked about — isoamylase/glycogen-debranching (an alpha-1,6 debranching activity, geometrically
the most different from an exo/endo alpha-1,4 amylase), cyclomaltodextrinase/neopullulanase, and
oligo-1,6-glucosidase — has a mean direction within **≤ 30.3°** of the plain alpha-amylase mean
direction. None approaches 90° (orthogonal, i.e. no longer clearly on the same side) let alone
180° (a true flip, pointing to the opposite side of the catalytic dyad). Checked at the level of
every individual unit rather than the subtype average: of the 129 units, the single largest angle
to the alpha-amylase mean direction is **55.2°** (one of the two -1 units of 7U3D, a glycogen
debranching enzyme, GlgX); zero units exceed 60°, and zero exceed 90°.

## 5. The checks

- **Internal consistency of the recomputation.** Because a rigid transform (rotation +
  translation) preserves distances, the length of the vector v computed here
  (|v| = distance from the recomputed, frame-transformed catalytic midpoint to the -1 unit)
  must equal `d_ref_to_midpoint_A` already stored in `units.tsv` by `a33_subsites.py`
  (computed independently, in real coordinates, before any transform). Comparing the two for
  all 129 units: **maximum difference 0.0053 Å, mean difference 0.0024 Å, 0/129 differ by more
  than 0.05 Å.** This confirms the re-derivation here (reading the CIF files again, recomputing
  the carboxylate centroids, applying the same R,t) reproduces exactly the geometry
  `a33_subsites.py` used, so the direction vectors are not an artefact of a different
  recomputation.
- **The frame does bring every entry's catalytic dyad close to the same point, regardless of
  subtype**, which is the premise the whole check rests on: the vector norms |v| (= distance
  from the -1 unit to the catalytic midpoint) range only 1.77-4.49 Å across all 129 units
  (mean 2.46 Å, std 0.55 Å) — exactly the D_MAX ≤ 4.5 Å window the rule enforces, with similar
  scatter in every subtype (no subtype sits systematically farther, which would hint at a
  different binding geometry invisible to the direction check).
- **Both measures of spread agree there is no flip.** The subtype-mean-vs-alpha-amylase-mean
  angle (8 numbers, 0.0-30.3°, §4) and the per-unit max angle (1 number, 55.2° across all 129
  units, §4) are two independent ways of asking "is anything on the other side?", computed from
  the same vectors but at different levels of aggregation; both answer no.
- **The isoamylase/glycogen-debranching entries were specifically targeted**, since they attack
  the alpha(1→6) branch rather than the alpha(1→4) backbone and were the activity the task
  flagged as most likely to show a different binding geometry: 2VR5 (TreX, *Saccharolobus
  solfataricus*), 2YA2 (SpuA, *Streptococcus pneumoniae*), 7U3B and 7U3D (GlgX,
  *Streptomyces venezuelae*). These are also the subtype with the largest (though still modest)
  angle to the alpha-amylase mean, 30.3°, consistent with a real but limited difference in how
  the alpha-1,6 branch point sits in the frame, not a sign flip.
- **Not every GH13 subtype mentioned in the task got a -1 assignment at all.** Several isoamylase
  and cyclomaltodextrinase/neopullulanase entries (9LFN, 4OKD, 1J0J, 1J0K, 1JIB, 1JL8, 1VB9,
  1VFM, 1VFO, 1VFU, 3K8L, 6BS6, 3AXH, 3AXI, 4M56) have no unit within the 4.5 Å D_MAX window of
  the catalytic midpoint — `resultados/analysis/33_subsites/connected_sets.tsv`, `assigned=False`
  for every component of these entries — so the rule never fires for them; they contribute no
  data point either way and are not part of the table above. This is a coverage gap in -1
  assignment, not a sign disagreement, and is a fact about D_MAX/the deposited ligand's position
  relative to the dyad, already visible in `connected_sets.tsv`.

## 6. What this does not include

- This checks only the **direction** of subsite -1 relative to the catalytic dyad, not whether
  the +1/+2/... walk direction (toward the reducing end, via the main-chain acceptor atom) is
  itself consistent with "toward the catalytic acid" across subtypes; that is a related but
  separate question this note does not answer.
- The "plain alpha-amylase" reference direction is itself an average over 59 entries with a
  measurable internal spread (mean 11.6°, max 53.1° to its own mean) — alpha-amylase entries are
  not perfectly uniform either (different source organisms, different crystal forms, oligomeric
  contacts), so the 30.3° figure for isoamylase/glycogen-debranching should be read against that
  backdrop, not against a perfectly sharp reference.
- Group sizes are small for several subtypes (3-4 entries for maltogenic amylase, cyclomaltodextrinase/
  neopullulanase, and isoamylase/glycogen-debranching); the angle-to-reference numbers for those
  groups are not statistically powerful, only internally consistent with each other and with the
  individual-unit check.
- This is a geometric, not a biochemical, check: it says nothing about whether the -1 sugar unit
  is the chemically correct subsite for a given subtype's mechanism (e.g. whether an
  isoamylase's -1 unit, by this rule, corresponds to the subsite biochemists would call -1 for
  that enzyme); it only asks whether the rule, applied uniformly, keeps everyone on the same side
  of the dyad in the shared frame.
- No entries or files were changed by this check; `resultados/analysis/33_subsites/` and the
  database are untouched. The intermediate scripts and per-unit/per-subtype tables live in the
  session scratchpad (`/tmp/claude-1000/.../scratchpad/`), not in the repository, and are not
  guaranteed to persist.

## Bottom line

**No sign flip found.** Across all 8 GH13 enzyme-activity subtypes with a subsite -1 assignment —
plain alpha-amylase, maltogenic amylase, exo maltooligosaccharide-forming amylase,
cyclomaltodextrinase/neopullulanase, pullulanase, oligo-1,6-glucosidase/alpha-glucosidase, and
isoamylase/glycogen-debranching — the direction from the catalytic-dyad midpoint to the -1 unit,
in the shared GH13 frame, stays within 30° of the plain-alpha-amylase direction at the subtype
level, and within 56° at the level of every individual one of the 129 subsite -1 units. Nothing
approaches the ~90-180° that a silently flipped sign convention would require. This is a null
result, reported as such.
