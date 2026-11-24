# Figure 2 structure panels: rendering variants (2026-11-24)

All atoms come from the deposited mmCIF files (resultados/census/cif/1G94.cif, 3WEO.cif; model 1, first altloc,
no H). Residue and sugar selections come from scripts/figures/fig2_v2_struct.py. A hydrogen bond is drawn when a
side-chain N/O lies within 3.5 A of a sugar O/N; for each side-chain atom only the closest partner is drawn.

## Zoomed pair (current candidate, the "fewer residues, larger" direction)
- `zoom_pair.png` (7 in, 300 dpi; `_hires` 600 dpi), `zoom_GH13.png`, `zoom_GH31.png` (3.45 in, 300 dpi; `_hires` 900 dpi).
- Script: `scripts/figures/fig2_zoom.py` (needs `fig2_style_common.py` and `fig2_v2_struct.py`).
- Content: the -1 and +1 sugar units only. GH13 1G94: Asp179, Glu204, Asp289, His288 (H), Tyr51 (Ar),
  Ala at the Phe180 position (Al). GH31 3WEO: Asp327, Asp443, Asp542 (H), Phe575 (Ar), Ile328 (Al).
- Rendering: ChimeraX 1.11.1, offscreen. Lighting soft (ambient occlusion), no depth cue, dull material,
  silhouettes 2 px dark grey. Side chains are drawn from CB (Ala CB as a small ball), stick radius 0.12,
  carbons coloured by class (Okabe-Ito), N and O in standard colours. The sugar is thin ball-and-stick with grey
  carbons. Hydrogen-bond dashes are thin and drawn only to the -1/+1 units. No cartoon. Mild perspective
  (16 deg). Labels are a matplotlib overlay with leader lines, at 8.5 pt bold.
- View: a search over 20,000 directions. It penalises crossing bonds, atoms of different groups that meet on the
  page, sugar atoms drawn over each other, and residue atoms inside a ring outline. It requires the -1 ring to be
  tilted and the aromatic ring not to be edge-on. GH31 uses the top-scoring view; Asp327 and His600 gave the same
  view, and Asp327 was kept because it makes two H-bonds. The GH13 view was picked by eye from the 8 best distinct
  directions, because the top score hid the -1 ring behind Asp289.
- Known residual overlap in GH13: Glu204 crosses the edge of the +1 ring, and the His288 dash passes behind
  Asp289. Six residues surround the -1 ring on all sides, and none of the searched views was completely free.

## Rejected first variants (whole site, all 10 GH13 residues; the author judged them worse than v2)
- `style_1.png`: PyMOL ray trace with outlines, thin sticks, sugar as ball-and-stick, faint backbone trace,
  depth cue. Script `fig2_style_1_pymol.py`.
- `style_2.png`: ChimeraX lighting soft with silhouettes, thin sticks, sugar as ball-and-stick, no cartoon.
  Script `fig2_style_2_chimerax.py`.
- Styles 3 (pocket surface, `fig2_style_3_pocket.py`) and 4 (flat line art, `fig2_style_4_lineart.py`) were
  stopped before rendering when the direction changed. The 2D interaction diagram was not made: PLIP and
  LigPlot+ are not installed, and RDKit is available.
