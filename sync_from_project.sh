#!/usr/bin/env bash
# Refresh code/ and data/ of this package from the project folder, so the copies never drift silently.
#
# Usage:  PROJECT=<path to Proyecto_F_Enzimas_Amiloliticas> ./sync_from_project.sh
#         (default PROJECT: the parent folder of this package)
#
# Every source -> destination pair is listed below. After copying, two transformations are applied to the copies
# only (never to the project files):
#   1. path shim (code/ only): the project layout scripts/ + resultados/ becomes code/ + data/ here, so
#      'resultados/' -> 'data/', 'scripts/{analysis,calc,figures}' -> 'code/{...}', figure outputs
#      'manuscript_ijbm/figures' and 'manuscript_short/figures' -> 'figures'; absolute project paths in the
#      grip_independence_check steps -> a path computed from the script location.
#   2. scrub (code/ and data/ text files): local and cluster home paths and cluster account names are replaced by
#      placeholders (<PROJECT_ROOT>, <LOCAL_DISK>, ~/, <cluster-account>).
# At the end the script lists any file in code/ or data/ that still matches a private pattern; it must print nothing.
set -euo pipefail
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT="${PROJECT:-$(cd "$PKG/.." && pwd)}"
S="$PROJECT/scripts"
R="$PROJECT/resultados"
C="$PKG/code"
D="$PKG/data"
[ -d "$S" ] && [ -d "$R" ] || { echo "PROJECT=$PROJECT has no scripts/ and resultados/"; exit 1; }

cp1() {  # cp1 SOURCE_FILE DEST_FILE
  mkdir -p "$(dirname "$2")"; cp -p "$1" "$2"; echo "  ${1#$PROJECT/} -> ${2#$PKG/}" >> "$PKG/.sync_manifest"
}
cpdir() {  # cpdir SOURCE_DIR DEST_DIR [rsync excludes...]
  local src="$1" dst="$2"; shift 2
  mkdir -p "$dst"; rsync -a --delete "$@" "$src/" "$dst/"; echo "  ${src#$PROJECT/}/ -> ${dst#$PKG/}/  ($*)" >> "$PKG/.sync_manifest"
}
: > "$PKG/.sync_manifest"

echo "== code"
# --- database and inventory builders (provenance of data/db/enzimas_experimental.sqlite)
for f in part1_table.py part2_structures.py part2_pdb_details.py ligands_inventory.py part3_measure.py \
         census_bound_linkages.py build_database.py add_sifts.py add_chembl.py; do
  cp1 "$S/$f" "$C/$f"
done
# --- supporting calculations: common helpers, frames (calc. 3), PROPKA (calc. 1), MM/GBSA series (calc. 9)
for f in common.py c03_frames.py c01_propka.py c01_collect.py c09_series.py c09_series_summary.py; do
  cp1 "$S/calc/$f" "$C/calc/$f"
done
# --- analysis: subsites, contact kinds, control, ring distortion, MM/GBSA per subsite, rule test, hydroxyl map, RCSB queries
for f in a33_subsites.py a40_contact_types.py a41_control_noncatalytic.py a42_ring_distortion_subsite.py \
         a43_mmgbsa_subsite.py a44a_build.py a44b_tests.py a44c_tables.py a45_hydroxyl_map.py rcsb_q.py rcsb_q2.py; do
  cp1 "$S/analysis/$f" "$C/analysis/$f"
done
# --- independence check (protein identity clusters, ident.json used by a44a_build.py); scripts live in resultados/ in the project
for f in step1.py step2.py step3.py step4.py step5.py; do
  cp1 "$R/analysis/grip_independence_check/$f" "$C/analysis/grip_independence_check/$f"
done
# --- figures
for f in fig1_protein_common.py fig1_protein_chem.py fig1_protein_frame.py fig1_protein_struct.py fig1_protein_layouts.py \
         fig2_v2_data.py fig2_v2_struct.py fig2_v2_layouts.py fig2_style_common.py fig2_style_1_pymol.py \
         fig2_style_2_chimerax.py fig2_style_3_pocket.py fig2_style_4_lineart.py fig2_zoom.py fig_contact_classes.py; do
  cp1 "$S/figures/$f" "$C/figures/$f"
done

echo "== data"
# --- database (41 MB, included whole) and its summaries
cp1 "$R/db/enzimas_experimental.sqlite" "$D/db/enzimas_experimental.sqlite"
for f in db_summary.json sifts_summary.json chembl_summary.json; do cp1 "$R/db/$f" "$D/db/$f"; done
# --- inventory inputs of the builders (small tables)
for f in cazy_20261111_activities.tsv cazy_20261111_families.tsv cazy_acarbose_families.tsv part1_starch_breaking_activities.tsv \
         part1_excluded_activities.tsv part2_structures_all.tsv part2_summary_family_ec.tsv part2_pdb_details.tsv; do
  cp1 "$R/$f" "$D/$f"
done
cp1 "$R/ligands/ligands_all.tsv" "$D/ligands/ligands_all.tsv"
cp1 "$R/ligands/ligands_summary.json" "$D/ligands/ligands_summary.json"
for f in mcsa_starch_entries.json part3_site_geometry.tsv part3_step1_catalytic_residues.tsv; do cp1 "$R/part3/$f" "$D/part3/$f"; done
# --- census of bound linkages (ring puckering input of a42)
cp1 "$R/census/bound_linkages.tsv" "$D/census/bound_linkages.tsv"
cp1 "$R/census/bound_linkages_checked.tsv" "$D/census/bound_linkages_checked.tsv"
# --- calculation 3: frames, catalytic residues, conserved water/ion sites, reference mmCIF of each frame
cpdir "$R/calc/03_frames" "$D/calc/03_frames"
# --- calculation 1: PROPKA tables (per-entry runs/ left out: 601 MB, regenerable with code/calc/c01_propka.py)
for f in README.md propka_all.tsv propka_near_ligand.tsv entries.txt guane_entries.txt; do
  cp1 "$R/calc/01_propka/$f" "$D/calc/01_propka/$f"
done
# --- calculation 9: MM/GBSA series, inputs and results only (topologies, coordinates, restart files left out)
cp1 "$R/calc/09_mmgbsa/README.md" "$D/calc/09_mmgbsa/README.md"
cp1 "$R/calc/09_mmgbsa/series/series_summary.tsv" "$D/calc/09_mmgbsa/series/series_summary.tsv"
for d in 5NN5_NOJ 5NN6_MIG 9GTC_A1IO7 9GTW_A1IPC; do
  for f in FINAL_RESULTS_MMPBSA.dat FINAL_DECOMP_MMPBSA.dat complex_in.pdb reduce_out.pdb build.log tleap.in min.in mmgbsa.in \
           lig_in.mol2 lig.mol2 lig.frcmod propka_in.pka; do
    [ -f "$R/calc/09_mmgbsa/series/$d/$f" ] && cp1 "$R/calc/09_mmgbsa/series/$d/$f" "$D/calc/09_mmgbsa/series/$d/$f"
  done
done
# --- analysis result tables
cpdir "$R/analysis/33_subsites" "$D/analysis/33_subsites"
cpdir "$R/analysis/40_contact_types" "$D/analysis/40_contact_types"
cpdir "$R/analysis/41_control_noncatalytic" "$D/analysis/41_control_noncatalytic"
cpdir "$R/analysis/42_ring_distortion" "$D/analysis/42_ring_distortion"
cpdir "$R/analysis/43_mmgbsa_subsite" "$D/analysis/43_mmgbsa_subsite"
cpdir "$R/analysis/44_rule_test" "$D/analysis/44_rule_test"
cpdir "$R/analysis/45_hydroxyl_map" "$D/analysis/45_hydroxyl_map"
cpdir "$R/analysis/grip_independence_check" "$D/analysis/grip_independence_check" --exclude='*.py'
# --- checks and literature reports cited for statements of the manuscript and SI
for f in carbpol_computational_survey.md substituted_starch_hydrolysis_literature.md conserved_sites_resolution_confound.md \
         ring_distortion_by_subsite.md propka_nucleophile_acid_check.md grip_independence_check.md grip_checks_2026-11-23.md \
         hydroxyl_ts_literature_check.md mutagenesis_literature_check.md gh13_gh31_clan_check.md family_identity_audit.md \
         subsite_sign_convention_check.md structure_activity_gap.md; do
  cp1 "$R/analysis/$f" "$D/analysis/$f"
done
# --- figure intermediates (renders *.png left out: regenerated by the figure scripts)
cpdir "$R/figures/fig2_v2_work" "$D/figures/fig2_v2_work" --exclude='*.png'
cpdir "$R/figures/fig2_style_work" "$D/figures/fig2_style_work" --exclude='*.png'
cpdir "$R/figures/fig2_work" "$D/figures/fig2_work" --exclude='*.png'

# --- derived listings: entries, mmCIF files to download, database schema
python3 - "$D" "$R" <<'PY'
import sqlite3, sys, os
D, R = sys.argv[1], sys.argv[2]
con = sqlite3.connect(f"file:{D}/db/enzimas_experimental.sqlite?mode=ro", uri=True)
with open(f"{D}/db/schema.sql", "w") as f:
    for (sql,) in con.execute("select sql from sqlite_master where sql is not null order by type desc, name"):
        f.write(sql + ";\n\n")
with open(f"{D}/pdb_entries.tsv", "w") as f:
    f.write("pdb\tfamily\tmethod\tresolution_A\trelease_date\n")
    for r in con.execute("select pdb, family, method, resolution_A, release_date from entry order by family, pdb"):
        f.write("\t".join("" if v is None else str(v) for v in r) + "\n")
cifs = sorted(x[:-4] for x in os.listdir(f"{R}/census/cif") if x.endswith(".cif"))
os.makedirs(f"{D}/census", exist_ok=True)
open(f"{D}/census/cif_entries.txt", "w").write("\n".join(cifs) + "\n")
print(f"  schema.sql, pdb_entries.tsv ({con.execute('select count(*) from entry').fetchone()[0]} entries), "
      f"census/cif_entries.txt ({len(cifs)} mmCIF files)")
PY

echo "== path shim (code/)"
find "$C" -name '*.py' -not -path '*/obsolete_build/*' -print0 | xargs -0 sed -i \
  -E -e "s#'/media/[^']*/Proyecto_F_Enzimas_Amiloliticas/#__import__('os').path.abspath(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', '..')) + '/#g" \
  -e "s#resultados/#data/#g" \
  -e "s#'resultados'#'data'#g" \
  -e "s#scripts/analysis#code/analysis#g; s#scripts/calc#code/calc#g; s#scripts/figures#code/figures#g" \
  -e "s#manuscript_ijbm/figures#figures#g; s#manuscript_short/figures#figures#g"

echo "== scrub (code/, data/)"
# Cluster account names are read from .scrub_accounts (one per line; kept out of the repository by .gitignore).
ACC=""; [ -f "$PKG/.scrub_accounts" ] && ACC=$(grep -v '^\s*$' "$PKG/.scrub_accounts" | paste -sd'|')
[ -n "$ACC" ] || echo "WARNING: no .scrub_accounts file; cluster account names are not scrubbed"
TEXT=( \( -name '*.py' -o -name '*.sh' -o -name '*.md' -o -name '*.tsv' -o -name '*.txt' -o -name '*.json' \
       -o -name '*.log' -o -name '*.out' -o -name '*.dat' -o -name '*.in' -o -name '*.pml' -o -name '*.cxc' -o -name '*.tex' \) )
find "$C" "$D" -type f "${TEXT[@]}" -print0 | xargs -0 sed -i -E \
  -e 's#/media/[^ ]*/Proyecto_F_Enzimas_Amiloliticas#<PROJECT_ROOT>#g' \
  -e 's#/media/[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+#<LOCAL_DISK>#g' \
  -e 's#/home/[A-Za-z0-9_.-]+/#~/#g'
[ -n "$ACC" ] && find "$C" "$D" -type f "${TEXT[@]}" -print0 | xargs -0 sed -i -E "s#\b($ACC)\b#<cluster-account>#g"
# mentions of the earlier submission and its reviews (not part of this record)
sed -i 's#^Both JACS Au referees said the draft.s claim is circular (2026-11-24 reviews in `docs/reviews/`):#The claim of the earlier draft was open to the objection that it is circular:#' \
  "$D/analysis/44_rule_test/REPORT.md"
sed -i 's#for the short draft (JACS Au, two-column width 7.0 in)#for the short draft (two-column width 7.0 in)#' "$C/figures/fig2_v2_layouts.py"

echo "== private-pattern check (must list nothing)"
PAT='/home/|/media/|referee|JACS Au|docs/reviews|[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[a-z]{2,}'
[ -n "$ACC" ] && PAT="$PAT|$ACC"
grep -rIE "$PAT" "$C" "$D" | grep -v 'aldo.combariza@unisucre.edu.co' | cut -c1-200 || true
echo "done; pairs in $PKG/.sync_manifest; size: $(du -sh "$PKG" --exclude=.git | cut -f1)"
