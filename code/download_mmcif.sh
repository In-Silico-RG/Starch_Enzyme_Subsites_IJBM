#!/usr/bin/env bash
# Download the deposited mmCIF files (and, with --validation, the wwPDB validation XML reports) of every entry
# listed in data/census/cif_entries.txt from the RCSB PDB, into data/census/cif/ and data/census/validation/.
# These files are not distributed with the package (711 MB and 31 MB). Note: the PDB re-issues files after
# remediation; the analysis used the versions downloaded 2026-11-11 to 2026-11-23.
set -euo pipefail
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LIST="$PKG/data/census/cif_entries.txt"
CIF="$PKG/data/census/cif"; VAL="$PKG/data/census/validation"
mkdir -p "$CIF"
while read -r p; do
  [ -z "$p" ] && continue
  [ -s "$CIF/$p.cif" ] || curl -sfL --retry 3 -o "$CIF/$p.cif" "https://files.rcsb.org/download/$p.cif" || echo "failed: $p"
  if [ "${1:-}" = "--validation" ]; then
    mkdir -p "$VAL"; l=$(echo "$p" | tr 'A-Z' 'a-z'); m=${l:1:2}
    [ -s "$VAL/$p.xml.gz" ] || curl -sfL --retry 3 -o "$VAL/$p.xml.gz" \
      "https://files.rcsb.org/pub/pdb/validation_reports/$m/$l/${l}_validation.xml.gz" || echo "failed validation: $p"
  fi
done < "$LIST"
echo "$(ls "$CIF" | wc -l) mmCIF files in $CIF"
