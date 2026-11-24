"""Part 2, step 2: details of every PDB entry of the Part 2 inventory, from the RCSB PDB (AFC 2026-11-11).

Problem: CAZy lists only carbohydrate ligands and resolution. We also need the method, the release date,
every ligand (including non-carbohydrate inhibitors) and the length of each bound oligosaccharide.
Data: RCSB Data API, GraphQL endpoint https://data.rcsb.org/graphql, queried in batches of 100 entries.
  - branched_entities: oligosaccharides deposited as branched entities, with their description and
    number of sugar residues (rcsb_branched_component_count)
  - nonpolymer_entities: single-residue ligands (ions, buffers, monosaccharides, inhibitors), chem_comp id/name/type
Output: ../data/part2_pdb_details.tsv, one row per PDB entry (all entries of part2_structures_all.tsv).
  longest_glucan = largest residue count among branched entities whose description names only
  alpha/beta-D-glucopyranose units (rule corrected 2026-11-23, see glucan_only); acarbose-type units are
  reported in the ligand list. `--rederive` recomputes this column from the stored descriptions.
"""
import csv
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

R = Path(__file__).resolve().parent.parent / "resultados"
Q = """{entries(entry_ids:%s){rcsb_id exptl{method} rcsb_accession_info{initial_release_date}
 rcsb_entry_info{resolution_combined}
 branched_entities{rcsb_branched_entity{pdbx_description} pdbx_entity_branch{rcsb_branched_component_count}}
 nonpolymer_entities{nonpolymer_comp{chem_comp{id name type}}}}}"""
COMMON = {"HOH", "SO4", "PO4", "CL", "NA", "K", "CA", "MG", "ZN", "MN", "ACT", "GOL", "EDO", "PEG", "PG4",
          "TRS", "MES", "EPE", "IMD", "FMT", "NI", "CO", "CD", "BR", "IOD", "NO3", "DMS", "1PE", "P6G",
          "MPD", "PGE", "BTB", "ACY", "CIT", "TLA", "MLI", "SCN", "AZI", "NH4", "CU", "FE", "FE2", "PLP"}


def post(ids):
    body = json.dumps({"query": Q % json.dumps(ids)}).encode()
    for a in range(5):
        try:
            req = urllib.request.Request("https://data.rcsb.org/graphql", data=body,
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)["data"]["entries"]
        except Exception:
            time.sleep(5 * (a + 1))
    raise RuntimeError(f"RCSB batch failed: {ids[:3]}...")


def glucan_only(desc):
    """True when every residue named in the description is alpha- or beta-D-glucopyranose.
    Corrected 2026-11-23: the first version split on every hyphen, so names such as
    '2-acetamido-2-deoxy-beta-D-glucopyranose' (NAG) or '4,6-dideoxy-4-...-alpha-D-glucopyranose' (acarbose unit)
    passed as glucose. Residues are now split at the linkage markers '(1-4)', '(1-6)' and at branch brackets;
    agrees with the chemical-component composition (data/ligands/ligands_all.tsv) in all 645 in-scope entries."""
    d = re.sub(r"(?i)^cyclic\s+", "", desc)
    toks = [t.strip(" ()").lower() for t in re.split(r"-?\(\d+-\d+\)-?|\[|\]", d)]
    toks = [t for t in toks if t]
    seen = False
    for t in toks:
        if t in ("alpha-d-glucopyranose", "beta-d-glucopyranose"):
            seen = True
        elif not re.fullmatch(r"cyclo\w*kis", t):
            return False
    return seen


def rederive():
    """Recompute longest_glucan_residues from the descriptions stored on 2026-11-11 (no new query)."""
    f = R / "part2_pdb_details.tsv"
    rows = list(csv.reader(open(f), delimiter="\t"))
    for r in rows[1:]:
        n = 0
        for part in [x for x in r[5].split(" | ") if x]:
            m = re.match(r"(.*) \[(\d+)\]$", part)
            if m and glucan_only(m.group(1)):
                n = max(n, int(m.group(2)))
        r[4] = str(n)
    with open(f, "w") as g:
        for r in rows:
            g.write("\t".join(r) + "\n")


def main():
    ids = sorted({r["pdb"] for r in csv.DictReader(open(R / "part2_structures_all.tsv"), delimiter="\t")})
    rows = []
    for i in range(0, len(ids), 100):
        for e in post(ids[i:i + 100]) or []:
            if e is None:
                continue
            br = e.get("branched_entities") or []
            sugars = [(b["rcsb_branched_entity"]["pdbx_description"] or "",
                       (b.get("pdbx_entity_branch") or {}).get("rcsb_branched_component_count") or 0) for b in br]
            gl = [n for d, n in sugars if glucan_only(d)]
            nps = [(n["nonpolymer_comp"]["chem_comp"]["id"], n["nonpolymer_comp"]["chem_comp"]["name"])
                   for n in (e.get("nonpolymer_entities") or [])]
            ligs = [f"{i_}:{nm}" for i_, nm in nps if i_ not in COMMON]
            res = (e.get("rcsb_entry_info") or {}).get("resolution_combined") or [""]
            rows.append([e["rcsb_id"], ";".join(x["method"] for x in e.get("exptl") or []),
                         (e.get("rcsb_accession_info") or {}).get("initial_release_date", "")[:10],
                         str(res[0]), str(max(gl) if gl else 0),
                         " | ".join(f"{d} [{n}]" for d, n in sugars), " | ".join(ligs)])
        print(min(i + 100, len(ids)), "/", len(ids), flush=True)
        time.sleep(0.5)
    with open(R / "part2_pdb_details.tsv", "w") as f:
        f.write("pdb\tmethod\trelease_date\tresolution_A\tlongest_glucan_residues\tbranched_oligosaccharides\t"
                "other_ligands\n")
        for r in rows:
            f.write("\t".join(" ".join(x.split()) for x in r) + "\n")
    print(len(rows), "entries written;", len(ids) - len(rows), "not returned")


if __name__ == "__main__":
    rederive() if "--rederive" in sys.argv else main()
