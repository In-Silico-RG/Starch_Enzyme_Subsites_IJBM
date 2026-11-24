"""Ligands of the in-scope enzyme structures: every ligand, in full detail (AFC 2026-11-23).

Problem: the Part 2 table lists ligands only by name, one column per entry. To know what is bound in these
enzymes, and later which bound sugars have linkages that can be measured (including non-glucose units such
as those of acarbose, which the census of 2026-11-11 skipped), each ligand must be recorded with its
chemical components, how they are linked, and where it sits.
Data: RCSB Data API, GraphQL https://data.rcsb.org/graphql, batches of 50 entries. In-scope entries as in the
census (part2_structures_all.tsv in_scope == yes, glycogen phosphorylases of GT35 excluded, option C).
  - branched entities (oligosaccharides): description, residue list (set of chem comp ids), number of residues, BIRD reference
    molecule (prd), descriptors (WURCS; LINUCS, which carries the linkages), chains (asym ids), number of copies;
  - non-polymer entities: chem comp id, name, type, formula, SMILES, chains, number of copies.
Classification (stated rules, not tuned):
  glucan          branched entity whose residues are all GLC/BGC
  sugar_other     branched entity with any other residue (acarbose units, deoxy/amino/thio sugars, NAG...)
  n_glycan        branched entity made only of NAG/NDG/MAN/BMA/FUC/GAL/XYS... containing NAG/NDG
  mono_sugar      non-polymer whose chem comp type contains "saccharide"
  common          buffer, ion, cryoprotectant (list COMMON below)
  other           any other non-polymer (inhibitors, nucleotides, cofactors...)
Output: ../data/ligands/ligands_all.tsv (one row per ligand entity per entry),
        ../data/ligands/ligands_summary.json, raw API answers in ../data/ligands/raw/ (cached).
"""
import csv
import json
import time
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

R = Path(__file__).resolve().parent.parent / "resultados"
OUT = R / "ligands"
RAW = OUT / "raw"
Q = """{entries(entry_ids:%s){rcsb_id
 branched_entities{rcsb_id rcsb_branched_entity{pdbx_description} pdbx_entity_branch{type rcsb_branched_component_count}
  chem_comp_monomers{chem_comp{id}} prd{pdbx_reference_molecule{prd_id name class type}}
  pdbx_entity_branch_descriptor{descriptor type}
  rcsb_branched_entity_container_identifiers{asym_ids auth_asym_ids}}
 nonpolymer_entities{rcsb_id nonpolymer_comp{chem_comp{id name type formula}
  rcsb_chem_comp_descriptor{SMILES}}
  rcsb_nonpolymer_entity_container_identifiers{asym_ids auth_asym_ids}}}}"""
COMMON = {"HOH", "SO4", "PO4", "CL", "NA", "K", "CA", "MG", "ZN", "MN", "ACT", "GOL", "EDO", "PEG", "PG4",
          "TRS", "MES", "EPE", "IMD", "FMT", "NI", "CO", "CD", "BR", "IOD", "NO3", "DMS", "1PE", "P6G",
          "MPD", "PGE", "BTB", "ACY", "CIT", "TLA", "MLI", "SCN", "AZI", "NH4", "CU", "FE", "FE2", "PLP",
          "FLC", "B3P", "DTT", "BME", "PE4", "PE8", "2PE", "12P", "15P", "P33", "SIN", "MOH", "IPA", "EOH",
          "BU3", "PGO", "PGR", "TAR", "NHE", "CXS", "CAC", "HEZ", "OXY", "MRD", "UNX", "UNL", "SR", "CS", "RB", "LI"}
GLUC = {"GLC", "BGC"}
NGLY = {"NAG", "NDG", "MAN", "BMA", "FUC", "FUL", "GAL", "GLA", "XYS", "XYP", "SIA", "A2G", "NGA"}


def in_scope():
    rows = [r for r in csv.DictReader(open(R / "part2_structures_all.tsv"), delimiter="\t") if r["in_scope"] == "yes"]
    gp = {r["pdb"] for r in rows if r["family"] == "GT35" and "glycogen" in r["protein"].lower()}
    det = {r["pdb"] for r in csv.DictReader(open(R / "part2_pdb_details.tsv"), delimiter="\t")}
    info = {}
    for r in rows:
        if r["pdb"] in gp or r["pdb"] not in det:
            continue
        info.setdefault(r["pdb"], dict(family=r["family"], EC=r["EC"], protein=r["protein"], organism=r["organism"]))
    return info


def post(ids):
    f = RAW / f"{ids[0]}_{ids[-1]}.json"
    if f.exists():
        return json.loads(f.read_text())
    body = json.dumps({"query": Q % json.dumps(ids)}).encode()
    for a in range(5):
        try:
            req = urllib.request.Request("https://data.rcsb.org/graphql", data=body,
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=180) as r:
                d = json.load(r)
            if "errors" in d and not d.get("data"):
                raise RuntimeError(d["errors"][:1])
            ent = d["data"]["entries"]; f.write_text(json.dumps(ent)); return ent
        except Exception as ex:  # noqa: BLE001
            print("retry", a, ex, flush=True); time.sleep(5 * (a + 1))
    raise RuntimeError(f"RCSB batch failed: {ids[:3]}")


def classify_branched(comps):
    s = set(comps)
    if s <= GLUC:
        return "glucan"
    if s <= NGLY and s & {"NAG", "NDG"}:
        return "n_glycan"
    return "sugar_other"


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    info = in_scope(); ids = sorted(info)
    rows = []; returned = set()
    for i in range(0, len(ids), 50):
        for e in post(ids[i:i + 50]) or []:
            if e is None:
                continue
            p = e["rcsb_id"]; returned.add(p); m = info[p]
            for b in e.get("branched_entities") or []:
                comps = sorted({x["chem_comp"]["id"] for x in (b.get("chem_comp_monomers") or [])})
                nres = (b.get("pdbx_entity_branch") or {}).get("rcsb_branched_component_count") or 0
                desc = {d["type"]: d["descriptor"] for d in (b.get("pdbx_entity_branch_descriptor") or [])}
                prd = ((b.get("prd") or {}).get("pdbx_reference_molecule") or {})
                links = desc.get("LINUCS", "")
                ci = b.get("rcsb_branched_entity_container_identifiers") or {}
                rows.append(dict(pdb=p, family=m["family"], EC=m["EC"], protein=m["protein"], organism=m["organism"],
                                 kind="branched", ligand_class=classify_branched(comps), entity=b["rcsb_id"],
                                 comp_ids="-".join(comps), n_residues=nres,
                                 prd=f'{prd.get("prd_id", "")} {prd.get("name", "")}'.strip(),
                                 name=(b.get("rcsb_branched_entity") or {}).get("pdbx_description") or "",
                                 links=" ".join(links), formula="", smiles="",
                                 wurcs=desc.get("WURCS", ""), copies=len(ci.get("asym_ids") or []),
                                 chains=",".join(ci.get("auth_asym_ids") or [])))
            for n in e.get("nonpolymer_entities") or []:
                cc = n["nonpolymer_comp"]["chem_comp"]; cid = cc["id"]
                smi = (n["nonpolymer_comp"].get("rcsb_chem_comp_descriptor") or {}).get("SMILES") or ""
                ci = n.get("rcsb_nonpolymer_entity_container_identifiers") or {}
                cls = "common" if cid in COMMON else ("mono_sugar" if "saccharide" in (cc["type"] or "").lower() else "other")
                rows.append(dict(pdb=p, family=m["family"], EC=m["EC"], protein=m["protein"], organism=m["organism"],
                                 kind="nonpolymer", ligand_class=cls, entity=n["rcsb_id"], comp_ids=cid, n_residues=1, prd="",
                                 name=cc["name"] or "", links="", formula=cc["formula"] or "", smiles=smi, wurcs="",
                                 copies=len(ci.get("asym_ids") or []), chains=",".join(ci.get("auth_asym_ids") or [])))
        print(min(i + 50, len(ids)), "/", len(ids), flush=True)
        time.sleep(0.3)
    cols = ["pdb", "family", "EC", "protein", "organism", "kind", "ligand_class", "entity", "comp_ids", "n_residues",
            "name", "prd", "links", "formula", "smiles", "wurcs", "copies", "chains"]
    with open(OUT / "ligands_all.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, cols, delimiter="\t"); w.writeheader()
        for r in rows:
            w.writerow({k: " ".join(str(v).split()) for k, v in r.items()})
    # summary
    by_class = Counter(r["ligand_class"] for r in rows)
    entries_with = {c: len({r["pdb"] for r in rows if r["ligand_class"] == c}) for c in by_class}
    comp_count = Counter(); comp_entries = defaultdict(set); comp_name = {}
    for r in rows:
        if r["ligand_class"] == "common":
            continue
        for c in r["comp_ids"].split("-"):
            comp_count[c] += 1; comp_entries[c].add(r["pdb"])
        if r["kind"] == "nonpolymer":
            comp_name[r["comp_ids"]] = r["name"]
    non_glc_units = Counter()
    for r in rows:
        if r["ligand_class"] == "sugar_other":
            for c in set(r["comp_ids"].split("-")) - GLUC:
                non_glc_units[c] += 1
    fam = defaultdict(Counter)
    for r in rows:
        fam[r["family"]][r["ligand_class"]] += 1
    summ = dict(source="RCSB Data API GraphQL, queried 2026-11-23", entries_in_scope=len(ids), entries_returned=len(returned),
                not_returned=sorted(set(ids) - returned), ligand_entities=len(rows), by_class=dict(by_class),
                entries_with_class=entries_with,
                components_by_entries=[(c, len(comp_entries[c]), comp_name.get(c, "")) for c, _ in
                                       sorted(comp_entries.items(), key=lambda x: -len(x[1]))[:60]],
                non_glucose_units_in_sugar_other=non_glc_units.most_common(),
                by_family={k: dict(v) for k, v in sorted(fam.items())})
    (OUT / "ligands_summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps({k: summ[k] for k in ("entries_in_scope", "entries_returned", "ligand_entities", "by_class",
                                           "entries_with_class")}, indent=1))


if __name__ == "__main__":
    main()
