"""Experimental database of Project F: every molecule bound or attached in the in-scope structures (AFC 2026-11-23).

Problem: the information gathered so far is spread over many tables, and the ligands were known only by name.
To work with what the crystal structures contain, every molecule other than the protein chain and water must be
registered where it sits, how well it is determined, what it is bonded to and what surrounds it, in one place
together with the tables already built (Part 1 activities, Part 2 entries, Part 3 catalytic sites, census).
Data:
  - in-scope entries (645): part2_structures_all.tsv in_scope == yes, glycogen phosphorylases excluded (option C),
    returned by the RCSB in Part 2;
  - coordinates: mmCIF from files.rcsb.org (cache data/census/cif/, not in git);
  - fit to the density: wwPDB validation reports (cache data/census/validation/, not in git), per residue
    RSCC (real-space correlation coefficient) and RSR (real-space R);
  - chemical components: RCSB Data API (name, type, formula, weight, SMILES, InChIKey).
Method (model 1 of each entry, all alternate locations kept for atom counts, first one used for geometry):
  1. a molecule instance is every residue that is not a standard amino acid of a protein chain and not water:
     non-polymer molecules (ions, buffers, inhibitors, cofactors, single sugars), every residue of a branched
     entity (oligosaccharide), and modified residues inside a protein chain;
  2. bonds: the deposited _struct_conn records (covalent 'covale', metal 'metalc', disulfide 'disulf' when a
     partner is an instance), with the deposited distance;
  3. contacts: every protein residue, other instance and water with a heavy atom within 4.0 Å of a heavy atom of
     the instance (minimum distance kept), symmetry copies included; for entries whose coordinates are not in the
     standard crystal frame, contacts within the deposited coordinates only (listed in db_summary.json);
  4. classes (stated rules): water excluded; 'ion' = component whose formula is one atom; 'common' = buffer/cryoprotectant
     list of ligands_inventory.py; 'sugar' = component type contains 'saccharide' or residue of a branched entity;
     'n_glycan' = sugar covalently bonded to Asn ND2 or in a branched entity containing NAG/NDG bonded to protein;
     'modified_residue' = non-standard residue inside a protein chain; 'other' = the rest (inhibitors, cofactors).
Output: ../data/db/enzimas_experimental.sqlite and ../data/db/db_summary.json
"""
import csv
import gzip
import json
import sqlite3
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligands_inventory import COMMON  # noqa: E402

R = Path(__file__).resolve().parent.parent / "resultados"
CIF = R / "census" / "cif"
VAL = R / "census" / "validation"
OUT = R / "db"
DB = OUT / "enzimas_experimental.sqlite"
CUT = 4.0
AA = set(gemmi.find_tabulated_residue(x).name for x in
         "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())


def get(url, f, tries=4):
    if f.exists() and f.stat().st_size > 0:
        return True
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                f.write_bytes(r.read())
            return True
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return False
            time.sleep(5 * (a + 1))
        except Exception:  # noqa: BLE001
            time.sleep(5 * (a + 1))
    return False


def entries():
    det = {r["pdb"]: r for r in csv.DictReader(open(R / "part2_pdb_details.tsv"), delimiter="\t")}
    rows = [r for r in csv.DictReader(open(R / "part2_structures_all.tsv"), delimiter="\t") if r["in_scope"] == "yes"]
    gp = {r["pdb"] for r in rows if r["family"] == "GT35" and "glycogen" in r["protein"].lower()}
    info = {}
    for r in rows:
        p = r["pdb"]
        if p in gp or p not in det:
            continue
        e = info.setdefault(p, dict(pdb=p, method=det[p]["method"], resolution_A=det[p]["resolution_A"],
                                    release_date=det[p]["release_date"], family=set(), EC=set(), protein=set(),
                                    organism=set(), uniprot=set()))
        for k in ("family", "EC", "protein", "organism", "uniprot"):
            if r.get(k):
                e[k].add(r[k])
    for e in info.values():
        for k in ("family", "EC", "protein", "organism", "uniprot"):
            e[k] = " | ".join(sorted(e[k]))
    return info


def validation(p):
    f = VAL / f"{p}.xml.gz"
    if not get(f"https://files.rcsb.org/pub/pdb/validation_reports/{p[1:3].lower()}/{p.lower()}/"
               f"{p.lower()}_validation.xml.gz", f):
        return {}
    out = {}
    try:
        root = ET.fromstring(gzip.open(f).read())
    except Exception:  # noqa: BLE001
        return {}
    for s in root.iter("ModelledSubgroup"):
        a = s.attrib
        if a.get("model", "1") != "1":
            continue
        k = (a.get("chain"), a.get("resnum"), (a.get("icode") or "").strip(), a.get("resname"))
        if k in out and out[k].get("rscc") is not None:
            continue
        def fl(x):
            try:
                return float(a[x])
            except (KeyError, ValueError):
                return None
        out[k] = dict(rscc=fl("rscc"), rsr=fl("rsr"), avgoccu=fl("avgoccu"))
    return out


def chem_comps(ids):
    q = """{chem_comps(comp_ids:%s){chem_comp{id name type formula formula_weight}
            rcsb_chem_comp_descriptor{SMILES InChIKey}}}"""
    out = {}
    ids = sorted(ids)
    for i in range(0, len(ids), 200):
        body = json.dumps({"query": q % json.dumps(ids[i:i + 200])}).encode()
        for a in range(5):
            try:
                req = urllib.request.Request("https://data.rcsb.org/graphql", data=body,
                                             headers={"Content-Type": "application/json"})
                with urllib.request.urlopen(req, timeout=180) as r:
                    d = json.load(r)["data"]["chem_comps"]
                break
            except Exception as ex:  # noqa: BLE001
                print("retry chem_comps", ex, flush=True); time.sleep(5 * (a + 1))
        for c in d or []:
            if not c:
                continue
            cc = c["chem_comp"]; ds = c.get("rcsb_chem_comp_descriptor") or {}
            out[cc["id"]] = dict(comp_id=cc["id"], name=cc["name"], type=cc["type"], formula=cc["formula"],
                                 formula_weight=cc["formula_weight"], smiles=ds.get("SMILES"), inchikey=ds.get("InChIKey"))
    return out


def single_atom(formula):
    """True when the chemical formula of the component is one atom, e.g. 'Ca 2' or 'Cl -1' (ions);
    decided on the component, not on the number of modelled atoms (a partly modelled sugar can have one atom)."""
    if not formula:
        return False
    parts = [t for t in formula.split() if t[0].isalpha()]
    if len(parts) != 1:
        return False
    num = "".join(ch for ch in parts[0] if ch.isdigit())
    return num in ("", "1")


def import_tsv(con, table, path, rename=None):
    rows = list(csv.reader(open(path), delimiter="\t"))
    head = [(rename or {}).get(h, h) for h in rows[0]]
    head = [h.replace("(", "_").replace(")", "").replace(" ", "_").replace("-", "_").replace("/", "_") for h in head]
    con.execute(f"DROP TABLE IF EXISTS {table}")
    con.execute(f"CREATE TABLE {table} ({', '.join(f'{h!r}' for h in head)})")
    con.executemany(f"INSERT INTO {table} VALUES ({', '.join('?' * len(head))})",
                    [r + [None] * (len(head) - len(r)) for r in rows[1:]])
    return len(rows) - 1


def main():
    OUT.mkdir(parents=True, exist_ok=True); CIF.mkdir(parents=True, exist_ok=True); VAL.mkdir(parents=True, exist_ok=True)
    info = entries(); ids = sorted(info)
    if len(sys.argv) > 1:
        ids = [x for x in ids if x in sys.argv[1:]]
    inst, links, cont, comps, miss = [], [], [], set(), []
    water = {}
    nosym = []
    for n, p in enumerate(ids, 1):
        f = CIF / f"{p}.cif"
        if not get(f"https://files.rcsb.org/download/{p}.cif", f):
            miss.append(p); continue
        st = gemmi.read_structure(str(f)); st.setup_entities(); st.remove_hydrogens()
        model = st[0]
        val = validation(p)
        ent_type = {}
        for ent in st.entities:
            for sc in ent.subchains:
                ent_type[sc] = (ent.entity_type.name, ent.name)
        # instances
        key2id = {}
        nwat = 0
        for ch in model:
            for res in ch:
                et, eid = ent_type.get(res.subchain, ("Unknown", ""))
                if res.is_water():
                    nwat += 1; continue
                if et == "Polymer":
                    if res.name in AA or gemmi.find_tabulated_residue(res.name) and \
                            gemmi.find_tabulated_residue(res.name).is_nucleic_acid():
                        continue
                    kind = "modified_residue"
                elif et == "Branched":
                    kind = "branched"
                else:
                    kind = "nonpolymer"
                atoms = list(res)
                alts = sorted({a.altloc for a in atoms if a.altloc != "\0"})
                first = [a for a in atoms if a.altloc in ("\0", alts[0] if alts else "\0")]
                c = gemmi.Position(0, 0, 0)
                for a in first:
                    c += a.pos
                c /= max(len(first), 1)
                ic = res.seqid.icode.strip()
                v = val.get((ch.name, str(res.seqid.num), ic, res.name), {})
                iid = len(inst) + 1
                key2id[(ch.name, res.seqid.num, ic, res.name)] = iid
                inst.append(dict(instance_id=iid, pdb=p, auth_asym=ch.name, label_asym=res.subchain, auth_seq=res.seqid.num,
                                 ins_code=ic, comp_id=res.name, entity_id=eid, kind=kind, n_atoms=len(first),
                                 altlocs="".join(alts), mean_occupancy=round(sum(a.occ for a in first) / max(len(first), 1), 3),
                                 mean_b=round(sum(a.b_iso for a in first) / max(len(first), 1), 2),
                                 rscc=v.get("rscc"), rsr=v.get("rsr"), x=round(c.x, 3), y=round(c.y, 3), z=round(c.z, 3)))
                comps.add(res.name)
        water[p] = nwat
        # bonds as deposited
        for con_ in st.connections:
            a1, a2 = con_.partner1, con_.partner2
            k1 = (a1.chain_name, a1.res_id.seqid.num, a1.res_id.seqid.icode.strip(), a1.res_id.name)
            k2 = (a2.chain_name, a2.res_id.seqid.num, a2.res_id.seqid.icode.strip(), a2.res_id.name)
            if k1 not in key2id and k2 not in key2id:
                continue
            links.append(dict(pdb=p, conn_type=con_.type.name, instance_id_1=key2id.get(k1), chain_1=k1[0], seq_1=k1[1],
                              comp_1=k1[3], atom_1=a1.atom_name, instance_id_2=key2id.get(k2), chain_2=k2[0], seq_2=k2[1],
                              comp_2=k2[3], atom_2=a2.atom_name,
                              distance_A=round(con_.reported_distance, 3) if con_.reported_distance > 0 else None))
        # contacts within 4.0 A
        try:
            ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
        except RuntimeError:  # non-standard crystal frame: contacts within the deposited coordinates only
            ns = gemmi.NeighborSearch(model, gemmi.UnitCell(), 5).populate(); nosym.append(p)
        best = {}
        for ch in model:
            for res in ch:
                k = (ch.name, res.seqid.num, res.seqid.icode.strip(), res.name)
                if k not in key2id:
                    continue
                iid = key2id[k]
                for a in res:
                    for m in ns.find_atoms(a.pos, "\0", radius=CUT):
                        cra = m.to_cra(model)
                        r2 = cra.residue
                        k2 = (cra.chain.name, r2.seqid.num, r2.seqid.icode.strip(), r2.name)
                        if k2 == k:
                            continue
                        d = cra.atom.pos.dist(a.pos)
                        if r2.is_water():
                            pk = "water"
                        elif k2 in key2id:
                            pk = "instance"
                        else:
                            pk = "protein"
                        kk = (iid, k2)
                        if kk not in best or d < best[kk][0]:
                            best[kk] = (d, pk, a.name, cra.atom.name)
        for (iid, k2), (d, pk, an, bn) in best.items():
            cont.append(dict(instance_id=iid, partner_kind=pk, partner_instance_id=key2id.get(k2), chain=k2[0], seq=k2[1],
                             ins_code=k2[2], comp=k2[3], min_distance_A=round(d, 3), atom=an, partner_atom=bn))
        if n % 25 == 0:
            print(n, "/", len(ids), len(inst), "instances", flush=True)
    # chemical components
    cc = chem_comps(comps)
    # classes
    bonded_asn = {(l["pdb"], l["instance_id_1"] or l["instance_id_2"]) for l in links
                  if l["conn_type"] == "Covale" and ((l["comp_1"] == "ASN" and l["atom_1"] == "ND2") or
                                                     (l["comp_2"] == "ASN" and l["atom_2"] == "ND2"))}
    ent_has_nag = defaultdict(set)
    for i in inst:
        if i["kind"] == "branched":
            ent_has_nag[(i["pdb"], i["label_asym"])].add(i["comp_id"])
    by_id = {i["instance_id"]: i for i in inst}
    asn_groups = {(p_, by_id[j]["label_asym"]) for p_, j in bonded_asn if j in by_id}
    for i in inst:
        c = cc.get(i["comp_id"], {})
        typ = (c.get("type") or "").lower()
        if i["kind"] == "modified_residue":
            cls = "modified_residue"
        elif single_atom(c.get("formula")):
            cls = "ion"
        elif i["comp_id"] in COMMON:
            cls = "common"
        elif "saccharide" in typ or i["kind"] == "branched":
            glyc = (i["pdb"], i["instance_id"]) in bonded_asn or (i["kind"] == "branched" and
                    ent_has_nag.get((i["pdb"], i["label_asym"]), set()) & {"NAG", "NDG"} and
                    (i["pdb"], i["label_asym"]) in asn_groups)
            cls = "n_glycan" if glyc else "sugar"
        else:
            cls = "other"
        i["ligand_class"] = cls
    # write
    if DB.exists():
        DB.unlink()
    con = sqlite3.connect(DB)
    def table(name, rows, pk=None):
        cols = list(rows[0].keys())
        con.execute(f"CREATE TABLE {name} ({', '.join(c + (' PRIMARY KEY' if c == pk else '') for c in cols)})")
        con.executemany(f"INSERT INTO {name} VALUES ({', '.join('?' * len(cols))})", [tuple(r[c] for c in cols) for r in rows])
    ent_rows = [dict(info[p], n_waters=water.get(p)) for p in ids if p not in miss]
    table("entry", ent_rows, "pdb")
    table("component", list(cc.values()), "comp_id")
    table("instance", inst, "instance_id")
    # de-duplicate: gemmi's st.connections can list the same deposited bond twice (one row per pair of
    # alternate locations, which this schema does not track), giving exact-duplicate rows; found 2026-11-23
    # (docs/04_curation_debt.md item 6), 24 duplicated rows across 14 entries, keyed on full row content
    seen = set(); dedup_links = []
    for r in links:
        k = tuple(r[c] for c in ("pdb", "conn_type", "chain_1", "seq_1", "comp_1", "atom_1",
                                  "chain_2", "seq_2", "comp_2", "atom_2"))
        if k in seen:
            continue
        seen.add(k); dedup_links.append(r)
    print(f"bond rows: {len(links)} before de-duplication, {len(dedup_links)} after "
          f"({len(links) - len(dedup_links)} exact duplicates removed)")
    table("bond", dedup_links)
    table("contact", cont)
    for t, f in (("part1_activity", R / "part1_starch_breaking_activities.tsv"),
                 ("part1_excluded_activity", R / "part1_excluded_activities.tsv"),
                 ("part2_structure", R / "part2_structures_all.tsv"),
                 ("part2_pdb_detail", R / "part2_pdb_details.tsv"),
                 ("part3_catalytic_residue", R / "part3" / "part3_step1_catalytic_residues.tsv"),
                 ("part3_site_geometry", R / "part3" / "part3_site_geometry.tsv"),
                 ("census_linkage", R / "census" / "bound_linkages_checked.tsv"),
                 ("census_on_maps", R / "census" / "census_on_maps.tsv"),
                 ("ligand_entity", R / "ligands" / "ligands_all.tsv")):
        import_tsv(con, t, f)
    for ix in ("CREATE INDEX i_inst_pdb ON instance(pdb)", "CREATE INDEX i_inst_comp ON instance(comp_id)",
               "CREATE INDEX i_cont_inst ON contact(instance_id)", "CREATE INDEX i_bond_pdb ON bond(pdb)"):
        con.execute(ix)
    con.commit()
    q = lambda s: con.execute(s).fetchall()
    summ = dict(source="RCSB mmCIF and wwPDB validation reports, RCSB Data API; built 2026-11-23",
                entries=len(ent_rows), entries_missing=miss, entries_contacts_without_symmetry=nosym, instances=len(inst), bonds=len(links), contacts=len(cont),
                components=len(cc), waters=sum(water.values()),
                instances_by_class=dict(q("SELECT ligand_class, COUNT(*) FROM instance GROUP BY 1 ORDER BY 2 DESC")),
                entries_by_class=dict(q("SELECT ligand_class, COUNT(DISTINCT pdb) FROM instance GROUP BY 1 ORDER BY 2 DESC")),
                bonds_by_type=dict(q("SELECT conn_type, COUNT(*) FROM bond GROUP BY 1 ORDER BY 2 DESC")),
                instances_with_rscc=q("SELECT COUNT(*) FROM instance WHERE rscc IS NOT NULL")[0][0])
    con.close()
    (OUT / "db_summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
