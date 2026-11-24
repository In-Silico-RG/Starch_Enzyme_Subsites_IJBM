"""Census of the alpha(1->4) glycosidic linkages bound in the starch-breaking enzymes (AFC 2026-11-11).

Problem: what phi/psi do real alpha(1->4) linkages take when they are bound in these enzymes, and does any
reach the anti-psi conformation that experiment cannot measure in solution?
Data: in-scope PDB entries of the Part 2 inventory (option C: glycogen phosphorylases excluded) that carry a
bound oligosaccharide of >= 2 residues containing glucose (data/ligands/ligands_all.tsv; 276 entries); coordinates as mmCIF from RCSB, cached in
../data/census/cif/ (not in git).
Method, per entry, per model 1, alternate location A (or blank):
  1. glucose rings = residues GLC or BGC with atoms O5 C1 C2 C3 C4 C5;
  2. an alpha(1->4) linkage = C1 of ring i within 1.75 Å of O4 of ring j (i != j), ring i being GLC (alpha);
  3. phi = O5(i)-C1(i)-O4(j)-C4(j), psi = C1(i)-O4(j)-C4(j)-C5(j)  (same definitions as Enlace_Alfa14_Conformacion);
  4. ring puckering (Cremer-Pople theta) of both rings;
  5. anti-psi criterion of Enlace_Alfa14_Conformacion (docs/05_outline.md): |psi + 117.8| > 90, angle difference taken on
     the circle.
Output: ../data/census/bound_linkages.tsv (one row per linkage) and a summary on stdout.
"""
import csv
import math
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from part3_measure import cremer_pople, dihedral  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent / "resultados"
OUT = ROOT / "census"
CIF = OUT / "cif"


def entries():
    det = {r["pdb"]: r for r in csv.DictReader(open(ROOT / "part2_pdb_details.tsv"), delimiter="\t")}
    rows = [r for r in csv.DictReader(open(ROOT / "part2_structures_all.tsv"), delimiter="\t") if r["in_scope"] == "yes"]
    gp = {r["pdb"] for r in rows if r["family"] == "GT35" and "glycogen" in r["protein"].lower()}
    # entries with a bound oligosaccharide of >= 2 residues that contains glucose (GLC/BGC), from the ligand inventory
    # (scripts/ligands_inventory.py). Corrected 2026-11-23: the first version used longest_glucan_residues >= 2 of
    # part2_pdb_details.tsv under a loose name rule (294 entries, some with N-glycans only); this set (276 entries)
    # contains every entry that gave a linkage, so the census is unchanged.
    withglc = {r["pdb"] for r in csv.DictReader(open(ROOT / "ligands" / "ligands_all.tsv"), delimiter="\t")
               if r["kind"] == "branched" and int(r["n_residues"]) >= 2 and set(r["comp_ids"].split("-")) & {"GLC", "BGC"}}
    info = {}
    for r in rows:
        p = r["pdb"]
        if p in gp or p not in det or p not in withglc:
            continue
        info.setdefault(p, dict(family=r["family"], EC=r["EC"], protein=r["protein"], organism=r["organism"],
                                resolution=det[p]["resolution_A"], method=det[p]["method"]))
    return info


def fetch(p):
    CIF.mkdir(parents=True, exist_ok=True)
    f = CIF / f"{p}.cif"
    if not f.exists():
        for a in range(4):
            try:
                with urllib.request.urlopen(f"https://files.rcsb.org/download/{p}.cif", timeout=180) as r:
                    f.write_bytes(r.read())
                break
            except Exception:
                time.sleep(5 * (a + 1))
    return f


def atoms_from_cif(path):
    cols, rows, inloop = [], [], False
    for line in open(path):
        if line.startswith("_atom_site."):
            cols.append(line.split(".", 1)[1].strip()); inloop = True; continue
        if inloop:
            if line.startswith(("ATOM", "HETATM")):
                rows.append(line.split())
            elif cols and not line.startswith("_atom_site."):
                if rows:
                    break
    ix = {c: i for i, c in enumerate(cols)}
    out = []
    for t in rows:
        if t[ix["pdbx_PDB_model_num"]] != "1" or t[ix["label_alt_id"]] not in (".", "A"):
            continue
        out.append(dict(name=t[ix["label_atom_id"]].strip('"'), resn=t[ix["label_comp_id"]],
                        chain=t[ix["auth_asym_id"]], resi=t[ix["auth_seq_id"]],
                        x=np.array([float(t[ix["Cartn_x"]]), float(t[ix["Cartn_y"]]), float(t[ix["Cartn_z"]])])))
    return out


def circ(a):
    return (a + 180) % 360 - 180


def main():
    info = entries()
    print(len(info), "entries")
    out = []
    for n, (p, meta) in enumerate(sorted(info.items()), 1):
        f = fetch(p)
        if not f.exists():
            print("missing", p); continue
        rings = {}
        for a in atoms_from_cif(f):
            if a["resn"] in ("GLC", "BGC"):
                rings.setdefault((a["chain"], a["resi"], a["resn"]), {})[a["name"]] = a["x"]
        rings = {k: v for k, v in rings.items() if all(x in v for x in ("O5", "C1", "C2", "C3", "C4", "C5", "O4"))}
        for ki, A in rings.items():
            if ki[2] != "GLC":
                continue
            for kj, B in rings.items():
                if kj == ki or float(np.linalg.norm(A["C1"] - B["O4"])) > 1.75:
                    continue
                phi = dihedral(A["O5"], A["C1"], B["O4"], B["C4"])
                psi = dihedral(A["C1"], B["O4"], B["C4"], B["C5"])
                thi = cremer_pople([A[k] for k in ("O5", "C1", "C2", "C3", "C4", "C5")])[1]
                thj = cremer_pople([B[k] for k in ("O5", "C1", "C2", "C3", "C4", "C5")])[1]
                anti = abs(circ(psi + 117.8)) > 90
                out.append([p, meta["family"], meta["EC"], meta["protein"][:60], meta["organism"][:50], meta["method"],
                            meta["resolution"], f"{ki[0]}:{ki[2]}{ki[1]}", f"{kj[0]}:{kj[2]}{kj[1]}",
                            round(phi, 1), round(psi, 1), round(thi, 1), round(thj, 1), "anti-psi" if anti else "normal"])
        if n % 25 == 0:
            print(n, "/", len(info), len(out), "linkages", flush=True)
    with open(OUT / "bound_linkages.tsv", "w") as fh:
        fh.write("pdb\tfamily\tEC\tprotein\torganism\tmethod\tresolution_A\tnonreducing_ring(C1)\treducing_ring(O4)\t"
                 "phi_deg\tpsi_deg\ttheta_nonred_deg\ttheta_red_deg\tconformation\n")
        for r in out:
            fh.write("\t".join(map(str, r)) + "\n")
    na = sum(1 for r in out if r[-1] == "anti-psi")
    print(len(out), "linkages in", len({r[0] for r in out}), "entries;", na, "anti-psi")


if __name__ == "__main__":
    main()
