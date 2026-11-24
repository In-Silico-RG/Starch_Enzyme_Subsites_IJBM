"""Add ChEMBL targets, activities and molecules to the experimental database (AFC 2026-11-23).

Problem: the structures say where molecules bind, not how strongly. ChEMBL holds potencies measured experimentally
(IC50, Ki, Kd, EC50) curated from the literature. Joining both gives structure, site and potency for the same
proteins.
Data: ChEMBL REST API, https://www.ebi.ac.uk/chembl/api/data/ (target, activity, molecule endpoints). The proteins are
the UniProt accessions mapped by SIFTS to the 645 entries (table sifts_residue). Raw answers cached in
../data/db/chembl_raw/ (not in git).
Method:
  1. for each accession, the ChEMBL targets that contain it (table chembl_target, with target type);
  2. for each SINGLE PROTEIN target, every activity with standard_type IC50, Ki, Kd or EC50 (table chembl_activity),
     keeping relation, value, units, pChEMBL (= -log10 of the molar value, when ChEMBL computes it), assay type
     (B binding/inhibition, F functional, A ADME...), assay description, document and ChEMBL validity comment;
  3. for each molecule, its structure and InChIKey (table chembl_molecule);
  4. table chembl_pdb_ligand: molecules of ChEMBL that are also deposited ligands of the database (same InChIKey as a
     row of table component).
Output: tables chembl_target, chembl_activity, chembl_molecule, chembl_pdb_ligand in
../data/db/enzimas_experimental.sqlite; summary ../data/db/chembl_summary.json.
"""
import hashlib
import json
import sqlite3
import time
import urllib.parse
import urllib.request
from pathlib import Path

R = Path(__file__).resolve().parent.parent / "resultados"
DB = R / "db" / "enzimas_experimental.sqlite"
RAW = R / "db" / "chembl_raw"
API = "https://www.ebi.ac.uk"
TYPES = "IC50,Ki,Kd,EC50"


def get_json(path):
    f = RAW / (hashlib.md5(path.encode()).hexdigest() + ".json")
    if f.exists():
        return json.loads(f.read_text())
    for a in range(6):
        try:
            req = urllib.request.Request(API + path, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=180) as r:
                d = json.load(r)
            f.write_text(json.dumps(d)); return d
        except Exception as ex:  # noqa: BLE001
            print("retry", a, path[:90], ex, flush=True); time.sleep(5 * (a + 1))
    raise RuntimeError(path)


FAILED = []


def paged(path, key):
    out = []
    while path:
        try:
            d = get_json(path)
        except RuntimeError:
            FAILED.append(path); return out
        out += d.get(key, [])
        path = (d.get("page_meta") or {}).get("next")
        time.sleep(0.2)
    return out


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(DB)
    accs = [r[0] for r in con.execute("SELECT DISTINCT uniprot FROM sifts_residue WHERE uniprot IS NOT NULL ORDER BY 1")]
    targets, acts = [], []
    for n, acc in enumerate(accs, 1):
        for t in paged(f"/chembl/api/data/target.json?target_components__accession={acc}&limit=100", "targets"):
            targets.append((acc, t["target_chembl_id"], t["target_type"], t["pref_name"], t["organism"]))
        if n % 25 == 0:
            print(n, "/", len(accs), "accessions", len(targets), "targets", flush=True)
    single = sorted({t[1] for t in targets if t[2] == "SINGLE PROTEIN"})
    for n, tid in enumerate(single, 1):
        for a in paged(f"/chembl/api/data/activity.json?target_chembl_id={tid}&standard_type__in={TYPES}&limit=1000",
                       "activities"):
            acts.append((a["activity_id"], tid, a["molecule_chembl_id"], a.get("canonical_smiles"), a["standard_type"],
                         a.get("standard_relation"), a.get("standard_value"), a.get("standard_units"), a.get("pchembl_value"),
                         a.get("assay_chembl_id"), a.get("assay_type"), a.get("assay_description"), a.get("document_chembl_id"),
                         a.get("document_year"), a.get("data_validity_comment"), a.get("target_organism")))
        print(tid, n, "/", len(single), len(acts), "activities", flush=True)
    mols = sorted({a[2] for a in acts})
    mrows = []
    for i in range(0, len(mols), 50):
        ids = ",".join(mols[i:i + 50])
        for m in paged(f"/chembl/api/data/molecule.json?molecule_chembl_id__in={urllib.parse.quote(ids)}&limit=50", "molecules"):
            s = m.get("molecule_structures") or {}; pr = m.get("molecule_properties") or {}
            mrows.append((m["molecule_chembl_id"], m.get("pref_name"), s.get("canonical_smiles"), s.get("standard_inchi_key"),
                          pr.get("full_mwt"), pr.get("alogp"), pr.get("hbd"), pr.get("hba"), pr.get("num_ro5_violations"),
                          m.get("max_phase")))
    for t in ("chembl_target", "chembl_activity", "chembl_molecule", "chembl_pdb_ligand"):
        con.execute(f"DROP TABLE IF EXISTS {t}")
    con.execute("CREATE TABLE chembl_target (uniprot, target_chembl_id, target_type, pref_name, organism)")
    con.executemany("INSERT INTO chembl_target VALUES (?,?,?,?,?)", targets)
    con.execute("""CREATE TABLE chembl_activity (activity_id, target_chembl_id, molecule_chembl_id, canonical_smiles,
                   standard_type, standard_relation, standard_value, standard_units, pchembl_value, assay_chembl_id,
                   assay_type, assay_description, document_chembl_id, document_year, data_validity_comment, target_organism)""")
    con.executemany(f"INSERT INTO chembl_activity VALUES ({','.join('?' * 16)})", acts)
    con.execute("""CREATE TABLE chembl_molecule (molecule_chembl_id PRIMARY KEY, pref_name, canonical_smiles, inchikey,
                   full_mwt, alogp, hbd, hba, num_ro5_violations, max_phase)""")
    con.executemany(f"INSERT OR IGNORE INTO chembl_molecule VALUES ({','.join('?' * 10)})", mrows)
    con.execute("""CREATE TABLE chembl_pdb_ligand AS SELECT c.comp_id, c.name AS pdb_name, m.molecule_chembl_id, m.pref_name
                   FROM component c JOIN chembl_molecule m ON m.inchikey = c.inchikey WHERE c.inchikey IS NOT NULL""")
    con.commit()
    q = lambda s: con.execute(s).fetchall()
    summ = dict(source="ChEMBL REST API, queried 2026-11-23", uniprot_accessions=len(accs),
                targets_by_type=q("SELECT target_type, COUNT(DISTINCT target_chembl_id) FROM chembl_target GROUP BY 1"),
                activities=len(acts), molecules=len(mrows),
                activities_by_assay_type=q("SELECT assay_type, COUNT(*) FROM chembl_activity GROUP BY 1 ORDER BY 2 DESC"),
                activities_by_type=q("SELECT standard_type, COUNT(*) FROM chembl_activity GROUP BY 1 ORDER BY 2 DESC"),
                with_pchembl=q("SELECT COUNT(*) FROM chembl_activity WHERE pchembl_value IS NOT NULL")[0][0],
                top_targets=q("""SELECT a.target_chembl_id, t.pref_name, t.organism, COUNT(*) FROM chembl_activity a
                                 JOIN (SELECT DISTINCT target_chembl_id, pref_name, organism FROM chembl_target) t USING(target_chembl_id)
                                 GROUP BY 1 ORDER BY 4 DESC LIMIT 20"""),
                pdb_ligands_in_chembl=q("SELECT comp_id, pdb_name, molecule_chembl_id FROM chembl_pdb_ligand"),
                failed_requests=FAILED)
    (R / "db" / "chembl_summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))
    con.close()


if __name__ == "__main__":
    main()
