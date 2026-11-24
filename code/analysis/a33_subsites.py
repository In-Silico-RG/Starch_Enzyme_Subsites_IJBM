"""3.3 Subsite maps and interaction fingerprints.

For each family whose frame carries a catalytic pair (GH13, GH14, GH15, GH31, GH57):
  1. units = sugar rings of the sugar and 'other' instances (one unit per ring; a residue holding two rings,
     such as AC1 = cyclitol + 4-amino-4,6-dideoxyglucose, gives two units); reference atom of a unit = its ring
     atom named C1 / C1<letter> (the anomeric or pseudo-anomeric carbon);
  2. links between units = the deposited covalent bonds between instances (table bond, conn_type Covale) and the
     inter-ring links inside one residue; the donor (non-reducing) unit of a link is the one whose reference atom,
     or an exocyclic atom bonded to it, is the bonded atom;
  3. per connected set of units: subsite -1 = the unit whose reference atom lies nearest the midpoint of the two
     catalytic carboxylate centroids (any accepted chain of the entry), if that distance is <= D_MAX;
  4. labels walk along the main-chain links (acceptor atom O4/N4/C4/S4): toward the reducing end +1, +2, ...;
     toward the non-reducing end -2, -3, ...; units reached through another acceptor atom (O6, O3, O2, O1)
     are labelled "<parent>~O6" etc. and not walked further;
  5. contacts: database table contact (protein partners, 4.0 A heavy atoms), ligand atom -> unit; each contact
     re-measured in the deposited asymmetric unit (contacts with symmetry copies dropped); UniProt number from
     contact_uniprot; reference-structure equivalent = reference Calpha within 2.0 A of the residue's Calpha
     after the catalytic chain's superposition (data/calc/03_frames/superposition.tsv).
Database opened read-only. Run:  nice -n 19 ~/miniconda3/bin/python code/analysis/a33_subsites.py
"""
import sqlite3, re, json, collections
from pathlib import Path
import numpy as np
import pandas as pd
import gemmi
import networkx as nx

P = Path(__file__).resolve().parents[2]
DB = P / 'data/db/enzimas_experimental.sqlite'
FR = P / 'data/calc/03_frames'
CIF = P / 'data/census/cif'
OUT = P / 'data/analysis/33_subsites'
OUT.mkdir(parents=True, exist_ok=True)
FAMS = ['GH13', 'GH14', 'GH15', 'GH31', 'GH57']
D_MAX = 4.5      # A, reference atom of the -1 unit to the carboxylate midpoint (set a priori, see README)
RSCC_MIN = 0.80
MAP_MAX = 2.0    # A, residue Calpha to reference Calpha in the frame
CONSERVED = 0.80
REF_RE = re.compile(r'^C1[A-Z]?$')
MAIN_ACC = re.compile(r'^[ONCS]4[A-Z]?$')


def read_struct(pdb):
    f = CIF / f'{pdb}.cif'
    if not f.exists():
        f = FR / 'ref_cif' / f'{pdb}.cif'
    st = gemmi.read_structure(str(f))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    return st


def residue_index(model):
    idx = {}
    for ch in model:
        for r in ch:
            idx[(ch.name, r.seqid.num, r.seqid.icode.strip())] = r
    return idx


def ring_units(res):
    """rings of one residue -> list of dict(atoms=set(names), ref=name or None); atom -> unit indices."""
    atoms = [a for a in res]
    G = nx.Graph()
    for a in atoms:
        G.add_node(a.name)
    for i, a in enumerate(atoms):
        for b in atoms[i + 1:]:
            d = a.pos.dist(b.pos)
            lim = 1.95 if ('S' in (a.element.name, b.element.name)) else 1.75
            if d < lim:
                G.add_edge(a.name, b.name)
    rings = [set(c) for c in nx.minimum_cycle_basis(G) if 5 <= len(c) <= 7]
    units = []
    for rg in rings:
        refs = sorted(n for n in rg if REF_RE.match(n))
        units.append(dict(atoms=rg, ref=refs[0] if refs else None))
    if not units:  # acyclic (e.g. xylitol) -> one unit
        units.append(dict(atoms=set(G.nodes), ref='C1' if 'C1' in G.nodes else None))
    # atom -> unit(s) by graph distance
    a2u = {}
    for n in G.nodes:
        best, who = 1e9, []
        for k, u in enumerate(units):
            if n in u['atoms']:
                dd = 0
            else:
                try:
                    dd = min(nx.shortest_path_length(G, n, m) for m in u['atoms'])
                except (nx.NetworkXNoPath, ValueError):
                    dd = 1e9
            if dd < best:
                best, who = dd, [k]
            elif dd == best:
                who.append(k)
        a2u[n] = who if best < 1e9 else []
    # inter-ring links inside the residue
    links = []
    for i in range(len(units)):
        for j in range(i + 1, len(units)):
            best = None
            for a in units[i]['atoms']:
                for b in units[j]['atoms']:
                    try:
                        p = nx.shortest_path(G, a, b)
                    except nx.NetworkXNoPath:
                        continue
                    if best is None or len(p) < len(best):
                        best = p
            if best is None or len(best) - 1 > 3:
                continue
            a, b, lk = best[0], best[-1], best[1:-1]
            if units[i]['ref'] == a:
                links.append((i, j, b if not lk else lk[0]))   # donor i, acceptor j, acceptor-side atom
            elif units[j]['ref'] == b:
                links.append((j, i, a if not lk else lk[-1]))
            else:
                links.append((None, (i, j), None))
    return units, a2u, links, G


def main():
    con = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
    entry = pd.read_sql('select pdb, family, uniprot, protein, resolution_A from entry', con).set_index('pdb')
    inst = pd.read_sql("""select i.instance_id, i.pdb, i.auth_asym, i.auth_seq, i.ins_code, i.comp_id, i.rscc,
                                 i.ligand_class, e.family from instance i join entry e using(pdb)
                          where i.ligand_class in ('sugar','other') and e.family in (%s)""" % ','.join(f"'{f}'" for f in FAMS), con)
    inst['ins_code'] = inst['ins_code'].fillna('').astype(str).str.strip()
    bonds = pd.read_sql("""select pdb, instance_id_1, atom_1, instance_id_2, atom_2, distance_A from bond
                           where conn_type='Covale' and instance_id_1 is not null and instance_id_2 is not null""", con)
    contacts = pd.read_sql("""select c.rowid contact_rowid, c.instance_id, c.chain, c.seq, c.ins_code, c.comp,
                                     c.min_distance_A, c.atom, c.partner_atom, u.uniprot, u.uniprot_resnum
                              from contact c left join contact_uniprot u on u.contact_rowid=c.rowid
                              where c.partner_kind='protein'""", con)
    contacts['ins_code'] = contacts['ins_code'].fillna('').astype(str).str.strip()
    sifts = pd.read_sql('select pdb, chain, auth_seq, ins_code, uniprot, uniprot_resnum from sifts_residue', con)
    sifts['ins_code'] = sifts['ins_code'].fillna('').astype(str).str.strip()
    sifts_map = {(r.pdb, r.chain, int(r.auth_seq), r.ins_code): (r.uniprot, r.uniprot_resnum) for r in sifts.itertuples()}

    sup = pd.read_csv(FR / 'superposition.tsv', sep='\t')
    sup = sup[sup.status.str.startswith('accepted')]
    mats = {}
    for r in sup.itertuples():
        v = [float(x) for x in str(r[-1]).split()]
        mats[(r.family, r.pdb, r.chain)] = (np.array(v[:9]).reshape(3, 3), np.array(v[9:12]))
    cat = pd.read_csv(FR / 'catalytic_residues_mapped.tsv', sep='\t', dtype={'ins_code': str})
    cat['ins_code'] = cat['ins_code'].fillna('').astype(str).str.strip()
    cat = cat[cat.family.isin(FAMS) & (cat.type_check == 'same type')]
    fams = {d['family']: d for d in json.load(open(FR / 'families.json'))}
    centres = pd.read_csv(FR / 'ligand_centres_in_frame.tsv', sep='\t')
    assigned_chain = dict(zip(centres.instance_id, centres.assigned_protein_chain))

    # reference structures: Calpha in the frame (reference chain is the frame itself) with UniProt numbers
    refs = {}
    for f in FAMS:
        rp, rc = fams[f]['reference'], fams[f]['reference_chain']
        st = read_struct(rp)
        cas = []
        for r in st[0][rc]:
            ca = r.find_atom('CA', '*')
            if ca is not None and gemmi.find_tabulated_residue(r.name) and gemmi.find_tabulated_residue(r.name).is_amino_acid():
                un = sifts_map.get((rp, rc, r.seqid.num, r.seqid.icode.strip()), (None, None))
                cas.append((r.name, r.seqid.num, r.seqid.icode.strip(), un[1], np.array(ca.pos.tolist())))
        refs[f] = dict(pdb=rp, chain=rc, res=cas, xyz=np.array([c[4] for c in cas]))

    unit_rows, contact_rows, link_log, comp_rows = [], [], [], []
    for (fam, pdb), grp in inst.groupby(['family', 'pdb']):
        st = read_struct(pdb)
        model = st[0]
        ridx = residue_index(model)
        # midpoints of the catalytic carboxylates per accepted chain
        mids = {}
        for ch, cg in cat[cat.pdb == pdb].groupby('chain'):
            if (fam, pdb, ch) not in mats:
                continue
            cents = {}
            for r in cg.itertuples():
                res = ridx.get((ch, int(r.auth_seq), r.ins_code))
                if res is None:
                    continue
                names = ('OD1', 'OD2') if res.name == 'ASP' else ('OE1', 'OE2') if res.name == 'GLU' else ()
                xs = [np.array(a.pos.tolist()) for a in res if a.name in names]
                if len(xs) == 2:
                    cents[r.role] = np.mean(xs, axis=0)
            if 'nuc_or_base' in cents and 'acid' in cents:
                mids[ch] = 0.5 * (cents['nuc_or_base'] + cents['acid'])
        # units
        U = {}          # (instance_id, k) -> dict
        a2u_all = {}
        Gu = nx.Graph(); Du = []   # directed links (donor, acceptor, acceptor_atom)
        for r in grp.itertuples():
            res = ridx.get((r.auth_asym, int(r.auth_seq), r.ins_code))
            if res is None or res.name != r.comp_id:
                link_log.append(dict(pdb=pdb, issue='instance not found in mmCIF', instance_id=r.instance_id, comp=r.comp_id))
                continue
            units, a2u, links, G = ring_units(res)
            a2u_all[r.instance_id] = (a2u, units, G)
            for k, u in enumerate(units):
                ref = u['ref']
                xyz = np.array(res.find_atom(ref, '*').pos.tolist()) if ref else None
                U[(r.instance_id, k)] = dict(instance_id=r.instance_id, ring=k, comp=r.comp_id, ref=ref, xyz=xyz,
                                             rscc=r.rscc, ligand_class=r.ligand_class, chain=r.auth_asym,
                                             seq=r.auth_seq, n_rings=len(units))
                Gu.add_node((r.instance_id, k))
            for d, a, acc_atom in links:
                if d is None:
                    link_log.append(dict(pdb=pdb, issue='intra-residue link, direction unknown', instance_id=r.instance_id, comp=r.comp_id))
                    Gu.add_edge((r.instance_id, a[0]), (r.instance_id, a[1]))
                    continue
                Gu.add_edge((r.instance_id, d), (r.instance_id, a))
                Du.append(((r.instance_id, d), (r.instance_id, a), acc_atom))
        ids = set(grp.instance_id)
        for b in bonds[(bonds.pdb == pdb) & bonds.instance_id_1.isin(ids) & bonds.instance_id_2.isin(ids)].itertuples():
            i1, i2 = int(b.instance_id_1), int(b.instance_id_2)
            if i1 not in a2u_all or i2 not in a2u_all:
                continue
            (m1, un1, G1), (m2, un2, G2) = a2u_all[i1], a2u_all[i2]
            u1s, u2s = m1.get(b.atom_1, []), m2.get(b.atom_2, [])
            if not u1s or not u2s:
                link_log.append(dict(pdb=pdb, issue=f'bond atom not in residue {b.atom_1}/{b.atom_2}', instance_id=i1, comp=''))
                continue
            u1, u2 = u1s[0], u2s[0]
            def is_donor(atom, units, uk, G):
                ref = units[uk]['ref']
                if ref is None:
                    return False
                if atom == ref:
                    return True
                return atom not in units[uk]['atoms'] and G.has_edge(atom, ref)
            n1, n2 = (i1, u1), (i2, u2)
            Gu.add_edge(n1, n2)
            if b.atom_1 == un1[u1]['ref'] or (b.atom_2 != un2[u2]['ref'] and is_donor(b.atom_1, un1, u1, G1)):
                Du.append((n1, n2, b.atom_2))
            elif b.atom_2 == un2[u2]['ref'] or is_donor(b.atom_2, un2, u2, G2):
                Du.append((n2, n1, b.atom_1))
            else:
                link_log.append(dict(pdb=pdb, issue=f'deposited bond, direction unknown {b.atom_1}-{b.atom_2}', instance_id=i1, comp=''))
        # per connected set
        for ci, comp in enumerate(nx.connected_components(Gu)):
            best = (1e9, None, None)
            for n in comp:
                if U[n]['xyz'] is None:
                    continue
                for ch, m in mids.items():
                    d = float(np.linalg.norm(U[n]['xyz'] - m))
                    if d < best[0]:
                        best = (d, n, ch)
            label = {}
            if best[1] is not None and best[0] <= D_MAX:
                label[best[1]] = -1
                # walk
                out = collections.defaultdict(list)
                for d_, a_, acc in Du:
                    if d_ in comp:
                        out[d_].append((a_, 'toward_reducing', acc))
                        out[a_].append((d_, 'toward_nonreducing', acc))
                queue = [best[1]]
                while queue:
                    u = queue.pop(0)
                    p = label[u]
                    if isinstance(p, str):
                        continue
                    for v, direction, acc in out[u]:
                        if MAIN_ACC.match(acc or ''):
                            q = (p + 1 if p != -1 else 1) if direction == 'toward_reducing' else (p - 1 if p != 1 else -1)
                        else:
                            q = f'{fmt(p)}~{acc}'
                        if v in label:
                            if label[v] != q:
                                link_log.append(dict(pdb=pdb, issue=f'label conflict {fmt(label[v])} vs {fmt(q)}', instance_id=v[0], comp=U[v]['comp']))
                            continue
                        label[v] = q
                        queue.append(v)
            cat_chain = best[2]
            comp_rows.append(dict(family=fam, pdb=pdb, component=ci, n_units=len(comp), d_min_A=round(best[0], 2) if best[1] else None,
                                  nearest_unit=f'{U[best[1]]["comp"]}:{U[best[1]]["chain"]}{U[best[1]]["seq"]}:{U[best[1]]["ref"]}' if best[1] else None,
                                  catalytic_chain=cat_chain, assigned=best[1] is not None and best[0] <= D_MAX))
            for n in comp:
                u = U[n]
                ch_for_frame = cat_chain if label else assigned_chain.get(u['instance_id'])
                fxyz = [None] * 3
                if u['xyz'] is not None and ch_for_frame and (fam, pdb, ch_for_frame) in mats:
                    R, t = mats[(fam, pdb, ch_for_frame)]
                    fxyz = list(np.round(R @ u['xyz'] + t, 3))
                dmid = min((float(np.linalg.norm(u['xyz'] - m)) for m in mids.values()), default=None) if u['xyz'] is not None else None
                unit_rows.append(dict(family=fam, pdb=pdb, component=ci, instance_id=u['instance_id'], ring=u['ring'],
                                      n_rings_in_residue=u['n_rings'], comp=u['comp'], chain=u['chain'], seq=u['seq'],
                                      ref_atom=u['ref'], ligand_class=u['ligand_class'], rscc=u['rscc'],
                                      subsite=fmt(label[n]) if n in label else '', frame_chain=ch_for_frame,
                                      d_ref_to_midpoint_A=round(dmid, 2) if dmid is not None else None,
                                      x_frame=fxyz[0], y_frame=fxyz[1], z_frame=fxyz[2]))
        # contacts of this entry's instances
        gidx = grp.set_index('instance_id')
        for c in contacts[contacts.instance_id.isin(ids)].itertuples():
            if c.instance_id not in a2u_all:
                continue
            m, units, G = a2u_all[c.instance_id]
            gi = gidx.loc[c.instance_id]
            lig_res = ridx[(gi['auth_asym'], int(gi['auth_seq']), gi['ins_code'])]
            pr = ridx.get((c.chain, int(c.seq), c.ins_code))
            la = lig_res.find_atom(c.atom, '*') if c.atom else None
            dasu = None
            if pr is not None and la is not None:
                dasu = min(la.pos.dist(a.pos) for a in pr)
            ca = pr.find_atom('CA', '*') if pr is not None else None
            contact_rows.append(dict(family=fam, pdb=pdb, instance_id=c.instance_id, ligand_atom=c.atom,
                                     rings=','.join(str(k) for k in m.get(c.atom, [])), chain=c.chain, seq=c.seq,
                                     ins_code=c.ins_code, resname=c.comp, uniprot=c.uniprot, uniprot_resnum=c.uniprot_resnum,
                                     min_distance_db_A=c.min_distance_A, distance_asu_A=round(dasu, 2) if dasu is not None else None,
                                     ca=np.array(ca.pos.tolist()) if ca is not None else None))
        print(pdb, fam, len(grp), flush=True)

    units = pd.DataFrame(unit_rows)
    comps = pd.DataFrame(comp_rows)
    cont = pd.DataFrame(contact_rows)
    logs = pd.DataFrame(link_log)
    comps.to_csv(OUT / 'connected_sets.tsv', sep='\t', index=False)
    logs.to_csv(OUT / 'link_issues.tsv', sep='\t', index=False)

    # contacts: drop symmetry contacts, map to reference positions in the frame of the unit's chain
    cont['asu'] = cont['distance_asu_A'].notna() & (cont['distance_asu_A'] <= 4.05)
    ufr = units.drop_duplicates('instance_id').set_index('instance_id')['frame_chain']
    ref_eq, ref_un, ref_d = [], [], []
    for c in cont.itertuples():
        fc = ufr.get(c.instance_id)
        key = (c.family, c.pdb, fc)
        if not c.asu or c.ca is None or key not in mats:
            ref_eq.append(''); ref_un.append(None); ref_d.append(None); continue
        R, t = mats[key]
        x = R @ c.ca + t
        dd = np.linalg.norm(refs[c.family]['xyz'] - x, axis=1)
        k = int(dd.argmin())
        if dd[k] <= MAP_MAX:
            rr = refs[c.family]['res'][k]
            ref_eq.append(f'{gemmi.find_tabulated_residue(rr[0]).one_letter_code.upper()}{rr[1]}{rr[2]}')
            ref_un.append(rr[3]); ref_d.append(round(float(dd[k]), 2))
        else:
            ref_eq.append(''); ref_un.append(None); ref_d.append(round(float(dd[k]), 2))
    cont['ref_position'] = ref_eq
    cont['ref_uniprot_resnum'] = ref_un
    cont['ref_ca_distance_A'] = ref_d
    cont['same_chain_as_frame'] = [c.chain == ufr.get(c.instance_id) for c in cont.itertuples()]
    cont.drop(columns=['ca']).to_csv(OUT / 'contacts_mapped.tsv', sep='\t', index=False)
    units.to_csv(OUT / 'units.tsv', sep='\t', index=False)

    # unit-level contact sets
    # units reported in the residue tables: sugar class, or 'other' class in a connected set that holds a sugar unit
    # (pseudo-sugar units bonded into an oligosaccharide, e.g. the valienamine HMC/ACI of acarbose analogues)
    has_sugar = units.groupby(['pdb', 'component'])['ligand_class'].apply(lambda s: (s == 'sugar').any())
    units['set_has_sugar'] = [has_sugar[(a, b)] for a, b in zip(units.pdb, units.component)]
    units['set_n_instances'] = units.groupby(['pdb', 'component'])['instance_id'].transform('nunique')
    units['in_tables'] = (units.ligand_class == 'sugar') | ((units.ligand_class == 'other') & units.set_has_sugar & (units.set_n_instances > 1))
    units.to_csv(OUT / 'units.tsv', sep='\t', index=False)
    rep = units[(units.subsite != '') & (units.rscc >= RSCC_MIN) & units.in_tables].copy()
    good = cont[cont.asu]
    def unit_contacts(iid, ring):
        g = good[good.instance_id == iid]
        g = g[g.rings.apply(lambda s: str(ring) in s.split(','))]
        return g
    occ_rows, prot_rows = [], []
    sub_entries = collections.defaultdict(set)
    res_entries = collections.defaultdict(set)
    ref_un_of = {}
    prot_sub_entries = collections.defaultdict(set)
    prot_res_entries = collections.defaultdict(set)
    for u in rep.itertuples():
        sub_entries[(u.family, u.subsite)].add(u.pdb)
        acc = entry.loc[u.pdb, 'uniprot']
        prot_sub_entries[(u.family, acc, u.subsite)].add(u.pdb)
        g = unit_contacts(u.instance_id, u.ring)
        for c in g.itertuples():
            if c.ref_position and c.same_chain_as_frame:
                res_entries[(u.family, u.subsite, c.ref_position)].add(u.pdb)
                ref_un_of[(u.family, c.ref_position)] = c.ref_uniprot_resnum
            if pd.notna(c.uniprot_resnum) and isinstance(c.uniprot, str):
                prot_res_entries[(u.family, c.uniprot, u.subsite, c.resname, int(c.uniprot_resnum))].add(u.pdb)
    for (f, s, rp), es in res_entries.items():
        run = ref_un_of.get((f, rp))
        n = len(sub_entries[(f, s)])
        occ_rows.append(dict(family=f, subsite=s, ref_position=rp, ref_uniprot_resnum=run,
                             n_entries_subsite=n, n_entries_contact=len(es), fraction=round(len(es) / n, 3),
                             category='conserved' if len(es) / n >= CONSERVED else 'variable'))
    occ = pd.DataFrame(occ_rows)
    occ['subsite_order'] = occ.subsite.map(order_key)
    occ = occ.sort_values(['family', 'subsite_order', 'fraction'], ascending=[True, True, False]).drop(columns='subsite_order')
    occ.to_csv(OUT / 'subsite_residues_reference_numbering.tsv', sep='\t', index=False)
    for (f, acc, s, rn, num), es in prot_res_entries.items():
        n = len(prot_sub_entries.get((f, acc, s), set())) or None
        prot_rows.append(dict(family=f, uniprot=acc, subsite=s, resname=rn, uniprot_resnum=num,
                              n_entries_subsite=n, n_entries_contact=len(es)))
    pd.DataFrame(prot_rows).sort_values(['family', 'uniprot', 'subsite', 'uniprot_resnum']).to_csv(
        OUT / 'subsite_residues_uniprot_numbering.tsv', sep='\t', index=False)

    # per family per subsite summary
    srows = []
    for (f, s), es in sub_entries.items():
        d = rep[(rep.family == f) & (rep.subsite == s)]
        xyz = d[['x_frame', 'y_frame', 'z_frame']].astype(float).values
        rms = float(np.sqrt(((xyz - xyz.mean(0)) ** 2).sum(1).mean())) if len(xyz) > 1 else 0.0
        comps_ = collections.Counter(d.comp)
        nc = occ[(occ.family == f) & (occ.subsite == s)]
        srows.append(dict(family=f, subsite=s, n_units=len(d), n_entries=len(es), n_proteins=d.pdb.map(entry.uniprot).nunique(),
                          components=';'.join(f'{k}:{v}' for k, v in comps_.most_common()),
                          frame_rms_radius_A=round(rms, 2),
                          frame_centroid=' '.join(f'{v:.2f}' for v in xyz.mean(0)),
                          n_residues_conserved=int((nc.category == 'conserved').sum()), n_residues_variable=int((nc.category == 'variable').sum()),
                          conserved=','.join(nc[nc.category == 'conserved'].ref_position)))
    ss = pd.DataFrame(srows)
    ss['o'] = ss.subsite.map(order_key)
    ss = ss.sort_values(['family', 'o']).drop(columns='o')
    ss.to_csv(OUT / 'subsite_summary.tsv', sep='\t', index=False)

    # fingerprints: per ligand instance (RSCC >= 0.80), bits over the family's site residues
    for f in FAMS:
        site = occ[occ.family == f][['ref_position', 'ref_uniprot_resnum']].drop_duplicates('ref_position')
        site['n'] = site.ref_position.str.extract(r'(\d+)').astype(int)
        site = site.sort_values('n')
        cols = list(site.ref_position)
        if not cols:
            continue
        ii = units[(units.family == f) & (units.rscc >= RSCC_MIN)].groupby('instance_id').agg(
            pdb=('pdb', 'first'), comp=('comp', 'first'), ligand_class=('ligand_class', 'first'), rscc=('rscc', 'first'),
            subsite=('subsite', lambda s: ','.join(x for x in s if x))).reset_index()
        g = good[(good.family == f) & good.same_chain_as_frame & (good.ref_position != '')]
        bits = g.groupby('instance_id')['ref_position'].apply(set)
        M = np.array([[int(c in bits.get(i, set())) for c in cols] for i in ii.instance_id], dtype=int)
        fp = pd.concat([ii, pd.DataFrame(M, columns=cols)], axis=1)
        fp['n_bits'] = M.sum(1)
        fp.to_csv(OUT / f'fingerprints_{f}.tsv', sep='\t', index=False)

    # figures
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for f in ['GH13', 'GH14', 'GH15', 'GH31', 'GH57']:
        d = occ[(occ.family == f) & (~occ.subsite.str.contains('~'))]
        if d.empty:
            continue
        keep = d.groupby('ref_position')['fraction'].max()
        keep = keep[keep >= 0.25].index
        d = d[d.ref_position.isin(keep)]
        subs = sorted(d.subsite.unique(), key=order_key)
        resn = sorted(d.ref_position.unique(), key=lambda s: int(re.search(r'\d+', s).group()))
        M = np.full((len(subs), len(resn)), np.nan)
        for r in d.itertuples():
            M[subs.index(r.subsite), resn.index(r.ref_position)] = r.fraction
        fig, ax = plt.subplots(figsize=(max(6, 0.32 * len(resn) + 2), 0.45 * len(subs) + 2))
        im = ax.imshow(M, aspect='auto', cmap='viridis', vmin=0, vmax=1)
        ax.set_xticks(range(len(resn))); ax.set_xticklabels(resn, rotation=90, fontsize=7)
        nsub = {s: len(sub_entries[(f, s)]) for s in subs}
        ax.set_yticks(range(len(subs))); ax.set_yticklabels([f'{s}  (n={nsub[s]})' for s in subs], fontsize=8)
        for i in range(len(subs)):
            for j in range(len(resn)):
                if not np.isnan(M[i, j]) and M[i, j] >= CONSERVED:
                    ax.text(j, i, '●', ha='center', va='center', fontsize=5, color='w')
        cb = fig.colorbar(im, ax=ax, fraction=0.03); cb.set_label('fraction of entries occupying the subsite\nwith the residue within 4.0 Å', fontsize=7)
        ax.set_xlabel(f'residue (reference {refs[f]["pdb"]} chain {refs[f]["chain"]} numbering; shown if ≥ 0.25 at some subsite)', fontsize=8)
        ax.set_ylabel('subsite (n = entries)', fontsize=8)
        ax.set_title(f'{f}: residues within 4.0 Å of the sugar unit at each subsite (dot: ≥ 0.80)', fontsize=9)
        fig.tight_layout()
        fig.savefig(OUT / f'subsite_map_{f}.pdf'); fig.savefig(OUT / f'subsite_map_{f}.png', dpi=150)
        plt.close(fig)

    summ = dict(gemmi=gemmi.__version__, networkx=nx.__version__, D_MAX_A=D_MAX, RSCC_MIN=RSCC_MIN, MAP_MAX_A=MAP_MAX,
                n_units=len(units), n_connected_sets=len(comps), n_sets_assigned=int(comps.assigned.sum()),
                n_units_labelled=int((units.subsite != '').sum()),
                n_units_labelled_rscc=int(len(rep)), n_link_issues=len(logs),
                n_contacts=len(cont), n_contacts_symmetry_dropped=int((~cont.asu).sum()),
                n_contacts_mapped_to_reference=int((cont.ref_position != '').sum()))
    json.dump(summ, open(OUT / 'summary.json', 'w'), indent=1)
    print(json.dumps(summ, indent=1))


def fmt(p):
    if isinstance(p, str):
        return p
    return f'+{p}' if p > 0 else f'{p}'


def order_key(s):
    m = re.match(r'^([+-]?\d+)', s)
    base = int(m.group(1)) if m else 99
    return (base, s)


if __name__ == '__main__':
    main()
