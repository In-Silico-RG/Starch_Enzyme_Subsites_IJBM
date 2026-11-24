"""Calculation 3: family frames, conserved waters and ions, ligand centres in the frame.

For each family with >= 2 entries:
 1. Reference structure: the Part 3 reference (first row of part3_site_geometry for the family),
    with its catalytic residues (nuc_or_base, acid). Families without a Part 3 reference: the
    best-resolution entry that holds a sugar/'other' ligand; its "site residues" are the protein
    residues within 4.0 A of those ligands (database table contact).
 2. Every protein chain (>= 30 amino acids with CA) of every entry is aligned by sequence to the
    reference chain (Biopython PairwiseAligner, global, BLOSUM62, gap open -10, extend -0.5,
    free end gaps). Core = aligned pairs whose reference CA lies within 12.0 A of any atom of the
    reference catalytic (site) residues.
 3. Kabsch on the CA of the core pairs; iterated: refit on pairs within 2.0 A after the
    previous fit, until stable (max 10 rounds). A chain is accepted when >= 20 pairs are
    retained and they are >= 40 % of the core pairs. RMSD reported on all core pairs and on the
    retained pairs.
 4. Waters (HOH) and ions (database class 'ion') and sugar/'other' instance centres are given to
    the nearest accepted chain (nearest protein heavy atom) and moved with that chain's
    transform into the reference frame.
 5. Waters and ions (separately) are clustered by greedy density clustering, radius 1.0 A
    around the seed (see cluster()); clusters holding members from >= 50 %
    of the family's superposed entries are reported with their distance to the catalytic
    residues and nearest reference residue.
Usage: c03_frames.py [FAMILY ...]   (default: all families with >= 2 entries)
"""
import os, sys, json, collections
import numpy as np
import gemmi
from scipy.spatial import cKDTree
from Bio import Align
from Bio.Align import substitution_matrices
sys.path.insert(0, os.path.dirname(__file__))
from common import load, CALC, db, is_protein_residue, CIF

OUT = os.path.join(CALC, '03_frames')
REFCIF = os.path.join(OUT, 'ref_cif')
CORE_R, TRIM, MIN_PAIRS, MIN_FRAC, CLUST_R, PRESENCE = 12.0, 2.0, 20, 0.40, 1.0, 0.50

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.end_gap_score = 0
aligner.mode = 'global'


def load_any(pdb):
    if os.path.exists(os.path.join(CIF, pdb + '.cif')):
        return load(pdb)
    st = gemmi.read_structure(os.path.join(REFCIF, pdb + '.cif'))
    st.setup_entities(); st.remove_hydrogens(); st.remove_alternative_conformations()
    for i in range(len(st) - 1, 0, -1):
        del st[i]
    return st


def one(res):
    info = gemmi.find_tabulated_residue(res.name)
    c = info.one_letter_code.upper() if info else 'X'
    return c if c.isalpha() else 'X'


def chains_of(st):
    """{chain: [(residue, CA position)]} for protein chains with >= 30 residues having CA."""
    out = {}
    for ch in st[0]:
        rl = []
        for r in ch:
            if is_protein_residue(r):
                ca = r.find_atom('CA', '*')
                if ca:
                    rl.append((r, np.array(ca.pos.tolist())))
        if len(rl) >= 30:
            out[ch.name] = rl
    return out


def kabsch(P, Q):
    """R, t minimising |R P + t - Q|."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def trim_fit(P, Q, R, t, cut, rounds=20):
    keep = None
    for _ in range(rounds):
        d = np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1))
        nk = d <= cut
        if nk.sum() < 3 or (keep is not None and (nk == keep).all()):
            break
        keep = nk
        R, t = kabsch(P[keep], Q[keep])
    return R, t


def icp(Rs, ts, refcore_idx, ref_ca, tgt_ca, cut=3.0, rounds=30):
    """Structural matching: each reference core CA paired with the nearest target CA within `cut`
    (mutual nearest), refit, repeat. Returns R, t, pairs [(i_ref, j_tgt)]."""
    R, t, prev = Rs, ts, None
    for _ in range(rounds):
        T = tgt_ca @ R.T + t
        tr = cKDTree(T)
        d, j = tr.query(ref_ca[refcore_idx])
        back = cKDTree(ref_ca).query(T[j])[1]
        pr = [(i, jj) for i, dd, jj, b in zip(refcore_idx, d, j, back) if dd <= cut and b == i]
        if len(pr) < 3 or pr == prev:
            break
        prev = pr
        R, t = kabsch(tgt_ca[[b for a, b in pr]], ref_ca[[a for a, b in pr]])
    return R, t, prev or []


def rmsd(P, Q):
    return float(np.sqrt(((P - Q) ** 2).sum(1).mean()))


def sifts_map(con, pdb):
    m = {}
    try:
        for ch, seq, ins, un, unum in con.execute(
                "select chain, auth_seq, ins_code, uniprot, uniprot_resnum from sifts_residue where pdb=?", (pdb,)):
            m[(ch, int(seq), (ins or '').strip())] = (un or '', '' if unum is None else str(unum))
    except Exception:
        pass
    return m


def rid(r):
    return (r.seqid.num, r.seqid.icode.strip())


def family_reference(con, fam):
    row = con.execute("select pdb, nuc_or_base, acid from part3_site_geometry where family=? limit 1", (fam,)).fetchone()
    if row:
        pdb, nb, ac = row
        site = []
        for role, s in (('nuc_or_base', nb), ('acid', ac)):
            name, num = s.split(':')[0][:3], int(s.split(':')[0][3:])
            site.append((role, name, num))
        return pdb, site, 'Part 3 reference (part3_site_geometry)'
    # no Part 3 reference: best-resolution entry with a sugar/'other' ligand
    rows = con.execute("""select e.pdb from entry e where e.family=? and exists
            (select 1 from instance i where i.pdb=e.pdb and i.ligand_class in ('sugar','other'))
            order by cast(e.resolution_A as real)""", (fam,)).fetchall()
    if not rows:
        rows = con.execute("select pdb from entry where family=? order by cast(resolution_A as real)", (fam,)).fetchall()
        return rows[0][0], [], 'best-resolution entry; no ligand: core = whole chain'
    pdb = rows[0][0]
    site = con.execute("""select distinct c.chain, c.seq, c.ins_code, c.comp from contact c join instance i
            on c.instance_id=i.instance_id where i.pdb=? and i.ligand_class in ('sugar','other')
            and c.partner_kind='protein'""", (pdb,)).fetchall()
    return pdb, [('site', comp, (ch, seq, ins or '')) for ch, seq, ins, comp in site], \
        'best-resolution entry with sugar/other ligand; site = residues within 4 A of them'


def run_family(fam, con):
    pdbs = [r[0] for r in con.execute("select pdb from entry where family=? order by pdb", (fam,))]
    refpdb, site, refnote = family_reference(con, fam)
    ref = load_any(refpdb)
    rch = chains_of(ref)
    # reference chain and catalytic residues
    cat = []  # (role, chain, num, icode, resname, residue)
    if site and site[0][0] != 'site':
        for chname, rl in rch.items():
            hits = []
            for role, name, num in site:
                for r, ca in rl:
                    if r.seqid.num == num and (r.name == name or True):
                        hits.append((role, chname, num, '', r.name, r))
            if len(hits) == len(site):
                cat = hits
                break
    elif site:
        chain_counts = collections.Counter(s[2][0] for s in site)
        chname = chain_counts.most_common(1)[0][0]
        want = {(s[2][1], s[2][2]) for s in site if s[2][0] == chname}
        for r, ca in rch.get(chname, []):
            if rid(r) in want:
                cat.append(('site', chname, r.seqid.num, r.seqid.icode.strip(), r.name, r))
    refchain = cat[0][1] if cat else sorted(rch, key=lambda c: -len(rch[c]))[0]
    rl_ref = rch[refchain]
    cat_atoms = np.array([a.pos.tolist() for c in cat for a in c[5]]) if cat else None
    ref_ca = np.array([ca for r, ca in rl_ref])
    if cat_atoms is not None:
        dmin = np.sqrt(((ref_ca[:, None, :] - cat_atoms[None]) ** 2).sum(-1)).min(1)
        core_mask = dmin <= CORE_R
    else:
        core_mask = np.ones(len(rl_ref), bool)
    refseq = ''.join(one(r) for r, _ in rl_ref)
    ref_sifts = sifts_map(con, refpdb)
    # functional atoms of catalytic residues (carboxylate O, else side-chain tip, else CA)
    func = []
    for role, ch, num, ic, name, r in cat:
        sel = [a for a in r if a.name in ('OD1', 'OD2', 'OE1', 'OE2', 'NZ', 'NE2', 'ND1', 'OH', 'OG1', 'NH1', 'NH2')]
        if not sel:
            sel = [a for a in r if a.name == 'CA']
        func.append(('%s %s%d%s' % (role, name, num, ic), np.array([a.pos.tolist() for a in sel])))
    ref_atoms, ref_atom_res = [], []
    for r, ca in rl_ref:
        for a in r:
            ref_atoms.append(a.pos.tolist()); ref_atom_res.append(r)
    ref_atoms = np.array(ref_atoms); ref_tree = cKDTree(ref_atoms)

    def where(x):
        """distance to each catalytic residue (functional atoms) and nearest reference residue."""
        out = []
        for lab, A in func:
            out.append('%s:%.2f' % (lab.replace(' ', '_'), np.sqrt(((A - x) ** 2).sum(1)).min()))
        d, i = ref_tree.query(x)
        r = ref_atom_res[i]
        un = ref_sifts.get((refchain, r.seqid.num, r.seqid.icode.strip()), ('', ''))
        return ';'.join(out), '%s:%s%d%s' % (refchain, r.name, r.seqid.num, r.seqid.icode.strip()), un[1], '%.2f' % d

    sup_rows, catmap_rows, wat, ion, lig_rows = [], [], [], [], []
    entries_ok = set()
    for pdb in pdbs:
        try:
            st = ref if pdb == refpdb else load_any(pdb)
        except Exception as e:
            sup_rows.append([fam, pdb, '', refpdb, refchain, 'load failed: %r' % e] + [''] * 7)
            continue
        sm = sifts_map(con, pdb)
        chs = chains_of(st)
        transforms = {}
        for chname, rl in chs.items():
            seq = ''.join(one(r) for r, _ in rl)
            aln = aligner.align(refseq, seq)[0]
            pairs = []
            ident = 0
            for (a0, a1), (b0, b1) in zip(*aln.aligned):
                for k in range(a1 - a0):
                    pairs.append((a0 + k, b0 + k))
                    ident += refseq[a0 + k] == seq[b0 + k]
            pid = ident / max(1, len(pairs))
            core = [(i, j) for i, j in pairs if core_mask[i]]
            if len(core) < MIN_PAIRS:
                sup_rows.append([fam, pdb, chname, refpdb, refchain, 'rejected: %d core pairs' % len(core),
                                 '%.3f' % pid, len(core), '', '', '', '', ''])
                continue
            P = np.array([rl[j][1] for i, j in core]); Q = np.array([ref_ca[i] for i, j in core])
            PA = np.array([rl[j][1] for i, j in pairs]); QA = np.array([ref_ca[i] for i, j in pairs])
            R0, t0 = kabsch(P, Q)
            rms_all0 = rmsd(P @ R0.T + t0, Q)
            best = None
            # starting fits: (a) the core pairs, (b) all aligned pairs of the chain, trimmed at 3 A
            RA, tA = trim_fit(PA, QA, *kabsch(PA, QA), 3.0)
            for Rs, ts in ((R0, t0), (RA, tA)):
                R, t = trim_fit(P, Q, Rs, ts, TRIM)
                keep = np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1)) <= TRIM
                if best is None or keep.sum() > best[2].sum():
                    best = (R, t, keep)
            R, t, keep = best
            method = 'sequence'
            if not (keep.sum() >= MIN_PAIRS and keep.mean() >= MIN_FRAC):
                # fallback: structural matching of CA in the core region, started from both fits
                tgt_ca = np.array([c for _, c in rl]); core_idx = np.where(core_mask)[0]
                cand = [icp(Rs, ts, core_idx, ref_ca, tgt_ca) for Rs, ts in ((R0, t0), (RA, tA))]
                Ri, ti, pri = max(cand, key=lambda c: len(c[2]))
                if len(pri) >= MIN_PAIRS and len(pri) >= MIN_FRAC * len(core_idx):
                    core = pri; pairs = pri
                    P = np.array([rl[j][1] for i, j in core]); Q = np.array([ref_ca[i] for i, j in core])
                    R, t = kabsch(P, Q)
                    R, t = trim_fit(P, Q, R, t, TRIM)
                    keep = np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1)) <= TRIM
                    rms_all0 = rmsd(P @ R0.T + t0, Q)
                    method = 'structural CA matching (sequence alignment failed)'
            Pt = P @ R.T + t
            rms_all, rms_keep = rmsd(Pt, Q), rmsd(Pt[keep], Q[keep])
            ok = keep.sum() >= MIN_PAIRS and keep.mean() >= MIN_FRAC
            sup_rows.append([fam, pdb, chname, refpdb, refchain, ('accepted' if ok else 'rejected: %d/%d pairs within 2 A' % (keep.sum(), len(core))) + '; ' + method,
                             '%.3f' % pid, len(core), int(keep.sum()), '%.3f' % rms_all0, '%.3f' % rms_all, '%.3f' % rms_keep,
                             ' '.join(str(x) for x in np.round(np.concatenate([R.ravel(), t]), 5))])
            if not ok:
                continue
            transforms[chname] = (R, t)
            # catalytic residues: through the alignment; if that gives no residue of the reference type
            # within 3.0 A (functional atoms, after superposition), the nearest residue of the reference
            # type within 3.0 A in space is taken instead (method column says which)
            m = dict(pairs)
            FA = ('OD1', 'OD2', 'OE1', 'OE2', 'NZ', 'NE2', 'ND1', 'OH', 'OG1', 'NH1', 'NH2')

            def fdist(r, A):
                sel = [a for a in r if a.name in FA] or [a for a in r if a.name == 'CA']
                X = np.array([a.pos.tolist() for a in sel]) @ R.T + t
                return float(np.sqrt(((X[:, None] - A[None]) ** 2).sum(-1)).min())
            for (role, rc, num, ic, name, rres), (lab, A) in zip(cat, func):
                i = next(k for k, (r, _) in enumerate(rl_ref) if r is rres or (r.seqid.num == num and r.seqid.icode.strip() == ic))
                j = m.get(i)
                ra = rl[j][0] if j is not None else None
                da = fdist(ra, A) if ra is not None else None
                how = 'alignment'
                r, dd = ra, da
                if ra is None or ra.name != name or da > 3.0:
                    cands = [(fdist(x, A), x) for x, _ in rl if x.name == name]
                    cands = [c for c in cands if c[0] <= 3.0]
                    if cands:
                        dd, r = min(cands, key=lambda c: c[0]); how = 'spatial (nearest %s within 3.0 A)' % name
                if r is None:
                    catmap_rows.append([fam, pdb, chname, role, '%s%d%s' % (name, num, ic), '', '', '', '', 'gap in alignment', '', ''])
                    continue
                un = sm.get((chname, r.seqid.num, r.seqid.icode.strip()), ('', ''))
                catmap_rows.append([fam, pdb, chname, role, '%s%d%s' % (name, num, ic), r.name, r.seqid.num,
                                    r.seqid.icode.strip(), un[1], 'same type' if r.name == name else 'type differs', '%.2f' % dd, how])
        if not transforms:
            continue
        entries_ok.add(pdb)
        # nearest accepted chain for every point
        pa, pc = [], []
        for chname in transforms:
            for r, _ in chs[chname]:
                for a in r:
                    pa.append(a.pos.tolist()); pc.append((chname, r))
        tree = cKDTree(np.array(pa))

        def to_frame(x):
            d, i = tree.query(x)
            chname, r = pc[i]
            R, t = transforms[chname]
            return R @ x + t, chname, r, d

        for ch in st[0]:
            for r in ch:
                if r.is_water():
                    a = r[0]
                    x, chname, near, d = to_frame(np.array(a.pos.tolist()))
                    wat.append((pdb, ch.name, r.seqid.num, r.seqid.icode.strip(), 'HOH', x, chname, near, d, sm))
        for iid, ch, seq, ins, comp, cls, x0, y0, z0, rscc in con.execute(
                "select instance_id, auth_asym, auth_seq, ins_code, comp_id, ligand_class, x, y, z, rscc from instance "
                "where pdb=? and ligand_class in ('ion','sugar','other')", (pdb,)):
            x, chname, near, d = to_frame(np.array([x0, y0, z0]))
            rec = (pdb, ch, seq, (ins or '').strip(), comp, x, chname, near, d, sm)
            if cls == 'ion':
                ion.append(rec)
            else:
                rel, nref, nun, dref = where(x)
                un = sm.get((chname, near.seqid.num, near.seqid.icode.strip()), ('', ''))
                lig_rows.append([fam, pdb, iid, ch, seq, (ins or '').strip(), comp, cls, rscc if rscc is not None else '',
                                 chname, '%s%d%s' % (near.name, near.seqid.num, near.seqid.icode.strip()), un[1], '%.2f' % d,
                                 '%.3f' % x[0], '%.3f' % x[1], '%.3f' % x[2], rel, nref, nun, dref])
    n_ok = len(entries_ok)

    def cluster(recs, label):
        if not recs:
            return [], []
        X = np.array([r[5] for r in recs])
        tree = cKDTree(X)
        # Greedy density clustering (single linkage at 1.0 A chained the GH13 waters into one
        # 139 214-member, 45 A cluster). Seed = the point whose 1.0 A neighbourhood holds the most
        # distinct entries; cluster = unassigned points within 1.0 A of the seed; repeat while the
        # best seed still reaches the presence threshold.
        import heapq
        ent = np.array([hash(r[0]) for r in recs])
        nb = tree.query_ball_point(X, CLUST_R)
        alive = np.ones(len(X), bool)
        need = PRESENCE * n_ok

        def cnt(k):
            return len({ent[j] for j in nb[k] if alive[j]})
        heap = [(-len({ent[j] for j in nb[k]}), k) for k in range(len(X))]
        heapq.heapify(heap)
        groups = {}
        while heap:
            negc, k = heapq.heappop(heap)
            if not alive[k]:
                continue
            c = cnt(k)
            if c < -negc:
                heapq.heappush(heap, (-c, k)); continue
            if c < need:
                break
            idx = [j for j in nb[k] if alive[j]]
            alive[idx] = False
            groups[len(groups)] = idx
        out, members = [], []
        for l, idx in groups.items():
            ents = {recs[k][0] for k in idx}
            if len(ents) < PRESENCE * n_ok:
                continue
            C = X[idx].mean(0)
            spread = float(np.sqrt(((X[idx] - C) ** 2).sum(1).mean()))
            rel, nref, nun, dref = where(C)
            comps = collections.Counter(recs[k][4] for k in idx)
            out.append([fam, label, None, len(ents), n_ok, '%.3f' % (len(ents) / n_ok), len(idx), '%.2f' % spread,
                        '%.2f' % np.sqrt(((X[idx] - C) ** 2).sum(1)).max(),
                        '%.3f' % C[0], '%.3f' % C[1], '%.3f' % C[2], rel, nref, nun, dref,
                        ';'.join('%s:%d' % kv for kv in comps.most_common())])
            members.append(idx)
        order = sorted(range(len(out)), key=lambda k: -out[k][3])
        out = [out[k] for k in order]; members = [members[k] for k in order]
        mem_rows = []
        for n, (row, idx) in enumerate(zip(out, members), 1):
            row[2] = '%s_%s%d' % (fam, label[0].upper(), n)
            for k in idx:
                pdb, ch, seq, ins, comp, x, chname, near, d, sm = recs[k]
                un = sm.get((chname, near.seqid.num, near.seqid.icode.strip()), ('', ''))
                mem_rows.append([fam, row[2], pdb, ch, seq, ins, comp, chname,
                                 '%s%d%s' % (near.name, near.seqid.num, near.seqid.icode.strip()), un[1], '%.2f' % d,
                                 '%.3f' % x[0], '%.3f' % x[1], '%.3f' % x[2]])
        return out, mem_rows

    wcl, wmem = cluster(wat, 'water')
    icl, imem = cluster(ion, 'ion')
    # single-linkage chaining check: largest cluster of all
    info = {'family': fam, 'reference': refpdb, 'reference_chain': refchain, 'reference_rule': refnote,
            'catalytic_or_site_residues': ['%s %s%d%s' % (c[0], c[4], c[2], c[3]) for c in cat],
            'n_core_residues_reference': int(core_mask.sum()), 'n_entries': len(pdbs), 'n_entries_superposed': n_ok,
            'n_waters_transformed': len(wat), 'n_ions_transformed': len(ion), 'n_ligand_centres': len(lig_rows),
            'n_water_clusters_ge50pct': len(wcl), 'n_ion_clusters_ge50pct': len(icl)}
    return info, sup_rows, catmap_rows, wcl + icl, wmem + imem, lig_rows


def write(path, header, rows, mode):
    new = mode == 'w' or not os.path.exists(path)
    with open(path, mode) as f:
        if new:
            f.write('\t'.join(header) + '\n')
        for r in rows:
            f.write('\t'.join('' if v is None else str(v) for v in r) + '\n')


H_SUP = ['family', 'pdb', 'chain', 'ref_pdb', 'ref_chain', 'status', 'seq_identity_aligned', 'n_core_pairs',
         'n_pairs_retained', 'rmsd_core_first_fit_A', 'rmsd_core_all_A', 'rmsd_retained_A', 'R11..R33_t1..t3']
H_CAT = ['family', 'pdb', 'chain', 'role', 'reference_residue', 'resname', 'auth_seq', 'ins_code', 'uniprot_resnum',
         'type_check', 'functional_atom_distance_to_reference_A', 'method']
H_CL = ['family', 'kind', 'cluster_id', 'n_entries', 'n_entries_superposed', 'fraction_entries', 'n_members',
        'rms_radius_A', 'max_radius_A', 'x', 'y', 'z', 'distance_to_catalytic_A', 'nearest_reference_residue',
        'nearest_reference_uniprot_resnum', 'nearest_reference_distance_A', 'components']
H_MEM = ['family', 'cluster_id', 'pdb', 'chain', 'auth_seq', 'ins_code', 'comp_id', 'assigned_protein_chain',
         'nearest_residue', 'nearest_residue_uniprot_resnum', 'nearest_distance_A', 'x_frame', 'y_frame', 'z_frame']
H_LIG = ['family', 'pdb', 'instance_id', 'chain', 'auth_seq', 'ins_code', 'comp_id', 'ligand_class', 'rscc',
         'assigned_protein_chain', 'nearest_residue', 'nearest_residue_uniprot_resnum', 'nearest_distance_A',
         'x_frame', 'y_frame', 'z_frame', 'distance_to_catalytic_A', 'nearest_reference_residue',
         'nearest_reference_uniprot_resnum', 'nearest_reference_distance_A']

if __name__ == '__main__':
    con = db()
    fams = sys.argv[1:] or [f for f, n in con.execute("select family, count(*) from entry group by family") if n >= 2]
    infos = []
    mode = 'w'
    for fam in fams:
        info, sup, catm, cl, mem, lig = run_family(fam, con)
        print(json.dumps(info), flush=True)
        infos.append(info)
        write(os.path.join(OUT, 'superposition.tsv'), H_SUP, sup, mode)
        write(os.path.join(OUT, 'catalytic_residues_mapped.tsv'), H_CAT, catm, mode)
        write(os.path.join(OUT, 'conserved_clusters.tsv'), H_CL, cl, mode)
        write(os.path.join(OUT, 'conserved_cluster_members.tsv'), H_MEM, mem, mode)
        write(os.path.join(OUT, 'ligand_centres_in_frame.tsv'), H_LIG, lig, mode)
        mode = 'a'
    json.dump(infos, open(os.path.join(OUT, 'families.json'), 'w'), indent=1)
