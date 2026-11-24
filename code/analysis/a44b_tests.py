"""Rule test, step 2: analyses A-H on the tables of a44a_build.py.

The rule under test (short draft, Table 1): at subsite positions of GH13 and GH31, positions whose side chain
hydrogen-bonds the sugar (H) keep their residue type in every protein; aromatic positions that touch the sugar
without a side-chain hydrogen bond (Ar) stay aromatic; aliphatic positions (Al) vary.

Definitions used throughout
  entry set of subsite s     E_s: entries with a density-supported unit (RSCC >= 0.80, in_tables) at s
  contact                    a ligand atom of the unit within 4.0 A of the residue (DB contact table, ASU,
                             residue in the frame chain), residue mapped to a reference position (1AMY / 2QMJ)
  observed class per entry   H if a side-chain N/O lies within the cutoff (3.5 A unless stated) of a ligand N/O
                             of any unit at the position's qualifying subsites; otherwise Ar if the residue is
                             F, W, Y or H; otherwise Al
  residue identity           residue at the reference position in the entry's frame chain, read from the
                             structural map (Calpha within 2.0 A), whether or not it contacts the sugar;
                             engineered mutations (_struct_ref_seq_dif) dropped
  cluster vote               modal residue of the cluster's entries (ties share the vote)
  modal fraction m_q         share of cluster votes carried by the most common residue type
Run: ~/miniconda3/bin/python code/analysis/a44b_tests.py
"""
import math, random
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[2]
D = P / 'data/analysis/44_rule_test'
FAMS = ['GH13', 'GH31']
SUBS = ['-3', '-2', '-1', '+1', '+2', '+3']
INNER, OUTER = ['-1', '+1'], ['-3', '-2', '+2', '+3']
CAT = {'GH13': {'D179', 'E204', 'D289'}, 'GH31': {'D443', 'D542'}}
AROM = {'PHE', 'TRP', 'TYR', 'HIS'}
AROM3 = {'PHE', 'TRP', 'TYR'}
ONE = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', ILE='I', LEU='L', LYS='K',
           MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V')
TYPE_RULE = {**{a: 'H' for a in 'DENQKRSTH'}, **{a: 'Ar' for a in 'FWY'}, **{a: 'Al' for a in 'AVLIMPGC'}}
DRAFT = {('GH13', '-1', q) for q in ('D179', 'D289', 'H288', 'R177', 'H92', 'Y51')} | \
    {('GH13', '+1', q) for q in ('E204', 'G183', 'F180')} | {('GH13', '+2', 'W206')} | \
    {('GH31', '-1', q) for q in ('D327', 'H600', 'D443', 'D542', 'W441', 'W539', 'F575', 'I328', 'I364')}
NPERM, NBOOT, NSPLIT = 10000, 2000, 1000
rng = np.random.default_rng(20261124)
random.seed(20261124)

units = pd.read_csv(D / 'units_used.tsv', sep='\t', dtype=str, keep_default_na=False)
geo = pd.read_csv(D / 'contacts_geom.tsv', sep='\t', dtype=str, keep_default_na=False)
georef = pd.read_csv(D / 'contacts_geom_ref.tsv', sep='\t', dtype=str, keep_default_na=False)
pm = pd.read_csv(D / 'posmap.tsv', sep='\t', dtype=str, keep_default_na=False)
mut = pd.read_csv(D / 'mutations.tsv', sep='\t', dtype=str, keep_default_na=False)
cb = pd.read_csv(D / 'clusters_bound.tsv', sep='\t', dtype=str, keep_default_na=False)
ch = pd.read_csv(D / 'clusters_heldout.tsv', sep='\t', dtype=str, keep_default_na=False)

# ------------------------------------------------------------------ entry-level dictionaries
ENT = {f: {} for f in FAMS}            # pdb -> dict(cl50, cl30, chain, subs)
for r in cb.itertuples():
    ENT[r.family][r.pdb] = dict(cl50=r.cl50, cl30=r.cl30, chain=r.frame_chain, key=r.protein_key, subs=set())
for r in units.itertuples():
    ENT[r.family][r.pdb]['subs'].add(r.subsite)
MUT = {(r.pdb, r.chain, r.seq) for r in mut.itertuples()}
RES, MUTATED, REFRES = {}, set(), {f: {} for f in FAMS}
chain_of = {(f, p): e['chain'] for f in FAMS for p, e in ENT[f].items()}
for r in pm.itertuples():
    if chain_of.get((r.family, r.pdb)) == r.chain:
        RES[(r.family, r.pdb, r.position)] = r.resname
        if (r.pdb, r.chain, r.seq) in MUT:
            MUTATED.add((r.family, r.pdb, r.position))
for f, (rp, rc) in {'GH13': ('1AMY', 'A'), 'GH31': ('2QMJ', 'A')}.items():
    for r in pm[(pm.pdb == rp) & (pm.chain == rc)].itertuples():
        REFRES[f][r.position] = r.resname
# held-out entries
HELD = {f: {} for f in FAMS}
for r in ch[ch.heldout == 'True'].itertuples():
    HELD[r.family][r.pdb] = dict(cl50=r.cl50, cl30=r.cl30)
held_chain = dict(zip(ch.pdb, ch.chain))
for r in pm.itertuples():
    if r.pdb in HELD[r.family] and held_chain[r.pdb] == r.chain:
        RES[(r.family, r.pdb, r.position)] = r.resname
        if (r.pdb, r.chain, r.seq) in MUT:
            MUTATED.add((r.family, r.pdb, r.position))


def fnum(x):
    try:
        return float(x)
    except ValueError:
        return math.inf


# contact rows: (fam, pdb, sub, q) -> min side-chain polar distance (plain, angle)
CON = {}
for r in geo.itertuples():
    k = (r.family, r.pdb, r.subsite, r.position)
    a, b = fnum(r.d_sc), fnum(r.d_sc_ang)
    if k in CON:
        a, b = min(a, CON[k][0]), min(b, CON[k][1])
    CON[k] = (a, b, r.resname)
CONREF = {}
for r in georef.itertuples():
    k = (r.family, r.position)
    a = fnum(r.d_sc)
    CONREF[k] = (min(a, CONREF.get(k, (math.inf,))[0]), r.resname)


def kind(fam, pdb, subs, q, cut=3.5, ang=False):
    """observed class of q in entry pdb pooled over subsites subs; None if not contacted there or mutated"""
    if (fam, pdb, q) in MUTATED:
        return None
    best, rn = None, None
    for s in subs:
        c = CON.get((fam, pdb, s, q))
        if c is None:
            continue
        rn = c[2]
        d = c[1] if ang else c[0]
        if d <= cut:
            return 'H'
        best = 'x'
    if best is None:
        return None
    return 'Ar' if rn in AROM else 'Al'


# ------------------------------------------------------------------ A: positions
def contact_table(fam, entries):
    rows = []
    for s in SUBS:
        Es = [p for p in entries if s in ENT[fam][p]['subs']]
        if not Es:
            continue
        qs = {q for (f, p, ss, q) in CON if f == fam and ss == s and p in entries}
        for q in qs:
            hit = [p for p in Es if (fam, p, s, q) in CON]
            cl = defaultdict(list)
            for p in Es:
                cl[ENT[fam][p]['cl50']].append((fam, p, s, q) in CON)
            rows.append(dict(family=fam, subsite=s, position=q, n_entries_subsite=len(Es), n_contact=len(hit),
                             frac_entry=len(hit) / len(Es), n_cl50_subsite=len(cl),
                             frac_cl50=np.mean([np.mean(v) for v in cl.values()])))
    return rows


A = pd.DataFrame([r for f in FAMS for r in contact_table(f, list(ENT[f]))])
A['catalytic'] = [q in CAT[f] for f, q in zip(A.family, A.position)]
A = A.sort_values(['family', 'subsite', 'frac_entry'], ascending=[True, True, False])
A.to_csv(D / 'A_positions_all.tsv', sep='\t', index=False, float_format='%.3f')


def position_set(fam, T, subs_allowed=SUBS, table=None, col='frac_entry'):
    t = A if table is None else table
    t = t[(t.family == fam) & (t[col] >= T - 1e-9) & t.subsite.isin(subs_allowed)]
    out = defaultdict(list)
    for r in t.itertuples():
        out[r.position].append(r.subsite)
    return dict(out)


# ------------------------------------------------------------------ votes
def mode_vote(values):
    c = Counter(values)
    if not c:
        return None
    top = max(c.values())
    w = [k for k, v in c.items() if v == top]
    return {k: 1 / len(w) for k in w}


def votes(fam, q, subs, clcol='cl50', entries=None, contacted_only=False, cut=3.5, ang=False, held=False):
    """per cluster: residue vote dict, class vote dict"""
    src = HELD[fam] if held else ENT[fam]
    pool = entries if entries is not None else list(src)
    byc_res, byc_cls = defaultdict(list), defaultdict(list)
    tmp_res, tmp_cls = defaultdict(list), defaultdict(list)
    for p in pool:
        e = src[p]
        if not held and not (set(subs) & e['subs']):
            continue
        k = None if held else kind(fam, p, subs, q, cut, ang)
        if k is not None:
            tmp_cls[(e[clcol], e['cl50'])].append(k)
        if (fam, p, q) in MUTATED:
            continue
        if contacted_only and k is None:
            continue
        rn = RES.get((fam, p, q), 'gap')
        tmp_res[(e[clcol], e['cl50'])].append(rn)
    for tmp, dst in ((tmp_res, byc_res), (tmp_cls, byc_cls)):
        for (c, c50), items in tmp.items():
            if clcol == 'cl50':
                dst[c].extend(items)
            else:
                dst[c].append((c50, items))
    if clcol == 'cl50':
        V = {c: mode_vote(v) for c, v in byc_res.items()}
        K = {c: mode_vote(v) for c, v in byc_cls.items()}
        return V, K
    # 30 % clusters: each 50 % cluster votes first (its modal residue), then the 30 % cluster takes the
    # modal of its 50 % votes, so that a 30 % cluster is not decided by one protein deposited many times
    V, K = {}, {}
    for src_d, out in ((byc_res, V), (byc_cls, K)):
        for c, lst in src_d.items():
            acc = Counter()
            for c50, items in lst:
                for k, x in mode_vote(items).items():
                    acc[k] += x
            top = max(acc.values())
            w = [k for k, v in acc.items() if abs(v - top) < 1e-9]
            out[c] = {k: 1 / len(w) for k in w}
    return V, K


def metrics(V, w=None, drop_gap=True):
    tot, acc, gaps = 0.0, Counter(), 0.0
    for c, v in V.items():
        wc = 1.0 if w is None else w.get(c, 0)
        if wc == 0:
            continue
        if drop_gap and 'gap' in v:
            gaps += wc * v['gap']
            v = {k: x for k, x in v.items() if k != 'gap'}
            if not v:
                continue
            s = sum(v.values())
            v = {k: x / s for k, x in v.items()}
        tot += wc
        for k, x in v.items():
            acc[k] += wc * x
    if tot == 0:
        return None
    p = np.array([x / tot for x in acc.values()])
    modal = max(acc, key=acc.get)
    return dict(n=tot, modal=modal, m=acc[modal] / tot, ntypes=len(acc), H=float(-(p * np.log2(p)).sum()),
                arom=sum(acc[k] for k in AROM) / tot, arom3=sum(acc[k] for k in AROM3) / tot, gaps=gaps,
                types=', '.join(f'{ONE.get(k, k)} {acc[k]:.3g}' for k in sorted(acc, key=lambda z: -acc[z])))


def majority(K, w=None):
    acc = Counter()
    for c, v in K.items():
        if v is None:
            continue
        wc = 1.0 if w is None else w.get(c, 0)
        for k, x in v.items():
            acc[k] += wc * x
    if not acc:
        return None
    top = max(acc.values())
    best = [k for k in ('H', 'Ar', 'Al') if acc.get(k, 0) == top]
    return best[0] if len(best) == 1 else 'tie:' + '/'.join(best)


def entry_majority(fam, q, subs, cut=3.5, ang=False):
    c = Counter(k for p in ENT[fam] if set(subs) & ENT[fam][p]['subs'] for k in [kind(fam, p, subs, q, cut, ang)] if k)
    return c.most_common(1)[0][0] if c else None, dict(c)


def ref_class(fam, q):
    ref = REFRES[fam].get(q)
    t = TYPE_RULE.get(ONE.get(ref, 'X'))
    c = CONREF.get((fam, q))
    if c is None:
        return None, t
    k = 'H' if c[0] <= 3.5 else ('Ar' if c[1] in AROM else 'Al')
    return k, t


# ------------------------------------------------------------------ statistics
def stat_vals(m, a, lab):
    lab = np.asarray(lab)
    out = {}
    for name, (g1, g2, use_arom) in {'H-Al': ('H', 'Al', False), 'H-Ar': ('H', 'Ar', False),
                                     'Ar(arom)-Al': ('Ar', 'Al', True)}.items():
        x1 = (a if use_arom else m)[lab == g1]
        x2 = m[lab == g2]
        out[name] = x1.mean() - x2.mean() if len(x1) and len(x2) else np.nan
    return out


def perm_test(m, a, lab, fam, stratified=False):
    obs = stat_vals(m, a, lab)
    lab = np.asarray(lab)
    fam = np.asarray(fam)
    ge = Counter()
    ge2 = Counter()
    for _ in range(NPERM):
        if stratified:
            L = lab.copy()
            for f in set(fam):
                idx = np.where(fam == f)[0]
                L[idx] = rng.permutation(lab[idx])
        else:
            L = rng.permutation(lab)
        s = stat_vals(m, a, L)
        for k in obs:
            if not np.isnan(s[k]) and s[k] >= obs[k] - 1e-12:
                ge[k] += 1
            if not np.isnan(s[k]) and abs(s[k]) >= abs(obs[k]) - 1e-12:
                ge2[k] += 1
    return obs, {k: (ge[k] + 1) / (NPERM + 1) for k in obs}, {k: (ge2[k] + 1) / (NPERM + 1) for k in obs}


def boot_ci(posrows, clcol='cl50'):
    """posrows: list of dict(fam, cls, V). Resample clusters within family."""
    allc = {f: sorted({c for r in posrows if r['fam'] == f for c in r['V']}) for f in FAMS}
    res = defaultdict(list)
    for _ in range(NBOOT):
        w = {}
        for f in FAMS:
            if allc[f]:
                draw = Counter(rng.choice(allc[f], size=len(allc[f]), replace=True))
                w.update(draw)
        m, a, lab = [], [], []
        for r in posrows:
            x = metrics(r['V'], w)
            if x is None:
                continue
            m.append(x['m']); a.append(x['arom']); lab.append(r['cls'])
        s = stat_vals(np.array(m), np.array(a), lab)
        for k, v in s.items():
            res[k].append(v)
    return {k: (np.nanpercentile(v, 2.5), np.nanpercentile(v, 97.5)) for k, v in res.items()}


def run_test(posrows, label, rows_out, do_boot=True):
    pr = [r for r in posrows if r['cls'] in ('H', 'Ar', 'Al') and r['met'] is not None]
    m = np.array([r['met']['m'] for r in pr]); a = np.array([r['met']['arom'] for r in pr])
    lab = [r['cls'] for r in pr]; fam = [r['fam'] for r in pr]
    n = Counter(lab)
    obs, p1, p2 = perm_test(m, a, lab, fam)
    _, p1s, p2s = perm_test(m, a, lab, fam, stratified=True)
    ci = boot_ci(pr) if do_boot else {k: (np.nan, np.nan) for k in obs}
    means = {g: (m[np.array(lab) == g].mean() if n[g] else np.nan) for g in ('H', 'Ar', 'Al')}
    arom_ar = a[np.array(lab) == 'Ar'].mean() if n['Ar'] else np.nan
    for k in obs:
        rows_out.append(dict(setting=label, n_H=n['H'], n_Ar=n['Ar'], n_Al=n['Al'], mean_m_H=means['H'],
                             mean_m_Ar=means['Ar'], mean_m_Al=means['Al'], mean_arom_Ar=arom_ar, statistic=k,
                             observed=obs[k], ci_lo=ci[k][0], ci_hi=ci[k][1], p_one_sided=p1[k], p_two_sided=p2[k],
                             p_one_sided_stratified=p1s[k], p_two_sided_stratified=p2s[k]))


def build_posrows(T, clcol='cl50', cut=3.5, ang=False, subs_allowed=SUBS, class_mode='cluster', held=False,
                  contacted_only=False, kind_subs=None):
    rows = []
    for f in FAMS:
        for q, subs in sorted(position_set(f, T, subs_allowed).items()):
            ks = kind_subs or subs
            V, K = votes(f, q, ks, clcol, cut=cut, ang=ang, contacted_only=contacted_only)
            cls_cl = majority(K)
            cls_en, cnt = entry_majority(f, q, ks, cut, ang)
            rk, rt = ref_class(f, q)
            cls = {'cluster': cls_cl, 'entry': cls_en, 'ref_entry': rk if rk else rt, 'ref_type': rt}[class_mode]
            if held:
                V, _ = votes(f, q, ks, clcol, held=True)
            rows.append(dict(fam=f, q=q, subs=subs, cls=cls, cls_cluster=cls_cl, cls_entry=cls_en, entry_counts=cnt,
                             cls_ref_entry=rk, cls_ref_type=rt, ref_res=REFRES[f].get(q, ''), V=V, met=metrics(V),
                             cat=q in CAT[f]))
    return rows


def ident_table(rows, rows30=None):
    out = []
    for i, r in enumerate(rows):
        x = r['met'] or {}
        d = dict(family=r['fam'], position=r['q'], subsites=','.join(r['subs']), catalytic=r['cat'],
                 ref_residue=r['ref_res'], class_cluster50=r['cls_cluster'], class_entry=r['cls_entry'],
                 entry_kind_counts=' '.join(f'{k}:{v}' for k, v in sorted(r['entry_counts'].items())),
                 class_ref_entry=r['cls_ref_entry'], class_ref_type=r['cls_ref_type'],
                 n_cl50=x.get('n'), modal50=ONE.get(x.get('modal'), x.get('modal')), m50=x.get('m'),
                 ntypes50=x.get('ntypes'), entropy50=x.get('H'), arom50=x.get('arom'), arom3_50=x.get('arom3'),
                 gap_votes50=x.get('gaps'), types50=x.get('types'))
        if rows30:
            y = rows30[i]['met'] or {}
            d.update(n_cl30=y.get('n'), modal30=ONE.get(y.get('modal'), y.get('modal')), m30=y.get('m'),
                     ntypes30=y.get('ntypes'), entropy30=y.get('H'), arom30=y.get('arom'), types30=y.get('types'))
        out.append(d)
    return pd.DataFrame(out)


TESTS = []
# ------------------------------------------------------------------ A + C + D (main)
for T in (0.2, 0.5, 0.8):
    r50 = build_posrows(T, 'cl50')
    r30 = build_posrows(T, 'cl30')
    ident_table(r50, r30).to_csv(D / f'C_identity_T{int(T*100)}.tsv', sep='\t', index=False, float_format='%.3f')
    for cm in ('cluster', 'entry', 'ref_entry', 'ref_type'):
        for cl, rr in (('cl50', r50), ('cl30', r30)):
            rows = [dict(r, cls={'cluster': r['cls_cluster'], 'entry': r['cls_entry'],
                                 'ref_entry': r['cls_ref_entry'] or r['cls_ref_type'], 'ref_type': r['cls_ref_type']}[cm])
                    for r in rr]
            boot = (cm in ('cluster', 'ref_entry') and cl == 'cl50') or T == 0.2
            run_test(rows, f'T{int(T*100)} {cl} class={cm} all', TESTS, do_boot=boot)
            run_test([r for r in rows if not r['cat']], f'T{int(T*100)} {cl} class={cm} no-catalytic', TESTS, do_boot=boot)
    print('done T', T, flush=True)

# contacted-only identity (draft-like) at T80/T20
for T in (0.2, 0.8):
    rows = build_posrows(T, 'cl50', contacted_only=True)
    run_test(rows, f'T{int(T*100)} cl50 class=cluster all identity-contacted-only', TESTS, do_boot=False)
    run_test([r for r in rows if not r['cat']], f'T{int(T*100)} cl50 class=cluster no-catalytic identity-contacted-only',
             TESTS, do_boot=False)

# ------------------------------------------------------------------ B(iii) held-out proteins
for T in (0.2, 0.5, 0.8):
    for cl in ('cl50', 'cl30'):
        rows = build_posrows(T, cl, held=True)
        if T == 0.2 and cl == 'cl50':
            ident_table(rows).to_csv(D / 'B3_heldout_identity_T20.tsv', sep='\t', index=False, float_format='%.3f')
        run_test(rows, f'HELDOUT T{int(T*100)} {cl} class=cluster(bound) all', TESTS, do_boot=(cl == 'cl50'))
        run_test([r for r in rows if not r['cat']], f'HELDOUT T{int(T*100)} {cl} class=cluster(bound) no-catalytic',
                 TESTS, do_boot=(cl == 'cl50'))
print('held-out done', flush=True)

# ------------------------------------------------------------------ E: control at the outer subsites
for T in (0.2, 0.5):
    rows = build_posrows(T, 'cl50', subs_allowed=OUTER)
    inner = {('GH13', q) for q in position_set('GH13', 0.2, INNER)} | {('GH31', q) for q in position_set('GH31', 0.2, INNER)}
    ident_table(rows).assign(also_inner=[(r['fam'], r['q']) in inner for r in rows]).to_csv(
        D / f'E_outer_identity_T{int(T*100)}.tsv', sep='\t', index=False, float_format='%.3f')
    run_test(rows, f'OUTER T{int(T*100)} cl50 class=cluster(outer contacts) all', TESTS)
    strict = [r for r in rows if (r['fam'], r['q']) not in inner]
    run_test(strict, f'OUTER T{int(T*100)} cl50 class=cluster(outer contacts) only-outer-positions', TESTS)
    # apolar pooled: H vs (Ar+Al)
    pooled = [dict(r, cls=('H' if r['cls'] == 'H' else ('Al' if r['cls'] in ('Ar', 'Al') else r['cls']))) for r in strict]
    run_test(pooled, f'OUTER T{int(T*100)} cl50 only-outer-positions H vs apolar(Ar+Al as Al)', TESTS)
# inner H vs outer H, same rules (T20)
inn = [r for r in build_posrows(0.2, 'cl50', subs_allowed=INNER) if r['cls'] == 'H']
out_ = [r for r in build_posrows(0.2, 'cl50', subs_allowed=OUTER) if r['cls'] == 'H']
innq = {(r['fam'], r['q']) for r in inn}
out_ = [r for r in out_ if (r['fam'], r['q']) not in innq]
EH = []
for tag, sel in (('all', lambda r: True), ('no-catalytic', lambda r: not r['cat'])):
    I = [r for r in inn if sel(r) and r['met']]; O = [r for r in out_ if r['met']]
    mi = np.array([r['met']['m'] for r in I]); mo = np.array([r['met']['m'] for r in O])
    obs = mi.mean() - mo.mean()
    allm = np.concatenate([mi, mo]); ge = 0
    for _ in range(NPERM):
        z = rng.permutation(allm)
        ge += (z[:len(mi)].mean() - z[len(mi):].mean()) >= obs - 1e-12
    EH.append(dict(comparison=f'H inner(-1,+1) vs H outer-only(-3,-2,+2,+3), T20, {tag}', n_inner=len(I), n_outer=len(O),
                   mean_m_inner=mi.mean(), mean_m_outer=mo.mean(), diff=obs, p_one_sided=(ge + 1) / (NPERM + 1),
                   inner_positions=' '.join(f"{r['fam']}:{r['q']}" for r in I),
                   outer_positions=' '.join(f"{r['fam']}:{r['q']}" for r in O)))
pd.DataFrame(EH).to_csv(D / 'E_inner_vs_outer_H.tsv', sep='\t', index=False, float_format='%.4f')
print('E done', flush=True)

# ------------------------------------------------------------------ F: cutoff / angle sensitivity
F_rows = []
base = {(r['fam'], r['q']): r['cls'] for r in build_posrows(0.2, 'cl50')}
for T in (0.2, 0.5, 0.8):
    base = {(r['fam'], r['q']): r['cls'] for r in build_posrows(T, 'cl50')}
    for cut, ang in ((3.2, False), (3.5, False), (4.0, False), (3.5, True)):
        rows = build_posrows(T, 'cl50', cut=cut, ang=ang)
        ch_ = [f"{r['fam']}:{r['q']} {base[(r['fam'], r['q'])]}->{r['cls']}" for r in rows if r['cls'] != base[(r['fam'], r['q'])]]
        F_rows.append(dict(T=T, cutoff=cut, angle=ang, n_positions=len(rows), n_changed_vs_3p5=len(ch_), changed='; '.join(ch_)))
        tag = f'SENS T{int(T*100)} cut{cut}{" angle" if ang else ""} cl50 class=cluster'
        run_test(rows, tag + ' all', TESTS, do_boot=(T == 0.2))
        run_test([r for r in rows if not r['cat']], tag + ' no-catalytic', TESTS, do_boot=(T == 0.2))
pd.DataFrame(F_rows).to_csv(D / 'F_class_changes.tsv', sep='\t', index=False)
print('F done', flush=True)

# ------------------------------------------------------------------ B(ii): split-half over 50 % clusters
cl_of = {f: defaultdict(list) for f in FAMS}
for f in FAMS:
    for p, e in ENT[f].items():
        cl_of[f][e['cl50']].append(p)
SPL, SPLPOS = [], []
for T in (0.2, 0.5, 0.8):
    posclass = defaultdict(Counter)
    for it in range(NSPLIT):
        posrows = []
        for f in FAMS:
            cls_ = sorted(cl_of[f]); perm = rng.permutation(cls_)
            hA, hB = set(perm[:len(perm) // 2]), set(perm[len(perm) // 2:])
            eA = [p for c in hA for p in cl_of[f][c]]; eB = [p for c in hB for p in cl_of[f][c]]
            tA = pd.DataFrame(contact_table(f, eA))
            if tA.empty:
                continue
            for q, subs in position_set(f, T, table=tA, col='frac_cl50').items():
                _, K = votes(f, q, subs, 'cl50', entries=eA)
                cls = majority(K)
                V, _ = votes(f, q, subs, 'cl50', entries=eB)
                posclass[(f, q)][cls] += 1
                posrows.append(dict(fam=f, q=q, cls=cls, V=V, met=metrics(V), cat=q in CAT[f]))
        for tag, pr in (('all', posrows), ('no-catalytic', [r for r in posrows if not r['cat']])):
            pr = [r for r in pr if r['cls'] in ('H', 'Ar', 'Al') and r['met']]
            m = np.array([r['met']['m'] for r in pr]); a = np.array([r['met']['arom'] for r in pr]); lab = [r['cls'] for r in pr]
            s = stat_vals(m, a, lab); n = Counter(lab)
            SPL.append(dict(T=T, subset=tag, it=it, n_H=n['H'], n_Ar=n['Ar'], n_Al=n['Al'], **s,
                            mean_m_H=m[np.array(lab) == 'H'].mean() if n['H'] else np.nan,
                            mean_m_Al=m[np.array(lab) == 'Al'].mean() if n['Al'] else np.nan,
                            all_H_invariant=bool(n['H']) and bool((m[np.array(lab) == 'H'] == 1).all())))
    for (f, q), c in posclass.items():
        SPLPOS.append(dict(T=T, family=f, position=q, n_splits_selected=sum(c.values()),
                           **{f'class_{k}': c.get(k, 0) for k in ('H', 'Ar', 'Al')},
                           other=sum(v for k, v in c.items() if k not in ('H', 'Ar', 'Al'))))
    print('split T', T, flush=True)
spl = pd.DataFrame(SPL)
spl.to_csv(D / 'B2_split_half_iterations.tsv', sep='\t', index=False, float_format='%.4f')
pd.DataFrame(SPLPOS).sort_values(['T', 'family', 'position']).to_csv(D / 'B2_split_half_position_classes.tsv', sep='\t', index=False)
summ = []
for (T, sub), g in spl.groupby(['T', 'subset']):
    for k in ('H-Al', 'H-Ar', 'Ar(arom)-Al'):
        v = g[k].dropna()
        summ.append(dict(T=T, subset=sub, statistic=k, n_iter=len(v), median=v.median(), p2_5=v.quantile(0.025),
                         p97_5=v.quantile(0.975), frac_positive=(v > 0).mean(), frac_le_zero=(v <= 0).mean(),
                         median_n_H=g.n_H.median(), median_n_Ar=g.n_Ar.median(), median_n_Al=g.n_Al.median(),
                         frac_iter_all_H_invariant=g.all_H_invariant.mean()))
pd.DataFrame(summ).to_csv(D / 'B2_split_half_summary.tsv', sep='\t', index=False, float_format='%.4f')

pd.DataFrame(TESTS).to_csv(D / 'D_tests.tsv', sep='\t', index=False, float_format='%.4f')

# ------------------------------------------------------------------ G: GH31 per protein
fun = {'A2RM80': 'alpha-xylosidase', 'A5ZY13': 'GH31 (B. obeum, 3PHA)', 'L0N7E5': 'alpha-glucosidase (sugar beet)',
       'O43451': 'alpha-glucosidase (human MGAM-C)', 'P10253': 'alpha-glucosidase (human GAA)',
       'Q9STC1': 'alpha-1,4-glucan lyase'}
G = []
for q, subs in sorted(position_set('GH31', 0.2).items()):
    row = dict(position=q, subsites=','.join(subs), catalytic=q in CAT['GH31'])
    for k, lab in fun.items():
        ents = [p for p, e in ENT['GH31'].items() if e['key'] == k and set(subs) & e['subs']]
        rs = Counter(ONE.get(RES.get(('GH31', p, q), 'gap'), '-') + ('*' if ('GH31', p, q) in MUTATED else '') for p in ents)
        row[lab] = ' '.join(f'{a}{"" if n == 1 else n}' for a, n in rs.most_common()) if ents else 'n/a'
        ks = Counter(kind('GH31', p, subs, q) for p in ents)
        row[lab + ' kind'] = ' '.join(f'{a}:{n}' for a, n in ks.items() if a)
    glc = [fun[k] for k in ('L0N7E5', 'O43451', 'P10253')]
    vals = {row[g].replace('*', '').split()[0][0] if row[g] not in ('n/a', '') else None for g in glc}
    vals.discard(None)
    row['varies_within_alpha_glucosidases'] = len(vals) > 1
    G.append(row)
pd.DataFrame(G).to_csv(D / 'G_gh31_per_protein.tsv', sep='\t', index=False)

# ------------------------------------------------------------------ H: corrected 80 % table (per entry, as the draft)
Hrows = []
for f in FAMS:
    a = A[(A.family == f) & (A.frac_entry >= 0.8 - 1e-9)]
    extra = A[(A.family == f) & (A.position == 'G183') & (A.subsite == '+1')]
    for r in pd.concat([a, extra]).itertuples():
        Es = [p for p in ENT[f] if r.subsite in ENT[f][p]['subs']]
        ks = Counter(); rt = Counter(); rt_all = Counter(); mutd = []
        for p in Es:
            k = kind(f, p, [r.subsite], r.position)
            rn = RES.get((f, p, r.position), 'gap')
            if (f, p, r.position) in MUTATED:
                mutd.append(p)
                continue
            rt_all[ONE.get(rn, rn)] += 1
            if k:
                ks[k] += 1; rt[ONE.get(rn, rn)] += 1
        Hrows.append(dict(family=f, subsite=r.subsite, position=r.position, catalytic=r.catalytic,
                          n_entries_subsite=r.n_entries_subsite, n_contact=r.n_contact, frac=r.frac_entry,
                          meets_80=r.frac_entry >= 0.8 - 1e-9, in_draft_table=(f, r.subsite, r.position) in DRAFT,
                          kinds_contacted=' '.join(f'{k}:{v}' for k, v in sorted(ks.items())),
                          majority_kind=ks.most_common(1)[0][0] if ks else '',
                          residues_contacted=' '.join(f'{k}{v}' for k, v in rt.most_common()),
                          residues_all_entries_at_subsite=' '.join(f'{k}{v}' for k, v in rt_all.most_common()),
                          mutant_entries_excluded=' '.join(mutd)))
pd.DataFrame(Hrows).to_csv(D / 'H_table80_corrected.tsv', sep='\t', index=False, float_format='%.3f')

# mutations touching any position used
mu = []
for (f, p, q) in sorted(MUTATED):
    if p in ENT[f]:
        mu.append(dict(family=f, pdb=p, position=q, protein=ENT[f][p]['key']))
pd.DataFrame(mu).to_csv(D / 'mutations_at_mapped_positions_bound.tsv', sep='\t', index=False)
print('all done')
