"""Revision analyses for the IJBM editorial reviews (2026-11-24 adjudication, items 1, 3, 4, 5, 10, 11).

All inputs are existing result tables; nothing is recomputed from structures.
  A. Inner/outer test of hydrogen-bonding (H) positions with ONE stated rule: a position is "at -1/+1" when the
     subsites listed for it in Table 1 (column `subsites` of C_identity_T20.tsv, the subsites at which it is
     contacted in >= 20% of the entries) include -1 or +1; otherwise it is "outer-only". With this rule GH13 Lys182
     (+1,+2) is inner. Permutation of the inner/outer labels over the H positions (10 000), plain and within family;
     one- and two-sided p. Same with the catalytic residues removed.
  A2. H-Al difference with the class labels permuted within location strata (inner / outer-only) and within
     family x location, and with positions seen in fewer than 5 sequence groups removed.
  A3. Summary of the H-Al variants already in D_tests.tsv (two-sided and family-stratified p), including the
     unfavourable ones (reference class at 50/80% thresholds, held-out without catalytic residues, 30% clusters).
  B. The two modal fractions of one position: m (Table 1; one vote per 50% group, groups of the entries with a unit at
     -3..+3) and m_map (map section; counted over all ligand-bound entries of the family), side by side.
  C. The map against experiment: score of the published map, of a trivial map (every oxygen at -1/+1 not
     tolerated, every other oxygen tolerated) and of maps with the states permuted within family (1000 permutations),
     on the literature cells of Table S6; agreement by state and by substituent.
  D. The 50% sequence groups (clusters) of GH13 and GH31: UniProt, EC, protein, organism, entries; GH13 subfamily
     from the local CAZy dump (cazy_data.zip) through UniProt -> NCBI protein accessions (UniProt REST, cached).
  E. Numbering concordance of the GH13 positions of Table 1: barley AMY2 reference (1AMY) <-> TAKA-amylase (7TAA,
     6TAA, 2TAA) <-> porcine pancreatic alpha-amylase (1PPI), from the structural position map (posmap.tsv), and the
     conserved sequence regions I-IV located by their motifs on the TAKA sequence (SIFTS table of the database).
Outputs: data/analysis/44_rule_test/revision/
Run: ~/miniconda3/bin/python code/analysis/a44d_revision.py   (about 1 min; D needs network once, then cached)
"""
import json, sqlite3, urllib.request, zipfile
from pathlib import Path
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[2]
R = P / 'data'
D = R / 'analysis/44_rule_test'
O = D / 'revision'
O.mkdir(exist_ok=True)
rng = np.random.default_rng(20261124)
NPERM = 10000
CAT = {'GH13': {'D179', 'E204', 'D289'}, 'GH31': {'D443', 'D542'}}
out_lines = []


def say(s=''):
    print(s)
    out_lines.append(s)


C = pd.read_csv(D / 'C_identity_T20.tsv', sep='\t', dtype={'subsites': str})
C['inner'] = C.subsites.apply(lambda s: bool({'-1', '+1'} & set(s.split(','))))
C['cat'] = [q in CAT[f] for f, q in zip(C.family, C.position)]


def perm_diff(x, lab, strata=None, n=NPERM):
    """mean(x[lab]) - mean(x[~lab]); labels permuted (within strata if given). Returns obs, p1, p2."""
    x = np.asarray(x, float); lab = np.asarray(lab, bool)
    obs = x[lab].mean() - x[~lab].mean()
    ge = ge2 = 0
    for _ in range(n):
        if strata is None:
            L = rng.permutation(lab)
        else:
            L = lab.copy()
            for s in set(strata):
                idx = np.where(np.asarray(strata) == s)[0]
                L[idx] = rng.permutation(lab[idx])
        if L.all() or (~L).all():
            continue
        d = x[L].mean() - x[~L].mean()
        ge += d >= obs - 1e-12
        ge2 += abs(d) >= abs(obs) - 1e-12
    return obs, (ge + 1) / (n + 1), (ge2 + 1) / (n + 1)


# ------------------------------------------------------------------ A: inner vs outer H, one rule
say('== A. H positions at -1/+1 (Table 1 subsites include -1 or +1) vs outer-only')
H = C[C.class_cluster50 == 'H'].copy()
rows = []
for tag, sub in (('all', H), ('no-catalytic', H[~H.cat])):
    obs, p1, p2 = perm_diff(sub.m50, sub.inner)
    _, p1s, p2s = perm_diff(sub.m50, sub.inner, strata=list(sub.family))
    I, Oo = sub[sub.inner], sub[~sub.inner]
    rows.append(dict(set=tag, n_inner=len(I), n_outer=len(Oo), mean_m_inner=I.m50.mean(), mean_m_outer=Oo.m50.mean(),
                     diff=obs, p_one_sided=p1, p_two_sided=p2, p_one_sided_family_strat=p1s,
                     p_two_sided_family_strat=p2s,
                     inner=' '.join(f'{f}:{q}({m:.2f})' for f, q, m in zip(I.family, I.position, I.m50)),
                     outer=' '.join(f'{f}:{q}({m:.2f})' for f, q, m in zip(Oo.family, Oo.position, Oo.m50))))
A = pd.DataFrame(rows)
A.to_csv(O / 'A_inner_vs_outer_H.tsv', sep='\t', index=False, float_format='%.4f')
for r in rows:
    say(f"  {r['set']}: inner {r['n_inner']} mean {r['mean_m_inner']:.3f}; outer {r['n_outer']} mean "
        f"{r['mean_m_outer']:.3f}; diff {r['diff']:+.3f}; p1 {r['p_one_sided']:.4f} p2 {r['p_two_sided']:.4f}; "
        f"family-strat p1 {r['p_one_sided_family_strat']:.4f} p2 {r['p_two_sided_family_strat']:.4f}")
inn = H[H.inner]
say(f"  H positions at -1/+1: {len(inn)}; m = 1: {(inn.m50 >= 0.9999).sum()}; "
    f"variable (m < 0.9): {', '.join(f'{f} {q} {m:.2f}' for f, q, m in zip(inn.family, inn.position, inn.m50) if m < 0.9)}")

# ------------------------------------------------------------------ A2: H-Al with location strata
say('\n== A2. H - Al, class labels permuted within strata')
HA = C[C.class_cluster50.isin(['H', 'Al'])].copy()
HA['loc'] = np.where(HA.inner, 'inner', 'outer')
rows = []
for tag, sub in (('all', HA), ('no-catalytic', HA[~HA.cat]), ('groups>=5', HA[HA.n_cl50 >= 5]),
                 ('groups>=5 no-catalytic', HA[(HA.n_cl50 >= 5) & ~HA.cat])):
    lab = (sub.class_cluster50 == 'H').values
    o0, p1, p2 = perm_diff(sub.m50, lab)
    _, p1f, p2f = perm_diff(sub.m50, lab, strata=list(sub.family))
    _, p1l, p2l = perm_diff(sub.m50, lab, strata=list(sub['loc']))
    _, p1fl, p2fl = perm_diff(sub.m50, lab, strata=[f + l for f, l in zip(sub.family, sub['loc'])])
    nH, nAl = lab.sum(), (~lab).sum()
    rows.append(dict(set=tag, n_H=nH, n_Al=nAl, diff=o0, p1=p1, p2=p2, p1_family=p1f, p2_family=p2f,
                     p1_location=p1l, p2_location=p2l, p1_family_location=p1fl, p2_family_location=p2fl))
    say(f"  {tag}: H {nH} Al {nAl} diff {o0:+.3f}; p2 {p2:.4f}; family p2 {p2f:.4f}; location p1 {p1l:.4f} "
        f"p2 {p2l:.4f}; family x location p2 {p2fl:.4f}")
pd.DataFrame(rows).to_csv(O / 'A2_H_Al_location_strata.tsv', sep='\t', index=False, float_format='%.4f')
# inner-only H vs Al (Al at -1/+1)
sub = HA[HA.inner]
lab = (sub.class_cluster50 == 'H').values
o0, p1, p2 = perm_diff(sub.m50, lab)
say(f"  at -1/+1 only: H {lab.sum()} Al {(~lab).sum()} diff {o0:+.3f} p1 {p1:.4f} p2 {p2:.4f}")
pd.DataFrame([dict(set='inner only', n_H=lab.sum(), n_Al=(~lab).sum(), diff=o0, p1=p1, p2=p2)]).to_csv(
    O / 'A2b_H_Al_inner_only.tsv', sep='\t', index=False, float_format='%.4f')

# ------------------------------------------------------------------ A3: variants from D_tests.tsv
say('\n== A3. H-Al variants (D_tests.tsv)')
T = pd.read_csv(D / 'D_tests.tsv', sep='\t')
T = T[T.statistic == 'H-Al']
pick = ['T20 cl50 class=cluster all', 'T20 cl50 class=cluster no-catalytic', 'T20 cl30 class=cluster all',
        'T50 cl50 class=cluster all', 'T80 cl50 class=cluster all',
        'T20 cl50 class=ref_entry all', 'T50 cl50 class=ref_entry all', 'T80 cl50 class=ref_entry all',
        'T50 cl50 class=ref_entry no-catalytic', 'T80 cl50 class=ref_entry no-catalytic',
        'HELDOUT T20 cl50 class=cluster(bound) all', 'HELDOUT T20 cl50 class=cluster(bound) no-catalytic',
        'HELDOUT T20 cl30 class=cluster(bound) all', 'OUTER T20 cl50 class=cluster(outer contacts) all',
        'SENS T20 cut3.2 cl50 class=cluster all', 'SENS T20 cut4.0 cl50 class=cluster all',
        'SENS T20 cut3.5 angle cl50 class=cluster all']
V = T[T.setting.isin(pick)].set_index('setting').loc[pick].reset_index()
V = V[['setting', 'n_H', 'n_Ar', 'n_Al', 'observed', 'ci_lo', 'ci_hi', 'p_one_sided', 'p_two_sided',
       'p_one_sided_stratified', 'p_two_sided_stratified']]
V.to_csv(O / 'A3_H_Al_variants.tsv', sep='\t', index=False, float_format='%.4f')
for r in V.itertuples():
    say(f"  {r.setting}: {r.observed:+.3f} p1 {r.p_one_sided:.4f} p2 {r.p_two_sided:.4f} strat p2 "
        f"{r.p_two_sided_stratified:.4f}")

# ------------------------------------------------------------------ B: the two modal fractions
say('\n== B. m (Table 1) and m_map (map section)')
PI = pd.read_csv(R / 'analysis/45_hydroxyl_map/position_identity.tsv', sep='\t')
B = PI[PI.m50_a44_T20.notna() & PI.family.isin(['GH13', 'GH31'])][
    ['family', 'position', 'n_cl50_a44', 'm50_a44_T20', 'n_cl50', 'm50', 'types']].rename(
    columns={'n_cl50_a44': 'groups_m', 'm50_a44_T20': 'm', 'n_cl50': 'groups_m_map', 'm50': 'm_map',
             'types': 'votes_m_map'})
B['abs_diff'] = (B.m - B.m_map).abs()
B.to_csv(O / 'B_m_vs_m_map.tsv', sep='\t', index=False, float_format='%.3f')
for r in B[B.abs_diff > 0.005].itertuples():
    say(f"  {r.family} {r.position}: m {r.m:.2f} ({r.groups_m:.0f} groups) / m_map {r.m_map:.2f} "
        f"({r.groups_m_map:.0f} groups)")

# ------------------------------------------------------------------ C: map against experiment, baseline
say('\n== C. Map against experiment: observed vs trivial map vs permuted states')
# Outcome per cell and substituent, transcribed from Table S6 (ijbm_SI.tex) / substituted_starch_hydrolysis_literature_v2.md.
# T: tolerated by every enzyme tested; N: tolerated by none; M: tolerated by some enzymes of the family, not by others.
# Low-level acceptance by one lineage only (Okada 2000, <= 3% of cuts / very slow reaction) is coded 'n' (not
# tolerated, single lineage: no enzyme contrast).
LIT = {
    ('GH13', '-3', 'O3'): {'P': 'T', 'glycosyl': 'T'},
    ('GH13', '-3', 'O6'): {'P': 'T', 'glycosyl': 'T', 'other': 'T'},          # 6-PA, 6-iodo, benzylidene
    ('GH13', '-2', 'O3'): {'P': 'M', 'glycosyl': 'n'},
    ('GH13', '-2', 'O6'): {'P': 'M', 'glycosyl': 'n'},
    ('GH13', '-1', 'O2'): {'HE': 'N', 'Me': 'N'},
    ('GH13', '-1', 'O3'): {'HE': 'N', 'P': 'N', 'Me': 'N'},
    ('GH13', '-1', 'O6'): {'HE': 'N', 'P': 'N', 'Me': 'N'},
    ('GH13', '+1', 'O3'): {'P': 'N'},
    ('GH13', '+2', 'O3'): {'P': 'M'},
    ('GH13', '+2', 'O6'): {'P': 'T', 'other': 'T'},                           # 6-amino-6-deoxy
    ('GH14', '-2', 'O3'): {'P': 'N'},
    ('GH14', '-2', 'O6'): {'P': 'N'},
    ('GH14', '-1', 'O6'): {'P': 'N'},
    ('GH14', '+1', 'O3'): {'P': 'N'},
    ('GH14', '+1', 'O6'): {'P': 'T'},
    ('GH15', '-1', 'O2'): {'Me': 'N'},
    ('GH15', '-1', 'O3'): {'Me': 'N', 'P': 'N'},
    ('GH15', '-1', 'O6'): {'Me': 'N', 'P': 'N', 'other': 'N'},                # 6'-F, 6'-Br
    ('GH15', '+1', 'O2'): {'Me': 'T'},
    ('GH15', '+1', 'O3'): {'P': 'N', 'Me': 'N'},
    ('GH15', '+2', 'O3'): {'P': 'T'},
    # cells without a prediction in the published map
    ('GH13', '+1', 'O6'): {'P': 'M'},
    ('GH13', '+3', 'O3'): {'P': 'T'},
    ('GH14', '-1', 'O3'): {'P': 'N'},
    ('GH15', '+1', 'O6'): {'P': 'T', 'other': 'M'},
}


def combine(d):
    v = set(d.values())
    if v <= {'N', 'n'}:
        return 'N'
    if v == {'T'}:
        return 'T'
    return 'M'


OUT = {k: combine(v) for k, v in LIT.items()}


def agree(state, outcome):
    """Scoring rule of the paper: INV agrees when no enzyme tolerates; EXP when at least one tolerates; VAR when enzymes
    differ or at least one tolerates. Other states make no prediction (None)."""
    if state == 'INV':
        return outcome == 'N'
    if state in ('EXP', 'VAR'):
        return outcome in ('T', 'M')
    return None


CS = pd.read_csv(R / 'analysis/45_hydroxyl_map/cell_summary.tsv', sep='\t', dtype={'subsite': str})
CS = CS[(CS.set == 'all_glucoselike') & CS.label.isin(['O2', 'O3', 'O6'])]
MAP = {(r.family, r.subsite, r.label): r.state for r in CS.itertuples()}
MAP[('GH15', '+1', 'O6')] = 'none'   # no O6 deposited at +1 in GH15


def score(m, cells):
    res = [agree(m.get(c, 'none'), OUT[c]) for c in cells]
    pred = [x for x in res if x is not None]
    return sum(pred), len(pred)


cells21 = [c for c in LIT if agree(MAP.get(c, 'none'), OUT[c]) is not None]
cells_all = list(LIT)
a_obs, n_obs = score(MAP, cells21)
say(f'  published map: {a_obs}/{n_obs} cells agree')
by_state = {}
for c in cells21:
    s = MAP[c]
    by_state.setdefault(s, [0, 0])
    by_state[s][0] += agree(s, OUT[c]); by_state[s][1] += 1
say('  by state: ' + '; '.join(f'{s} {a}/{n}' for s, (a, n) in by_state.items()))
TRIV = {c: ('INV' if c[1] in ('-1', '+1') else 'EXP') for c in MAP}
TRIV[('GH15', '+1', 'O6')] = 'INV'
a_t21, n_t21 = score(TRIV, cells21)
a_tall, n_tall = score(TRIV, cells_all)
say(f'  trivial map (-1/+1 blocked, rest tolerated): {a_t21}/{n_t21} on the 21 cells; {a_tall}/{n_tall} on all 25')
# INV-type cells of the trivial map
tri_inv = [c for c in cells21 if TRIV[c] == 'INV']
say(f'  trivial map, "blocked" cells among the 21: {sum(OUT[c] == "N" for c in tri_inv)}/{len(tri_inv)}')
# permutation of states within family
fam_cells = {}
for c, s in MAP.items():
    if s != 'none':
        fam_cells.setdefault(c[0], []).append(c)
perm_frac, perm_agree, perm_n = [], [], []
for _ in range(1000):
    M = dict(MAP)
    for f, cs in fam_cells.items():
        st = rng.permutation([MAP[c] for c in cs])
        for c, s in zip(cs, st):
            M[c] = s
    a, n = score(M, cells_all)
    perm_agree.append(a); perm_n.append(n); perm_frac.append(a / n if n else np.nan)
perm_frac = np.array(perm_frac); perm_agree = np.array(perm_agree)
obs_frac = a_obs / n_obs
p_perm = (np.sum(perm_frac >= obs_frac - 1e-12) + 1) / (len(perm_frac) + 1)
say(f'  permuted states (1000, within family, scored on the 25 literature cells): agreement fraction mean '
    f'{np.nanmean(perm_frac):.3f}, 95% range {np.nanpercentile(perm_frac, 2.5):.3f}-{np.nanpercentile(perm_frac, 97.5):.3f}; '
    f'cells with a prediction mean {np.mean(perm_n):.1f}; P(fraction >= {obs_frac:.3f}) = {p_perm:.3f}')
# falsifiable (INV) part under permutation: fraction of INV-assigned literature cells that are N
inv_frac = []
for _ in range(1000):
    M = dict(MAP)
    for f, cs in fam_cells.items():
        st = rng.permutation([MAP[c] for c in cs])
        for c, s in zip(cs, st):
            M[c] = s
    ii = [c for c in cells_all if M.get(c) == 'INV']
    if ii:
        inv_frac.append(np.mean([OUT[c] == 'N' for c in ii]))
inv_frac = np.array(inv_frac)
obs_inv = by_state['INV'][0] / by_state['INV'][1]
p_inv = (np.sum(inv_frac >= obs_inv - 1e-12) + 1) / (len(inv_frac) + 1)
say(f'  INV cells only: observed {by_state["INV"][0]}/{by_state["INV"][1]} = {obs_inv:.3f}; permuted mean '
    f'{inv_frac.mean():.3f}, 95% range {np.percentile(inv_frac, 2.5):.3f}-{np.percentile(inv_frac, 97.5):.3f}; '
    f'P(>= observed) = {p_inv:.3f}')
# base rate: share of the 25 literature cells that are not tolerated
say(f'  base rate: {sum(v == "N" for v in OUT.values())} of {len(OUT)} literature cells not tolerated')
# per substituent
rows = []
for sub in ('P', 'Me', 'HE', 'glycosyl', 'other'):
    cs = [c for c in cells21 if sub in LIT[c]]
    res = []
    for c in cs:
        o = LIT[c][sub]
        o2 = 'N' if o in ('N', 'n') else o
        s = MAP[c]
        if s == 'VAR' and o == 'n':
            res.append(None)          # one lineage only, no enzyme contrast: uninformative for a VAR cell
        else:
            res.append(agree(s, o2))
    a = sum(1 for x in res if x); n = sum(1 for x in res if x is not None)
    only = [c for c in cs if set(LIT[c]) == {sub}]
    rows.append(dict(substituent=sub, cells=len(cs), informative=n, agree=a, cells_resting_on_it_alone=len(only)))
    say(f'  {sub}: {len(cs)} cells, {a}/{n} informative agree, {len(only)} rest on it alone')
# without phosphate data
res = []
for c in cells21:
    d = {k: v for k, v in LIT[c].items() if k != 'P'}
    if not d:
        continue
    s = MAP[c]
    if s == 'VAR' and set(d.values()) <= {'n'}:
        res.append(None); continue
    res.append(agree(s, combine(d)))
a = sum(1 for x in res if x); n = sum(1 for x in res if x is not None)
say(f'  without phosphate data: {len(res)} cells keep a result, {a}/{n} informative agree')
rows.append(dict(substituent='all but phosphate', cells=len(res), informative=n, agree=a, cells_resting_on_it_alone=''))
pd.DataFrame(rows).to_csv(O / 'C_agreement_by_substituent.tsv', sep='\t', index=False)
pd.DataFrame([dict(map='published', cells=n_obs, agree=a_obs),
              dict(map='trivial (-1/+1 blocked), 21 cells', cells=n_t21, agree=a_t21),
              dict(map='trivial (-1/+1 blocked), 25 cells', cells=n_tall, agree=a_tall),
              dict(map='permuted within family, mean fraction', cells=np.mean(perm_n), agree=np.nanmean(perm_frac)),
              dict(map='permuted, P(fraction >= observed)', cells='', agree=p_perm),
              dict(map='INV cells observed fraction', cells=by_state['INV'][1], agree=obs_inv),
              dict(map='INV cells permuted mean fraction', cells='', agree=inv_frac.mean()),
              dict(map='INV cells permuted P(>= observed)', cells='', agree=p_inv)]).to_csv(
    O / 'C_baseline.tsv', sep='\t', index=False, float_format='%.3f')
pd.DataFrame([dict(family=c[0], subsite=c[1], oxygen=c[2], state=MAP.get(c, 'none'), trivial=TRIV.get(c),
                   outcome=OUT[c], evidence=' '.join(f'{k}:{v}' for k, v in LIT[c].items()),
                   agree_published=agree(MAP.get(c, 'none'), OUT[c]), agree_trivial=agree(TRIV.get(c), OUT[c]))
              for c in cells_all]).to_csv(O / 'C_cells.tsv', sep='\t', index=False)

# ------------------------------------------------------------------ D: clusters
say('\n== D. 50% sequence groups of GH13 and GH31')
cb = pd.read_csv(D / 'clusters_bound.tsv', sep='\t', dtype=str)
con = sqlite3.connect(f'file:{R}/db/enzimas_experimental.sqlite?mode=ro', uri=True)
ent = pd.read_sql('select pdb, EC, protein, organism, uniprot from entry', con)
cb = cb.merge(ent, on='pdb', how='left')
cache_f = O / 'uniprot_ncbi_cache.json'
cache = json.loads(cache_f.read_text()) if cache_f.exists() else {}


def ncbi_of(u):
    if u in cache:
        return cache[u]
    try:
        with urllib.request.urlopen(f'https://rest.uniprot.org/uniprotkb/{u}.json?fields=xref_embl,xref_refseq',
                                    timeout=30) as r:
            d = json.load(r)
        ids = []
        for x in d.get('uniProtKBCrossReferences', []):
            if x['database'] == 'EMBL':
                ids += [p['value'] for p in x.get('properties', []) if p['key'] == 'ProteinId' and p['value'] != '-']
            elif x['database'] == 'RefSeq':
                ids.append(x['id'])
        cache[u] = ids
    except Exception as e:     # network failure: leave the subfamily empty
        say(f'  UniProt lookup failed for {u}: {e}')
        return None
    return cache[u]


keys = sorted(set(cb.protein_key))
acc2key = {}
for k in keys:
    ids = ncbi_of(k)
    for a in ids or []:
        acc2key.setdefault(a, set()).add(k)
cache_f.write_text(json.dumps(cache, indent=0))
sub_of = {}
zf = R.parent / 'cazy_data.zip'
if zf.exists() and acc2key:
    with zipfile.ZipFile(zf) as z, z.open(z.namelist()[0]) as f:
        for line in f:
            p = line.decode('utf-8', 'replace').rstrip('\n').split('\t')
            if len(p) >= 4 and p[3] in acc2key and (p[0].startswith('GH13') or p[0].startswith('GH31')):
                for k in acc2key[p[3]]:
                    sub_of.setdefault(k, set()).add(p[0])
rows = []
for (fam, cl), g in cb.groupby(['family', 'cl50']):
    prots = sorted(set(g.protein_key))
    first = g.iloc[0]
    subs = sorted(set().union(*[sub_of.get(k, set()) for k in prots]) - {fam}) or sorted(
        set().union(*[sub_of.get(k, set()) for k in prots]))
    names = '; '.join(sorted(set(str(x).split(' (')[0] for x in g.protein)))
    orgs = '; '.join(sorted(set(str(x) for x in g.organism)))
    ecs = '; '.join(sorted(set(str(x) for x in g.EC)))
    rows.append(dict(family=fam, group=cl, group_30=';'.join(sorted(set(g.cl30))), uniprot=' '.join(prots),
                     EC=ecs, protein=names, organism=orgs, cazy_subfamily=' '.join(subs) if subs else '',
                     entries=len(g), pdb=' '.join(sorted(g.pdb))))
CL = pd.DataFrame(rows)
CL['n'] = CL.group.str.extract(r'_(\d+)$').astype(int)
CL = CL.sort_values(['family', 'n']).drop(columns='n')
CL.to_csv(O / 'D_clusters.tsv', sep='\t', index=False)
say(f'  {len(CL)} groups ({(CL.family == "GH13").sum()} GH13, {(CL.family == "GH31").sum()} GH31); with a CAZy '
    f'subfamily: {(CL.cazy_subfamily != "").sum()}; largest group: {CL.entries.max()} entries')

# ------------------------------------------------------------------ E: numbering concordance
say('\n== E. Numbering concordance (GH13)')
pm = pd.read_csv(D / 'posmap.tsv', sep='\t', dtype=str)
sifts = pd.read_sql("select pdb, chain, auth_seq, resname from sifts_residue where pdb in ('7TAA','1PPI')", con)
ONE = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', ILE='I', LEU='L', LYS='K',
           MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V')
taka = sifts[(sifts.pdb == '7TAA') & (sifts.chain == 'A')].copy()
taka['n'] = pd.to_numeric(taka.auth_seq, errors='coerce')
taka = taka.dropna(subset=['n']).drop_duplicates('n').set_index('n')
# CSR I-IV as TAKA-amylase ranges (MacGregor 2001, Janecek 2014); checked against the motif read from 7TAA
CSR = {'I': (117, 122, 'DVVANH'), 'II': (202, 210, 'GLRIDTVKH'), 'III': (230, 233, 'EVLD'), 'IV': (292, 297, 'FVENHD')}
for k, (a, b, motif) in CSR.items():
    seq = ''.join(ONE.get(taka.resname.get(i, ''), 'x') for i in range(a, b + 1))
    say(f'  CSR {k} TAKA {a}-{b}: 7TAA reads {seq}; expected {motif}; {"OK" if seq == motif else "MISMATCH"}')
CSR_OK = {k: v for k, v in CSR.items()
          if ''.join(ONE.get(taka.resname.get(i, ''), 'x') for i in range(v[0], v[1] + 1)) == v[2]}
gh13 = C[C.family == 'GH13']
rows = []
for r in gh13.itertuples():
    q = r.position
    def look(pdb):
        x = pm[(pm.pdb == pdb) & (pm.chain == 'A') & (pm.position == q)]
        return (x.resname.iloc[0].capitalize() + x.seq.iloc[0], float(x.d_ca.iloc[0])) if len(x) else ('--', np.nan)
    t7, d7 = look('7TAA'); t6, _ = look('6TAA'); t2, _ = look('2TAA'); pp, dp = look('1PPI')
    tk = t7 if t7 != '--' else (t6 if t6 != '--' else t2)
    csr = ''
    if tk != '--':
        n = int(''.join(ch for ch in tk if ch.isdigit()))
        csr = ' '.join(k for k, (a, b, _) in CSR_OK.items() if a <= n <= b)
    rows.append(dict(position_1AMY=q, cls=r.class_cluster50, subsites=r.subsites, m=r.m50, TAKA=tk, d_ca_7TAA=d7,
                     PPA_1PPI=pp, d_ca_1PPI=dp, CSR=csr))
E = pd.DataFrame(rows)
E.to_csv(O / 'E_numbering_concordance.tsv', sep='\t', index=False, float_format='%.2f')
for r in E.itertuples():
    say(f'  {r.position_1AMY}: TAKA {r.TAKA}, PPA {r.PPA_1PPI}, CSR {r.CSR}')

(O / 'REPORT.txt').write_text('\n'.join(out_lines) + '\n')
print('done')
