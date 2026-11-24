"""Round-2 revision analyses for the IJBM manuscript (second round of editorial reviews, 2026-11-24).

Inputs are existing result tables (44_rule_test, 45_hydroxyl_map, 47_md_hydroxyl); nothing is recomputed from
structures. Outputs: data/analysis/44_rule_test/revision2/
  A. Every H-Al (and H-Ar, H-H) permutation test reported in the manuscript, rerun with ONE fixed seed (20261124)
     and 100 000 permutations, so that one p value is quoted. Null: the class labels (H, Ar, Al) are permuted over
     all positions that have a class in the setting (as in a44b_tests.py), plainly, within family and within
     location (-1/+1 vs outer-only; a position is at -1/+1 when its listed subsites include -1 or +1). Statistic:
     mean m of H positions minus mean m of Al (or Ar) positions. The H-H row permutes the location label over the
     H positions. Settings other than the main one are rebuilt with the functions of a44b_tests.py (definitions
     executed, its output-writing lines skipped).                                      -> A_tests_100k.tsv
  B. A simple model of m: logistic regression of (m = 1) on location (-1/+1 vs outer-only), family (GH31 vs GH13)
     and contact class (Ar, Al vs H), statsmodels GLM binomial; and an ordinal (proportional-odds logit) model of
     m in three levels (m = 1; 0.6 <= m < 1; m < 0.6), statsmodels OrderedModel.          -> B_model_m.tsv
  C. Map against experiment: exact McNemar (binomial) test published vs trivial map on the 21 cells; INV agreement
     per family; score without the cells resting only on abstracts; score with only full-text evidence;
     permuted maps rerun with 10 000 permutations (fixed seed).          -> C_mcnemar_family_abstract.tsv, C_cells_ev.tsv
  D. GH13 map restricted to the 13 sequence groups that contain an endo-alpha-amylase (EC 3.2.1.1 in
     revision/D_clusters.tsv): cell states re-aggregated from atoms.tsv with the rules of a45_hydroxyl_map.py
     (checked: the same code on all 29 groups reproduces cell_summary.tsv).               -> D_alpha_amylase_map.tsv
  E. Molecular dynamics, the three systems with the glucan at the catalytic subsites of a family with a crystal map
     (7DCG, 6A0K: GH13; 1J18: GH14): frames, ns, cells keeping the crystal state.         -> E_md_three_systems.tsv
Run: ~/miniconda3/bin/python code/analysis/a44e_round2.py   (a few minutes)
"""
from pathlib import Path
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[2]
R = P / 'data'
D = R / 'analysis/44_rule_test'
O = D / 'revision2'
O.mkdir(exist_ok=True)
SEED = 20261124
NPERM = 100_000
rng = np.random.default_rng(SEED)
lines = []


def say(s=''):
    print(s, flush=True)
    lines.append(s)


# ------------------------------------------------------------------ permutation machinery (vectorised)
def perm_labels(lab, strata, n):
    """n x len(lab) matrix of permuted labels, permuted within strata (None = one stratum)."""
    lab = np.asarray(lab)
    L = np.empty((n, len(lab)), dtype=lab.dtype)
    st = np.zeros(len(lab), int) if strata is None else pd.factorize(np.asarray(strata))[0]
    for s in np.unique(st):
        idx = np.where(st == s)[0]
        keys = rng.random((n, len(idx)))
        L[:, idx] = lab[idx][np.argsort(keys, axis=1)]
    return L


def diff_of(x, L, a, b):
    A = (L == a); B = (L == b)
    na, nb = A.sum(1), B.sum(1)
    ok = (na > 0) & (nb > 0)
    with np.errstate(invalid='ignore', divide='ignore'):
        d = (A * x).sum(1) / na - (B * x).sum(1) / nb
    return d[ok]


def ptest(x, lab, a, b, strata=None, n=NPERM):
    x = np.asarray(x, float); lab = np.asarray(lab)
    obs = x[lab == a].mean() - x[lab == b].mean()
    codes = np.where(lab == a, 0, np.where(lab == b, 1, 2)).astype(np.int8)
    d = diff_of(x, perm_labels(codes, strata, n), 0, 1)
    p1 = (np.sum(d >= obs - 1e-12) + 1) / (len(d) + 1)
    p2 = (np.sum(np.abs(d) >= abs(obs) - 1e-12) + 1) / (len(d) + 1)
    return obs, p1, p2


# ------------------------------------------------------------------ A. tests, 100 000 permutations
say(f'== A. Permutation tests, seed {SEED}, {NPERM} permutations')
src = (P / 'code/analysis/a44b_tests.py').read_text()
defs = src.split('TESTS = []')[0].replace("A.to_csv(D / 'A_positions_all.tsv', sep='\\t', index=False, float_format='%.3f')",
                                          'pass')
NS = {'__file__': str(P / 'code/analysis/a44b_tests.py')}
exec(compile(defs, 'a44b_defs', 'exec'), NS)
build_posrows, INNER_S, OUTER_S = NS['build_posrows'], NS['INNER'], NS['OUTER']


def frame(rows, class_key='cls'):
    out = []
    for r in rows:
        if r['met'] is None or r[class_key] not in ('H', 'Ar', 'Al'):
            continue
        out.append(dict(family=r['fam'], position=r['q'], cls=r[class_key], m=r['met']['m'], n=r['met']['n'],
                        cat=r['cat'], inner=bool({'-1', '+1'} & set(r['subs']))))
    return pd.DataFrame(out)


def ref_rows(rows):
    return [dict(r, cls=r['cls_ref_entry'] or r['cls_ref_type']) for r in rows]


SET = {}
r20 = build_posrows(0.2, 'cl50')
SET['All positions'] = frame(r20)
# check against Table 1
C = pd.read_csv(D / 'C_identity_T20.tsv', sep='\t', dtype={'subsites': str})
chk = SET['All positions'].merge(C[['family', 'position', 'm50']], on=['family', 'position'])
assert np.allclose(chk.m, chk.m50, atol=1e-3), 'rebuilt m differs from C_identity_T20.tsv'
F0 = SET['All positions']
SET['Without catalytic'] = F0[~F0.cat]
SET['Positions with >=5 groups'] = F0[F0.n >= 5]
SET['At -1/+1 only'] = F0[F0.inner]
SET['30% identity groups'] = frame(build_posrows(0.2, 'cl30'))
SET['Contact criterion 50%'] = frame(build_posrows(0.5, 'cl50'))
SET['Contact criterion 80%'] = frame(build_posrows(0.8, 'cl50'))
SET['Class from reference, 20%'] = frame(ref_rows(r20))
r50 = build_posrows(0.5, 'cl50'); r80 = build_posrows(0.8, 'cl50')
SET['Class from reference, 50%'] = frame(ref_rows(r50))
SET['Class from reference, 80%'] = frame(ref_rows(r80))
SET['Class from reference, 50%, without catalytic'] = SET['Class from reference, 50%'].query('~cat')
SET['Class from reference, 80%, without catalytic'] = SET['Class from reference, 80%'].query('~cat')
hz = frame(build_posrows(0.2, 'cl50', held=True))
SET['Held-out proteins'] = hz
SET['Held-out, without catalytic'] = hz[~hz.cat]
SET['Held-out, 30% identity groups'] = frame(build_posrows(0.2, 'cl30', held=True))
SET['Outer subsites (outer contacts)'] = frame(build_posrows(0.2, 'cl50', subs_allowed=OUTER_S))
SET['H-bond 3.2 A'] = frame(build_posrows(0.2, 'cl50', cut=3.2))
SET['H-bond 4.0 A'] = frame(build_posrows(0.2, 'cl50', cut=4.0))
SET['H-bond 3.5 A + angle'] = frame(build_posrows(0.2, 'cl50', ang=True))

rows = []
for name, F in SET.items():
    stats = [('H-Al', 'H', 'Al')] + ([('H-Ar', 'H', 'Ar')] if name == 'All positions' else [])
    for st, a, b in stats:
        obs, p1, p2 = ptest(F.m, F.cls, a, b)
        _, _, p2f = ptest(F.m, F.cls, a, b, strata=F.family)
        if F.inner.nunique() > 1:
            _, _, p2l = ptest(F.m, F.cls, a, b, strata=F.inner)
            _, _, p2fl = ptest(F.m, F.cls, a, b, strata=F.family + F.inner.astype(str))
        else:
            p2l = p2fl = np.nan
        rows.append(dict(setting=name, statistic=st, n_a=int((F.cls == a).sum()), n_b=int((F.cls == b).sum()),
                         n_positions=len(F), diff=obs, p_one_sided=p1, p_two_sided=p2, p_two_sided_family=p2f,
                         p_two_sided_location=p2l, p_two_sided_family_location=p2fl))
        say(f'  {name:45s} {st}: {int((F.cls == a).sum())}/{int((F.cls == b).sum())} diff {obs:+.3f} p2 {p2:.4f} '
            f'fam {p2f:.4f} loc {p2l:.4f}')
# H at -1/+1 vs outer-only (location label permuted over H positions)
Hh = F0[F0.cls == 'H']
for tag, sub in (('H at -1/+1 vs outer-only', Hh), ('H at -1/+1 vs outer-only, without catalytic', Hh[~Hh.cat])):
    lab = np.where(sub.inner, 'in', 'out')
    obs, p1, p2 = ptest(sub.m, lab, 'in', 'out')
    _, _, p2f = ptest(sub.m, lab, 'in', 'out', strata=sub.family)
    rows.append(dict(setting=tag, statistic='H-H', n_a=int(sub.inner.sum()), n_b=int((~sub.inner).sum()),
                     n_positions=len(sub), diff=obs, p_one_sided=p1, p_two_sided=p2, p_two_sided_family=p2f,
                     p_two_sided_location=np.nan, p_two_sided_family_location=np.nan))
    say(f'  {tag}: {sub.inner.sum()}/{(~sub.inner).sum()} means {sub.m[sub.inner].mean():.3f} '
        f'{sub.m[~sub.inner].mean():.3f} diff {obs:+.3f} p2 {p2:.4f} fam {p2f:.4f}')
AT = pd.DataFrame(rows)
AT.to_csv(O / 'A_tests_100k.tsv', sep='\t', index=False, float_format='%.4f')

# ------------------------------------------------------------------ B. model of m
say('\n== B. Model of m on location, family and contact class (35 positions with a class)')
import statsmodels.api as sm
from statsmodels.miscmodels.ordinal_model import OrderedModel
M = F0.copy()
M['m1'] = (M.m >= 0.9999).astype(int)
M['outer'] = (~M.inner).astype(int)
M['GH31'] = (M.family == 'GH31').astype(int)
M['Ar'] = (M.cls == 'Ar').astype(int)
M['Al'] = (M.cls == 'Al').astype(int)
X = sm.add_constant(M[['outer', 'GH31', 'Ar', 'Al']].astype(float))
brows = []
try:
    lg = sm.GLM(M.m1, X, family=sm.families.Binomial()).fit()
    say(f'  logistic (m = 1), n = {len(M)}, events {M.m1.sum()}; deviance {lg.deviance:.2f}')
    for k in X.columns:
        lo, hi = lg.conf_int().loc[k]
        brows.append(dict(model='logistic, m = 1', term=k, coef=lg.params[k], se=lg.bse[k], ci_lo=lo, ci_hi=hi,
                          odds_ratio=np.exp(lg.params[k]), p=lg.pvalues[k], n=len(M), events=int(M.m1.sum())))
        say(f'    {k:6s} coef {lg.params[k]:+.2f} se {lg.bse[k]:.2f} OR {np.exp(lg.params[k]):.2f} '
            f'[{np.exp(lo):.2f}, {np.exp(hi):.2f}] p {lg.pvalues[k]:.3f}')
except Exception as e:
    say(f'  logistic failed: {e}')
lev = np.where(M.m >= 0.9999, 2, np.where(M.m >= 0.6, 1, 0))
try:
    om = OrderedModel(pd.Series(pd.Categorical(lev, ordered=True)), M[['outer', 'GH31', 'Ar', 'Al']].astype(float),
                      distr='logit').fit(method='bfgs', disp=False)
    say(f'  ordinal (m < 0.6 / 0.6-1 / 1; counts {np.bincount(lev).tolist()})')
    ci = om.conf_int()
    for k in ['outer', 'GH31', 'Ar', 'Al']:
        brows.append(dict(model='ordinal logit, 3 levels', term=k, coef=om.params[k], se=om.bse[k], ci_lo=ci.loc[k, 0],
                          ci_hi=ci.loc[k, 1], odds_ratio=np.exp(om.params[k]), p=om.pvalues[k], n=len(M),
                          events=''))
        say(f'    {k:6s} coef {om.params[k]:+.2f} se {om.bse[k]:.2f} OR {np.exp(om.params[k]):.2f} p {om.pvalues[k]:.3f}')
except Exception as e:
    say(f'  ordinal failed: {e}')
# cross-tabulation behind the model
xt = M.groupby(['inner', 'cls']).agg(n=('m', 'size'), m1=('m1', 'sum'), mean_m=('m', 'mean')).reset_index()
say('  cross-tab (inner, class): ' + '; '.join(f'{"in" if r.inner else "out"} {r.cls}: {r.m1}/{r.n} m=1, mean {r.mean_m:.2f}'
                                         for r in xt.itertuples()))
xt.to_csv(O / 'B_crosstab.tsv', sep='\t', index=False, float_format='%.3f')
pd.DataFrame(brows).to_csv(O / 'B_model_m.tsv', sep='\t', index=False, float_format='%.4f')

# ------------------------------------------------------------------ C. map against experiment
say('\n== C. Map against experiment')
cells = pd.read_csv(D / 'revision/C_cells.tsv', sep='\t', dtype={'subsite': str})
# Evidence level per cell and substituent, from Table S6 (ijbm_SI.tex): FT if at least one full-text source gives the
# result for that substituent, ABS if only abstracts do.
EV = {
    ('GH13', '-3', 'O3'): {'P': 'FT', 'glycosyl': 'FT'},
    ('GH13', '-3', 'O6'): {'P': 'FT', 'glycosyl': 'FT', 'other': 'FT'},
    ('GH13', '-2', 'O3'): {'P': 'FT', 'glycosyl': 'FT'},
    ('GH13', '-2', 'O6'): {'P': 'FT', 'glycosyl': 'FT'},
    ('GH13', '-1', 'O2'): {'HE': 'ABS', 'Me': 'FT'},
    ('GH13', '-1', 'O3'): {'HE': 'ABS', 'P': 'FT', 'Me': 'FT'},
    ('GH13', '-1', 'O6'): {'HE': 'ABS', 'P': 'FT', 'Me': 'FT'},
    ('GH13', '+1', 'O3'): {'P': 'FT'},
    ('GH13', '+2', 'O3'): {'P': 'FT'},
    ('GH13', '+2', 'O6'): {'P': 'ABS', 'other': 'ABS'},
    ('GH14', '-2', 'O3'): {'P': 'FT'}, ('GH14', '-2', 'O6'): {'P': 'FT'}, ('GH14', '-1', 'O6'): {'P': 'FT'},
    ('GH14', '+1', 'O3'): {'P': 'FT'}, ('GH14', '+1', 'O6'): {'P': 'FT'},
    ('GH15', '-1', 'O2'): {'Me': 'ABS'},
    ('GH15', '-1', 'O3'): {'Me': 'ABS', 'P': 'FT'},
    ('GH15', '-1', 'O6'): {'Me': 'ABS', 'P': 'FT', 'other': 'FT'},
    ('GH15', '+1', 'O2'): {'Me': 'ABS'},
    ('GH15', '+1', 'O3'): {'P': 'FT', 'Me': 'ABS'},
    ('GH15', '+2', 'O3'): {'P': 'FT'},
    ('GH13', '+1', 'O6'): {'P': 'FT'}, ('GH13', '+3', 'O3'): {'P': 'FT'}, ('GH14', '-1', 'O3'): {'P': 'FT'},
    ('GH15', '+1', 'O6'): {'P': 'FT', 'other': 'FT'},
}


def parse_ev(s):
    return dict(x.split(':') for x in s.split())


def combine(vals):
    v = set(vals)
    if v <= {'N', 'n'}:
        return 'N'
    if v == {'T'}:
        return 'T'
    return 'M'


def agree(state, outcome):
    if state == 'INV':
        return outcome == 'N'
    if state in ('EXP', 'VAR'):
        return outcome in ('T', 'M')
    return None


cells['key'] = list(zip(cells.family, cells.subsite, cells.oxygen))
cells['evid'] = cells.evidence.apply(parse_ev)
cells['ev_level'] = [' '.join(f'{k}:{EV[c][k]}' for k in e) for c, e in zip(cells.key, cells.evid)]
cells['abstract_only'] = [all(EV[c][k] == 'ABS' for k in e) for c, e in zip(cells.key, cells.evid)]
ft = []
for c, e, s in zip(cells.key, cells.evid, cells.state):
    d = {k: v for k, v in e.items() if EV[c][k] == 'FT'}
    if not d:
        ft.append(None); continue
    if s == 'VAR' and set(d.values()) <= {'n'}:
        ft.append(None); continue
    ft.append(agree(s, combine(d.values())))
cells['agree_fulltext_only'] = ft
pred = cells[cells.agree_published.notna()].copy()
pred['okp'] = pred.agree_published.astype(str) == 'True'
pred['okt'] = pred.agree_trivial.astype(str) == 'True'
b = int((pred.okp & ~pred.okt).sum()); c_ = int((~pred.okp & pred.okt).sum())
from scipy.stats import binomtest
pm = binomtest(b, b + c_, 0.5).pvalue if b + c_ else 1.0
say(f'  21 cells: published {pred.okp.sum()}/{len(pred)}, trivial {pred.okt.sum()}/{len(pred)}; discordant: published '
    f'right/trivial wrong {b}, reverse {c_}; exact McNemar (binomial) two-sided p = {pm:.3f}')
inv = pred[pred.state == 'INV']
tri_blk = pred[pred.trivial == 'INV']
say(f'  INV cells {inv.okp.sum()}/{len(inv)}; trivial "blocked" cells {(tri_blk.outcome == "N").sum()}/{len(tri_blk)}')
fam_rows = []
for f, g in inv.groupby('family'):
    fam_rows.append(dict(item='INV agreement', family=f, agree=int(g.okp.sum()), cells=len(g)))
    say(f'  INV {f}: {g.okp.sum()}/{len(g)}')
na = pred[~pred.abstract_only]
nai = na[na.state == 'INV']
say(f'  without abstract-only cells ({", ".join(" ".join(k) for k in pred[pred.abstract_only].key)}): '
    f'{na.okp.sum()}/{len(na)}; INV {nai.okp.sum()}/{len(nai)}')
ftp = pred[pred.agree_fulltext_only.notna()]
fti = ftp[ftp.state == 'INV']
say(f'  full-text evidence only: {int(ftp.agree_fulltext_only.astype(bool).sum())}/{len(ftp)}; INV '
    f'{int(fti.agree_fulltext_only.astype(bool).sum())}/{len(fti)}')
# permuted maps, 10 000
CS = pd.read_csv(R / 'analysis/45_hydroxyl_map/cell_summary.tsv', sep='\t', dtype={'subsite': str})
CS = CS[(CS.set == 'all_glucoselike') & CS.label.isin(['O2', 'O3', 'O6'])]
MAP = {(r.family, r.subsite, r.label): r.state for r in CS.itertuples()}
OUTC = dict(zip(cells.key, cells.outcome))
fam_cells = {}
for k, s in MAP.items():
    fam_cells.setdefault(k[0], []).append(k)
NP2 = 10_000
fr, invf = [], []
for _ in range(NP2):
    Mp = dict(MAP)
    for f, cs in fam_cells.items():
        st = rng.permutation([MAP[k] for k in cs])
        Mp.update(dict(zip(cs, st)))
    res = [agree(Mp.get(k, 'none'), OUTC[k]) for k in OUTC]
    pr = [x for x in res if x is not None]
    fr.append(np.mean(pr) if pr else np.nan)
    ii = [OUTC[k] == 'N' for k in OUTC if Mp.get(k) == 'INV']
    if ii:
        invf.append(np.mean(ii))
fr, invf = np.array(fr), np.array(invf)
obs_f, obs_i = pred.okp.mean(), inv.okp.mean()
p_f = (np.sum(fr >= obs_f - 1e-12) + 1) / (np.isfinite(fr).sum() + 1)
p_i = (np.sum(invf >= obs_i - 1e-12) + 1) / (len(invf) + 1)
say(f'  permuted maps ({NP2}): mean fraction {np.nanmean(fr):.3f} (95% {np.nanpercentile(fr, 2.5):.3f}-'
    f'{np.nanpercentile(fr, 97.5):.3f}); P(>= {obs_f:.3f}) = {p_f:.4f}; INV mean {invf.mean():.3f} (95% '
    f'{np.percentile(invf, 2.5):.3f}-{np.percentile(invf, 97.5):.3f}); P(>= {obs_i:.3f}) = {p_i:.4f}')
summ = fam_rows + [
    dict(item='published vs trivial, 21 cells', family='all', agree=int(pred.okp.sum()), cells=len(pred)),
    dict(item='trivial, 21 cells', family='all', agree=int(pred.okt.sum()), cells=len(pred)),
    dict(item='discordant published right / trivial wrong', family='all', agree=b, cells=''),
    dict(item='discordant trivial right / published wrong', family='all', agree=c_, cells=''),
    dict(item='exact McNemar two-sided p', family='all', agree=round(pm, 4), cells=''),
    dict(item='trivial blocked cells not tolerated', family='all', agree=int((tri_blk.outcome == 'N').sum()),
         cells=len(tri_blk)),
    dict(item='without abstract-only cells', family='all', agree=int(na.okp.sum()), cells=len(na)),
    dict(item='without abstract-only cells, INV', family='all', agree=int(nai.okp.sum()), cells=len(nai)),
    dict(item='full-text evidence only', family='all', agree=int(ftp.agree_fulltext_only.astype(bool).sum()),
         cells=len(ftp)),
    dict(item='full-text evidence only, INV', family='all', agree=int(fti.agree_fulltext_only.astype(bool).sum()),
         cells=len(fti)),
    dict(item=f'permuted maps ({NP2}) mean fraction', family='all', agree=round(float(np.nanmean(fr)), 3), cells=''),
    dict(item='permuted maps 2.5-97.5 percentile', family='all',
         agree=f'{np.nanpercentile(fr, 2.5):.3f}-{np.nanpercentile(fr, 97.5):.3f}', cells=''),
    dict(item='permuted maps P(fraction >= published)', family='all', agree=round(float(p_f), 4), cells=''),
    dict(item='permuted maps INV mean fraction', family='all', agree=round(float(invf.mean()), 3), cells=''),
    dict(item='permuted maps INV 2.5-97.5 percentile', family='all',
         agree=f'{np.percentile(invf, 2.5):.3f}-{np.percentile(invf, 97.5):.3f}', cells=''),
    dict(item='permuted maps P(INV fraction >= published)', family='all', agree=round(float(p_i), 4), cells=''),
]
pd.DataFrame(summ).to_csv(O / 'C_mcnemar_family_abstract.tsv', sep='\t', index=False)
cells.drop(columns=['key', 'evid']).to_csv(O / 'C_cells_ev.tsv', sep='\t', index=False)

# ------------------------------------------------------------------ D. alpha-amylase-only GH13 map
say('\n== D. GH13 map over the endo-alpha-amylase groups')
CL = pd.read_csv(D / 'revision/D_clusters.tsv', sep='\t')
AAG = set(CL[(CL.family == 'GH13') & CL.EC.astype(str).str.contains(r'3\.2\.1\.1(?![0-9])')].group)
say(f'  alpha-amylase groups: {len(AAG)}; entries {int(CL[CL.group.isin(AAG)].entries.sum())}')
AT_ = pd.read_csv(R / 'analysis/45_hydroxyl_map/atoms.tsv', sep='\t', dtype={'subsite': str})
REL_EXPOSED, SASA_EXPOSED = 0.25, 5.0


def cell_state(gg):
    e = gg.groupby('pdb').agg(held=('held', 'mean'), sc_inv=('sc_inv', 'mean'), sc_var=('sc_var', 'mean'),
                              mc=('mc', 'mean'), rel=('rel', 'median'), sasa=('sasa_cplx', 'median'),
                              cl50=('cl50', 'first'))
    for k in ['held', 'sc_inv', 'sc_var', 'mc']:
        e[k] = e[k] >= 0.5
    e['mc_only'] = e.mc & ~e.sc_inv & ~e.sc_var
    c = e.groupby('cl50').agg(held=('held', 'mean'), sc_inv=('sc_inv', 'mean'), sc_var=('sc_var', 'mean'),
                              mc_only=('mc_only', 'mean'))
    for k in ['held', 'sc_inv', 'sc_var', 'mc_only']:
        c[k] = c[k] >= 0.5
    h, i, v, m = c.held.mean(), c.sc_inv.mean(), c.sc_var.mean(), c.mc_only.mean()
    relm = float(np.nanmedian(e.rel)) if e.rel.notna().any() else np.nan
    if h >= 0.5:
        st = 'INV' if i >= 0.5 else ('VAR' if v >= m else 'MC')
    else:
        st = 'EXP' if (relm >= REL_EXPOSED and float(e.sasa.median()) >= SASA_EXPOSED) else 'BUR'
    return st, len(e), len(c), h, float(e.sasa.median())


OXS = ['O2', 'O3', 'O6']
SUBS = ['-3', '-2', '-1', '+1', '+2', '+3']
g13 = AT_[(AT_.family == 'GH13') & AT_.label.isin(OXS)]
full = CS[CS.family == 'GH13'].set_index(['subsite', 'label'])
drows, mism = [], 0
for s in SUBS:
    for lab in OXS:
        gg = g13[(g13.subsite == s) & (g13.label == lab)]
        if gg.empty:
            continue
        st_all, *_ = cell_state(gg)
        if st_all != full.loc[(s, lab), 'state']:
            mism += 1
        ga = gg[gg.cl50.isin(AAG)]
        if ga.empty:
            drows.append(dict(subsite=s, oxygen=lab, state_29_groups=st_all, groups_29=int(full.loc[(s, lab), 'n_cl50']),
                              state_alpha_amylase='no unit', groups_alpha=0, entries_alpha=0, frac_groups_held_alpha=np.nan,
                              median_sasa_alpha=np.nan, changes=True))
            continue
        st, ne, nc, h, sa = cell_state(ga)
        drows.append(dict(subsite=s, oxygen=lab, state_29_groups=st_all, groups_29=int(full.loc[(s, lab), 'n_cl50']),
                          state_alpha_amylase=st, groups_alpha=nc, entries_alpha=ne, frac_groups_held_alpha=round(h, 2),
                          median_sasa_alpha=round(sa, 1), changes=st != st_all))
DA = pd.DataFrame(drows)
say(f'  check: re-aggregated states differing from cell_summary.tsv over all 29 groups: {mism}')
for r in DA.itertuples():
    say(f'  {r.subsite} {r.oxygen}: all {r.state_29_groups} ({r.groups_29}) -> alpha-amylase {r.state_alpha_amylase} '
        f'({r.groups_alpha} groups, held {r.frac_groups_held_alpha}){"  CHANGES" if r.changes else ""}')
DA.to_csv(O / 'D_alpha_amylase_map.tsv', sep='\t', index=False)

# ------------------------------------------------------------------ E. MD, three compared systems
say('\n== E. MD, the three systems compared with a crystal map')
SY = pd.read_csv(R / 'analysis/47_md_hydroxyl/systems.tsv', sep='\t')
MC = pd.read_csv(R / 'analysis/47_md_hydroxyl/md_cells.tsv', sep='\t', dtype={'subsite': str})
three = ['S06_7DCG', 'S07_6A0K', 'S08_1J18']
erows = []
for s in three:
    y = SY[SY.system == s].iloc[0]
    erows.append(dict(item='frames, ns', system=s, value=f'{y.frames}, {y.t_last_ns}'))
    say(f'  {s}: {y.frames} frames, {y.t_last_ns} ns')
mc = MC[MC.system.isin(three) & MC.crystal_all_glucoselike.notna() & ~MC.subsite.fillna('').str.contains('~')]
g13m1 = mc[(mc.family == 'GH13') & (mc.subsite == '-1')]
say(f'  GH13 cells at -1 (O1-O6 of the -1 units, 7DCG and 6A0K): {len(g13m1)}; same state {int(g13m1.same_state_all.sum())}')
say(f'  all cells with a crystal state (GH13 and GH14): {len(mc)}; same state {int(mc.same_state_all.sum())}; '
    f'same held/not held {int(mc.same_held_all.sum())}')
o4 = mc[(mc.system == 'S06_7DCG') & (mc.subsite == '+1') & (mc.label == 'O4')].iloc[0]
say(f'  7DCG +1 O4 (scissile): held in {o4.occ_held:.3f} of frames ({o4.main_partners}); crystal state '
    f'{o4.crystal_all_glucoselike}, held in {o4.crystal_all_glucoselike_held:.2f} of groups')
erows += [dict(item='GH13 cells at -1', system='7DCG+6A0K', value=f'{int(g13m1.same_state_all.sum())}/{len(g13m1)}'),
          dict(item='cells same state', system='three', value=f'{int(mc.same_state_all.sum())}/{len(mc)}'),
          dict(item='cells same held/not held', system='three', value=f'{int(mc.same_held_all.sum())}/{len(mc)}'),
          dict(item='7DCG +1 O4 held fraction of frames', system='S06_7DCG', value=f'{o4.occ_held:.3f}')]
pd.DataFrame(erows).to_csv(O / 'E_md_three_systems.tsv', sep='\t', index=False)

(O / 'REPORT.txt').write_text('\n'.join(lines) + '\n')
print('done')
