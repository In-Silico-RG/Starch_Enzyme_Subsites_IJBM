import json, sys, itertools
import numpy as np, pandas as pd
S=sys.argv[1]
J=json.load(open(S+'/ident.json')); meta=J['meta']; keys=J['keys']; L=J['seqlen']
d=pd.read_csv(S+'/entries2.tsv',sep='\t',keep_default_na=False)
cm=pd.read_csv(S+'/contact_matrix.tsv',sep='\t')
POS={'GH13':['D289','H92','D179','R177','Y51','H288'],'GH31':['D327','H600','D443','D542','W539','W441','I328','F575','I364']}
def single_link(names,M,thr):
    par={n:n for n in names}
    def f(x):
        while par[x]!=x: x=par[x]
        return x
    for i,j in itertools.combinations(range(len(names)),2):
        if M[i][j]>=thr:
            a,b=f(names[i]),f(names[j])
            if a!=b: par[a]=b
    roots={}; return {n:roots.setdefault(f(n),len(roots)+1) for n in names}
for fam in POS:
    o=J['out'][fam]; names=o['names']; M=np.array(o['M']); Ms=np.array(o['Ms']); Lc=np.array(o['L'])
    dd=d[d.family==fam]; npe=dd.groupby('protein_key').size()
    print(f'\n##### {fam}: entries {len(dd)}, protein keys {len(names)}')
    tab=pd.DataFrame([dict(key=k,n_entries=int(npe[k]),pdbs=' '.join(sorted(dd[dd.protein_key==k].pdb)),len=L[k],**meta[k]) for k in names]).sort_values('n_entries',ascending=False)
    print(tab.drop(columns=['kingdom']).to_string(index=False))
    print('entries per protein distribution', dict(sorted(npe.value_counts().items())))
    # close pairs
    print('pairs >= 50% (ident over aligned cols / over shorter seq / aligned cols):')
    for i,j in itertools.combinations(range(len(names)),2):
        if M[i,j]>=0.30: print(f'  {names[i]:12s} {names[j]:12s} {M[i,j]:.3f} {Ms[i,j]:.3f} {Lc[i,j]}')
    off=M[np.triu_indices(len(names),1)]; print('identity distribution over all pairs: median %.3f, min %.3f, max %.3f'%(np.median(off),off.min(),off.max()))
    for t in [0.95,0.90,0.50,0.30]:
        cl=single_link(names,M,t); cls=single_link(names,Ms,t)
        print(f'threshold {t}: clusters (aligned cols) {len(set(cl.values()))}; (shorter seq) {len(set(cls.values()))}')
        # entry-level chain seq check
        E=np.array(o['E']); ce=single_link(o['ents'],E,t); print(f'   entry-level chain-seq clusters {len(set(ce.values()))}')
        if t in (0.50,0.30):
            grp={}
            for k,v in cl.items(): grp.setdefault(v,[]).append(k)
            for v,ks in grp.items():
                print(f'   cluster {v}: {len(ks)} proteins, {int(sum(npe[k] for k in ks))} entries: '+', '.join(f"{k}({meta[k]['org'].split(' (')[0]})" for k in ks))
            # votes
            dk=dd.set_index('pdb').protein_key; ce_=cm[cm.family==fam].copy(); ce_['cl']=ce_.pdb.map(dk).map(cl)
            g=ce_.groupby('cl')[POS[fam]]
            mean=g.mean(); anyc=g.max(); allc=g.min()
            nC=len(mean)
            print(f'   position votes (n clusters = {nC}): mean-fraction-per-cluster / clusters with any entry / clusters with all entries / clusters with majority')
            for q in POS[fam]:
                print(f'     {q}: entries {int(ce_[q].sum())}/{len(ce_)}; cluster-weighted fraction {mean[q].mean():.3f}; any {int(anyc[q].sum())}/{nC}; all {int(allc[q].sum())}/{nC}; majority {int((mean[q]>0.5).sum())}/{nC}')
            # which clusters miss
            for q in POS[fam]:
                miss=ce_[ce_[q]==0]
                if len(miss): print(f'     {q} missed by: '+', '.join(f"{p}(cl{c})" for p,c in zip(miss.pdb,miss.cl)))
    # taxonomy
    t=tab.copy(); print('superkingdom', t.superkingdom.value_counts().to_dict()); print('kingdom', t.kingdom.value_counts().to_dict()); print('phylum', t.phylum.value_counts().to_dict())
