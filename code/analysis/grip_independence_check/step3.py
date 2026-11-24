"""Per entry x position contact matrix, rebuilt exactly as a33_subsites.py (res_entries)."""
import sys, pandas as pd
P=__import__('os').path.abspath(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', '..')) + '/data/analysis/33_subsites/'
S=sys.argv[1]
u=pd.read_csv(P+'units.tsv',sep='\t',dtype={'subsite':str}); u['subsite']=u.subsite.fillna('')
rep=u[(u.subsite=='-1')&(u.rscc>=0.80)&u.in_tables&u.family.isin(['GH13','GH31'])]
c=pd.read_csv(P+'contacts_mapped.tsv',sep='\t',dtype={'rings':str,'ref_position':str})
c['rings']=c.rings.fillna(''); c['ref_position']=c.ref_position.fillna('')
good=c[c.asu==True]
POS={'GH13':['D289','H92','D179','R177','Y51','H288'],'GH31':['D327','H600','D443','D542','W539','W441','I328','F575','I364']}
rows=[]
for f in POS:
    r=rep[rep.family==f]
    hit={}
    for x in r.itertuples():
        g=good[(good.instance_id==x.instance_id)]
        g=g[g.rings.apply(lambda s: str(x.ring) in s.split(','))]
        g=g[(g.ref_position!='')&(g.same_chain_as_frame==True)]
        hit.setdefault(x.pdb,set()).update(g.ref_position)
    for p,h in hit.items():
        rows.append(dict(family=f,pdb=p,**{q:int(q in h) for q in POS[f]}))
m=pd.DataFrame(rows); m.to_csv(S+'/contact_matrix.tsv',sep='\t',index=False)
for f in POS:
    mm=m[m.family==f]; print(f,len(mm),{q:int(mm[q].sum()) for q in POS[f]})
