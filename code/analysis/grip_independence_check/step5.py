import json,sys,itertools,random
import numpy as np,pandas as pd
sys.path.insert(0,sys.argv[1]+'/../code') if False else None
S=sys.argv[1]
J=json.load(open(S+'/ident.json')); meta=J['meta']
d=pd.read_csv(S+'/entries2.tsv',sep='\t',keep_default_na=False)
cm=pd.read_csv(S+'/contact_matrix.tsv',sep='\t')
POS={'GH13':['D289','H92','D179','R177','Y51','H288'],'GH31':['D327','H600','D443','D542','W539','W441','I328','F575','I364']}
TAXN={'PDB:6GXV':['Bacteria','Bacillota'],'PDB:8DL2':['Bacteria','Bacteroidota'],'PDB:1W9X':['Bacteria','Bacillota'],
      'PDB:24ZQ':['Bacteria','Pseudomonadota'],'PDB:7P01':['Eukaryota','Fungi','Ascomycota'],'PDB:7P07':['Eukaryota','Fungi','Ascomycota']}
def lin(k):
    if k in TAXN: return TAXN[k]
    return json.load(open(f'{S}/dl/{k}.json'))['organism']['lineage']
def group(l):
    if l[0]=='Bacteria': return 'Bacteria / '+l[1]
    if l[0]=='Archaea': return 'Archaea / '+l[1]
    if 'Metazoa' in l: return 'Eukaryota / Metazoa / '+l[l.index('Metazoa')+1]
    if 'Fungi' in l: return 'Eukaryota / Fungi / '+[x for x in l if x.endswith('mycota')][0]
    if 'Viridiplantae' in l: return 'Eukaryota / Viridiplantae / Streptophyta'
    return 'Eukaryota / '+l[1]
def sl(names,M,t):
    par={n:n for n in names}
    def f(x):
        while par[x]!=x: x=par[x]
        return x
    for i,j in itertools.combinations(range(len(names)),2):
        if M[i][j]>=t:
            a,b=f(names[i]),f(names[j])
            if a!=b: par[a]=b
    r={}; return {n:r.setdefault(f(n),len(r)+1) for n in names}
res={}
for fam in POS:
    o=J['out'][fam]; names=o['names']; M=np.array(o['M']); Ms=np.array(o['Ms'])
    dd=d[d.family==fam]; dk=dd.set_index('pdb').protein_key
    tx=pd.Series({k:group(lin(k)) for k in names})
    print(fam,'taxonomy (protein keys):',tx.value_counts().to_dict())
    print(fam,'lineage[0] per entries:',dd.protein_key.map(lambda k: lin(k)[0]).value_counts().to_dict())
    for t,mat,lab in [(0.5,M,'aligned'),(0.3,Ms,'shorter'),(0.4,M,'aligned'),(0.4,Ms,'shorter')]:
        cl=sl(names,mat,t); grp={}
        for k,v in cl.items(): grp.setdefault(v,[]).append(k)
        c=cm[cm.family==fam].copy(); c['cl']=c.pdb.map(dk).map(cl); g=c.groupby('cl')[POS[fam]]
        mean,anyc,allc=g.mean(),g.max(),g.min(); n=len(mean)
        print(f' {lab} {t}: {n} clusters; sizes(entries)', sorted(c.cl.value_counts().tolist(),reverse=True))
        if t==0.3:
            for v,ks in grp.items(): print('   ',v,len(ks),int(sum((dk==k).sum() for k in ks)),', '.join(ks))
        for q in POS[fam]:
            print(f'   {q}: entries {int(c[q].sum())}/{len(c)}; clusters any {int(anyc[q].sum())}/{n}; all {int(allc[q].sum())}/{n}; mean {mean[q].mean():.3f}')
# null: shuffled sequences
from Bio import Align
from Bio.Align import substitution_matrices
from Bio import SeqIO
al=Align.PairwiseAligner(); al.mode='global'; al.substitution_matrix=substitution_matrices.load('BLOSUM62')
al.open_gap_score=-10; al.extend_gap_score=-0.5; al.end_gap_score=0
random.seed(1)
seqs=[]
for k in ['P04746','P07206','Q9FYY0','P38940','O06994','P53051','P00693','Q7LYT7','L0N7E5','P10253']:
    s=json.load(open(f'{S}/dl/{k}.json'))['sequence']['value']; seqs.append(s[:700])
vals=[];vs=[]
for i in range(12):
    a,b=random.sample(seqs,2); b=''.join(random.sample(b,len(b)))
    aln=al.align(a,b)[0]; cols=[(x,y) for x,y in zip(aln[0],aln[1]) if x!='-' and y!='-']
    nid=sum(x==y for x,y in cols); vals.append(nid/len(cols)); vs.append(nid/min(len(a),len(b)))
print('NULL shuffled: aligned-col identity mean %.3f range %.3f-%.3f; shorter-seq mean %.3f range %.3f-%.3f'%(np.mean(vals),min(vals),max(vals),np.mean(vs),min(vs),max(vs)))
