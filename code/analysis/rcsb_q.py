import json, urllib.request, re
URL='https://search.rcsb.org/rcsbsearch/v2/query'
def q(body):
    req=urllib.request.Request(URL,data=json.dumps(body).encode(),headers={'Content-Type':'application/json'})
    with urllib.request.urlopen(req,timeout=60) as r:
        d=r.read()
    return json.loads(d) if d.strip() else {'result_set':[],'total_count':0}
def T(attr,op,val,svc='text_chem'): return {"type":"terminal","service":svc,"parameters":{"attribute":attr,"operator":op,"value":val}}
subs=['methyl','acetyl','hydroxypropyl','carboxymethyl','succinyl','octenyl','sulfo','phosphono','ethyl','propyl','benzyl','caffeoyl','coumaroyl','feruloyl','galloyl']
comps={}
for g in ['glucopyranose','glucopyranoside','glucopyranosyl','glucose']:
    for s in subs:
        body={"query":{"type":"group","logical_operator":"and","nodes":[T("chem_comp.name","contains_words",g),T("chem_comp.name","contains_words",s)]},
              "return_type":"mol_definition","request_options":{"return_all_hits":True}}
        try: j=q(body)
        except Exception as e: print('ERR',g,s,e); continue
        for x in j.get('result_set',[]): comps.setdefault(x['identifier'],set()).add(s)
print('comps',len(comps))
ids=sorted(c for c in comps if not c.startswith('PRD'))
pf=["PF00128","PF01055","PF01373","PF00723","PF02806","PF16863"]
fam_node={"type":"group","logical_operator":"or","nodes":[T("rcsb_polymer_entity_annotation.annotation_id","exact_match",p,'text') for p in pf]}
chem_node={"type":"group","logical_operator":"or","nodes":[T("rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id","in",ids,'text'),
                                                       T("rcsb_branched_entity_container_identifiers.chem_comp_monomers","in",ids,'text')]}
body={"query":{"type":"group","logical_operator":"and","nodes":[fam_node,chem_node]},"return_type":"entry","request_options":{"return_all_hits":True}}
j=q(body)
ents=[x['identifier'] for x in j.get('result_set',[])]
print('entries',len(ents),ents)
json.dump(dict(comps={k:sorted(v) for k,v in comps.items()},entries=ents),open('rcsb_hits.json','w'),indent=1)
