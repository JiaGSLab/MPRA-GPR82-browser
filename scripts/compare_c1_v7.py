"""Compare actual v7 sequences with the frozen C1-2000 ID-selection proposal."""
import csv,json,hashlib,argparse
from pathlib import Path
from collections import Counter,defaultdict
ap=argparse.ArgumentParser();ap.add_argument('--project-root',type=Path,required=True);args=ap.parse_args()
R=args.project_root;OUT=Path(__file__).resolve().parents[1]/'metadata';OUT.mkdir(exist_ok=True)
def read(p):return list(csv.DictReader(open(R/p),delimiter='\t'))
def rev(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def union(xs):
    out=[]
    for a,b in sorted(xs):
        if out and a<=out[-1][1]:out[-1][1]=max(out[-1][1],b)
        else:out.append([a,b])
    return out
def covered(xs,a,b):return sum(d-c for c,d in union([(max(a,c),min(b,d)) for c,d in xs if a<d and c<b]))
o7=read('results/v7_2000/oligos.tsv');o6=read('results/v6_6000/oligos.tsv');oi={o['oligo_id']:o for o in o7}
c7=read('results/v7_2000/contrast_readings.tsv');c6=read('results/v6_6000/contrast_readings.tsv')
v7={x['variant_id']:x for x in c7};v6={x['variant_id']:x for x in c6}
c1=read('results/design_review_20261002/proposal_2000_planned_contrasts_NOT_ORDER.tsv');s1={x['variant_id'] for x in c1};s7=set(v7)
ids1={x['oligo_id'] for x in read('results/design_review_20261002/proposal_2000_existing_ids_NOT_ORDER.tsv')}
tile7=[(int(o['start0']),int(o['end0'])) for o in o7 if o['category'] in ['region_tile','locus_survey_tile']]
tile1=[(int(o['start0']),int(o['end0'])) for o in o6 if o['oligo_id'] in ids1 and o['category'] in ['region_tile','locus_survey_tile']]
broad=json.load(open(R/'results/codex_recheck_20261002/region_breakdown.json'))['broad_regions']
regions=json.load(open(R/'results/design_review_20261002/region_evidence.json'))
regionrows=[]
for rg in regions:
    a,b=rg['start0'],rg['end0'];regionrows.append(dict(region=rg['region'],interval=rg['interval_1based'],
      c1_variants=sum(a<=int(v6[k]['pos0'])<b for k in s1),v7_variants=sum(a<=int(v7[k]['pos0'])<b for k in s7),
      c1_scan_bp=covered(tile1,a,b),v7_scan_bp=covered(tile7,a,b),region_bp=b-a))
coverage1=sum(covered(tile1,r['start0'],r['end0']) for r in broad)
coverage7=sum(covered(tile7,r['start0'],r['end0']) for r in broad)
errors=[]
assert len(o7)==len(oi)==len({o['oligo_sequence'] for o in o7})==2000
assert {len(o['oligo_sequence']) for o in o7}=={200}
assert {len(o['insert_sequence']) for o in o7}=={170}
for r in c7:
    rf=oi[r['ref_oligo_id']]['insert_sequence'];alt=oi[r['alt_oligo_id']]['insert_sequence']
    if r['orientation']=='-':rf,alt=rev(rf),rev(alt)
    idx=int(r['pos0'])-int(r['window_start0'])
    if not (sum(a!=b for a,b in zip(rf,alt))==1 and rf[idx]==r['ref'] and alt[idx]==r['alt']):errors.append(r['variant_id'])
assert not errors
ref=json.load(open(R/'data/processed/v2_variant_scope.json'))['reference38']
for r in c7:
    a,b=int(r['window_start0']),int(r['window_end0']);seq=oi[r['ref_oligo_id']]['insert_sequence']
    if r['orientation']=='-':seq=rev(seq)
    assert seq==ref['sequence'][a-ref['start']:b-ref['start']]
haps=[o for o in o7 if o['category'].startswith('archaic_haplotype')]
rf=next(o for o in haps if o['category'].endswith('REF'));ab=next(o for o in haps if o['category'].endswith('ALT'))
diff=[i for i,(a,b) in enumerate(zip(rf['insert_sequence'],ab['insert_sequence'])) if a!=b]
assert diff==[9,160]
seqmap={o['insert_sequence']:o['oligo_id'] for o in o7};hapstates={}
for state,idxs in [('REF',[]),('A',[diff[0]]),('B',[diff[1]]),('AB',diff)]:
    seq=list(rf['insert_sequence'])
    for idx in idxs:seq[idx]=ab['insert_sequence'][idx]
    seq=''.join(seq);hapstates[state]={'forward_id':seqmap.get(seq),'reverse_id':seqmap.get(rev(seq))}
def count_contexts(cr,vs):
    by=defaultdict(set)
    for r in cr:by[r['variant_id']].add((r['window_start0'],r['window_end0']))
    return dict(Counter(len(by[k]) for k in vs))
manual7={k for k,v in v7.items() if v['tier']=='0'}
arc7={k for k,v in v7.items() if json.loads(v['archaic_ALT_samples'])}
focus7={k for k,v in v7.items() if v['manual_focus_Sheet3']=='True'}
ancestry={k for k,v in v6.items() if float(v.get('gnomAD_AF_grpmax') or 0)>=.01}
ann=json.load(open(R/'results/design_review_20261002/evidence_tracks.json'))
badannotations=[]
for k,v in v7.items():
    pos=int(v['pos0']);dn=sum(any(a<=pos<b for a,b in xs) for source,xs in ann.items() if source in ['ENCFF580ICE','ENCFF174FFX','ENCFF619TQP','ENCFF586HMB','ENCFF948HMV'])
    h3=any(a<=pos<b for a,b in ann['macrophage_H3K27ac'])
    if dn!=len(json.loads(v['macrophage_DNase_sources'])) or h3!=(v['macrophage_H3K27ac']=='True'):
        badannotations.append(dict(variant_id=k,old_dnase_count=len(json.loads(v['macrophage_DNase_sources'])),corrected_dnase_count=dn,old_H3=v['macrophage_H3K27ac'],corrected_H3=h3))
bedbad=[]
for line in (R/'results/v7_2000/deliverables/igv/GPR82_v7_negative_and_activity_controls.bed').read_text().splitlines():
    if line.startswith(('track','browser','#')):continue
    f=line.split('\t')
    if len(f)>3 and int(f[2])-int(f[1])!=170:bedbad.append(f[:4])
summary=dict(c1_version='C1-2000',c1_basis='v6 frozen pool; 1998 existing sequences + 2 undesigned slots',v7_version='Claude-v7-2000-170bp',
    c1_insert_bp=200,v7_insert_bp=170,c1_oligo_nt=230,v7_oligo_nt=200,c1_variants=len(s1),v7_variants=len(s7),
    shared_variants=len(s1&s7),c1_only=len(s1-s7),v7_only=len(s7-s1),v7_variants_outside_original_v6=sorted(s7-set(v6)),
    c1_tiling_oligos=len(tile1),v7_tiling_oligos=len(tile7),c1_covered_bp=coverage1,v7_covered_bp=coverage7,design_scope_bp=sum(r['bp'] for r in broad),
    v7_controls=sum(n for cat,n in Counter(o['category'] for o in o7).items() if 'control' in cat),
    v7_categories=dict(Counter(o['category'] for o in o7)),v7_manual=len(manual7),v7_archaic_ALT=len(arc7),v7_focus=len(focus7),
    c1_manual_genomic_backgrounds=count_contexts(c1,manual7),v7_manual_genomic_backgrounds=count_contexts(c7,manual7),
    c1_ancestry_common_retained=len(s1&ancestry),v7_ancestry_common_retained=len(s7&ancestry),ancestry_common_in_v6=len(ancestry),
    v7_haplotype_states=hapstates,v7_haplotype_flanks=[9,9],v7_annotation_discrepancy_count=len(badannotations),
    v7_export_control_span_errors=len(bedbad),v7_pairs_checked=len(c7),v7_unique_readings=len({(r['variant_id'],r['window_start0'],r['orientation']) for r in c7}),
    v7_pool_sha256=hashlib.sha256((R/'results/v7_2000/oligos.fasta').read_bytes()).hexdigest(),
    v7_summary_reported_oligo_nt=json.load(open(R/'results/v7_2000/design_summary.json'))['unique_complete_length_nt'])
for name,obj in [('C1_vs_v7_summary.json',summary),('C1_vs_v7_regions.json',regionrows),('v7_annotation_discrepancies.json',badannotations),('v7_control_bed_span_errors.json',bedbad)]:
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
with open(OUT/'C1_vs_v7_variant_membership.tsv','w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['variant_id','in_C1_2000','in_Claude_v7','original_tier'])
    for k in sorted(s1|s7):w.writerow([k,k in s1,k in s7,(v6.get(k) or v7[k])['tier']])
print(json.dumps(summary,ensure_ascii=False,indent=2));print('REGIONS',regionrows)
