"""Export hg38 browser annotations; does not modify or export synthesis sequences."""
import argparse,csv,json,re,hashlib,shutil,subprocess,urllib.parse,zipfile
from pathlib import Path
from collections import Counter
import xml.etree.ElementTree as ET

ap=argparse.ArgumentParser();ap.add_argument('--project-root',type=Path,required=True)
ap.add_argument('--bedToBigBed',type=Path,required=True);args=ap.parse_args()
ROOT=args.project_root.resolve();OUT=Path(__file__).resolve().parents[1]
REPO='https://github.com/JiaGSLab/MPRA-GPR82-browser'
RAW='https://raw.githubusercontent.com/JiaGSLab/MPRA-GPR82-browser/main'
for sub in ['tracks','hub/hg38','igv','ucsc','figures','metadata']:(OUT/sub).mkdir(parents=True,exist_ok=True)
def load(p):return json.loads((ROOT/p).read_text())
def read(p):return list(csv.DictReader(open(ROOT/p),delimiter='\t'))
def clean(s):return re.sub(r'[^A-Za-z0-9_:.+|=,;()/~-]','_',str(s))[:250]
def bed(a,b,name,color='40,110,159',chrom='chrX',strand='.',score=0):
    return [chrom,int(a),int(b),clean(name),int(score),strand,int(a),int(b),color]
def savejson(p,obj):(OUT/p).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
review='results/design_review_20261002/'
e=load(review+'region_evidence.json');evtracks=load(review+'evidence_tracks.json')
locus=load('data/processed/locus.json');summary=load(review+'proposal_2000_summary.json')
broad=load('results/codex_recheck_20261002/region_breakdown.json')['broad_regions']
oligos=read('results/v6_6000/oligos.tsv')
variants={r['variant_id']:r for r in read('results/v6_6000/contrast_readings.tsv')}
planned=read(review+'proposal_2000_planned_contrasts_NOT_ORDER.tsv')
keptvar={r['variant_id'] for r in planned}
keepids={r['oligo_id'] for r in read(review+'proposal_2000_existing_ids_NOT_ORDER.tsv')}
assert len(variants)==1574 and len(keptvar)==691 and len(keepids)==1998
stat3=read('results/stat3/STAT3_peaks_in_locus_hg38.tsv')
tracks=[]
def add(key,title,description,rows,vis='dense',fields=9):
    rows=sorted(rows,key=lambda r:(r[0],int(r[1]),int(r[2]),r[3]))
    if not rows:return
    path=OUT/'tracks'/f'{key}.bed'
    path.write_text(''.join('\t'.join(map(str,r))+'\n' for r in rows))
    tracks.append(dict(id=key,title=title,description=description,count=len(rows),visibility=vis,bed_fields=fields,file=f'tracks/{key}.bed'))

add('design_scope','Design scope','Five approved design intervals; labels give unique current oligo counts',
    [bed(r['start0'],r['end0'],f'{"R%02d"%(i+1) if i<4 else "MAIN"}|oligos={r["unique_oligos"]}|variants={r["selected_variants"]}') for i,r in enumerate(broad)],'pack')
add('evidence_regions','13 evidence regions','Merged evidence regions; region-wide marks do not support every base',
    [bed(r['start0'],r['end0'],f'{r["region"]}|DNase={r["dnase_source_count"]}|H3bp={r["macrophage_H3K27ac_bp"]}|rE2G={r["max_rE2G"] if r["max_rE2G"] else "NA"}',
         '133,99,166' if r['region'] in ['R06','R09'] else '40,110,159') for r in e],'pack')
add('gene_spans','Nearby gene spans','Cached Ensembl gene spans; these are NOT exon models',
    [bed(g['start']-1,g['end'],g.get('external_name') or g['id'],'117,132,147',strand='+' if g['strand']==1 else '-') for g in locus['neighbors']],'pack')
models=[]
for t in locus['gene']['Transcript']:
    ex=sorted(t.get('Exon',[]),key=lambda x:x['start'])
    if not ex:continue
    a=t['start']-1;z=t['end'];tr=t.get('Translation')
    row=bed(a,z,t['id']+'.'+str(t['version'])+'|'+t.get('display_name','GPR82'),'33,133,121',strand='+')
    row[6:8]=[tr['start']-1,tr['end']] if tr else [a,a]
    models.append(row+[len(ex),','.join(str(x['end']-x['start']+1) for x in ex)+',',','.join(str(x['start']-1-a) for x in ex)+','])
add('gpr82_transcripts','GPR82 transcripts','Cached Ensembl exon models; thick blocks are translated coding sequence',models,'pack',12)
add('gpr82_promoter','Design promoter','Design TSS=chrX:41724175; promoter is -2000 to +500 bp',
    [bed(41722174,41724674,'GPR82_promoter_-2000_+500','204,115,51',strand='+')],'pack')
add('re2g','rE2G to GPR82','Predicted linked element intervals; scores are not proof of causality',
    [bed(r['start'],r['end'],'GPR82_link|score='+r['Score'],'133,99,166',score=round(float(r['Score'])*1000)) for r in locus['links']],'pack')
add('ccre','Cached cCREs','ENCODE cCRE annotation; cache covers chrX:41624174-41830130 only; distal absence is not negative evidence',
    [bed(r['chromStart'],r['chromEnd'],r['name']+'|'+r['cCRE_class'],'184,140,58') for r in locus['ccres']],'dense')
dnase_names={'ENCFF580ICE':'inflammatory_male21_baseline','ENCFF174FFX':'inflammatory_male40_baseline',
    'ENCFF619TQP':'suppressor_male21_baseline','ENCFF586HMB':'suppressor_male40_baseline','ENCFF948HMV':'inflammatory_male40_LPS24h'}
for source,condition in dnase_names.items():
    add('dnase_'+source,source+' DNase',condition+'; called peak intervals, not signal intensity',
        [bed(a,z,source+'|'+condition,'43,143,99') for a,z in evtracks[source]])
for key,label,color in [('macrophage_H3K27ac','Macrophage H3K27ac','204,115,51'),('ENCFF560PQQ','Monocyte H3 proxy','218,169,89')]:
    add('h3_'+key,label,'GSE80727 control/IL10 union, reciprocal hg19-to-hg38 lift' if key.startswith('mac') else 'ENCFF560PQQ monocyte proxy; not macrophage',
        [bed(a,z,label,color) for a,z in evtracks[key]])
for condition in sorted({r['condition'] for r in stat3}):
    add('stat3_'+clean(condition),'STAT3 '+condition,'GSE120943 via ReMap2022; '+condition+' called peaks, not differential binding',
        [bed(r['start0'],r['end0'],'STAT3|'+condition+'|score='+r['raw_bed_score'],'156,87,142') for r in stat3 if r['condition']==condition])
tiercolor={'0':'214,39,40','1':'148,103,189','2':'31,78,160','3':'44,127,184','4':'127,180,214','5':'160,160,160'}
def vbed(v,color=None):
    name=v['variant_id']+'|tier='+v['tier']
    rs=json.loads(v['rsids']);arc=json.loads(v['archaic_ALT_samples'])
    if rs:name+='|'+','.join(rs)
    if v['manual_source_cell']:name+='|manual='+v['manual_source_cell']
    if arc:name+='|archaic='+','.join(arc)
    return bed(v['pos0'],int(v['pos0'])+1,name,color or tiercolor[v['tier']])
add('variants_6000','6000: 1574 SNPs','Current library distinct substitutions; colours show original tier',[vbed(v) for v in variants.values()])
add('manual_38','Manual: 38 SNPs','All 38 currently included user-specified variants',[vbed(v,'214,39,40') for v in variants.values() if v['tier']=='0'],'pack')
add('archaic_26','Archaic ALT: 26','ALT relative to hg38 in archaic samples; not proven lineage-specific or introgressed',[vbed(v,'148,103,189') for v in variants.values() if json.loads(v['archaic_ALT_samples'])],'pack')
add('sheet3_6','Sheet3 focus: 6','Six manually prioritized Sheet3 variants',[vbed(v,'230,137,25') for v in variants.values() if v['manual_focus_Sheet3']=='True'],'pack')
add('variants_plan2000','C1-2000: 691 SNPs','Codex C1-2000 planning proposal only; 691 variants, not synthesis-ready',[vbed(v,'33,133,121') for k,v in variants.items() if k in keptvar])
add('variants_not_retained','Not retained: 883','Current variants absent from the 2000-slot proposal',[vbed(v,'180,180,180') for k,v in variants.items() if k not in keptvar],'hide')
v7={r['variant_id']:r for r in read('results/v7_2000/contrast_readings.tsv')}
o7=read('results/v7_2000/oligos.tsv')
add('claude_v7_variants','Claude v7: 806 SNPs','Claude v7 actual 2000 oligos, 170-bp inserts; compare with C1-2000 proposal',
    [vbed(v,'193,112,43') for v in v7.values()])
add('c1_only_variants','C1-only SNPs','Variants in C1-2000 but absent from Claude v7',
    [vbed(variants[k],'33,133,121') for k in keptvar-set(v7)],'hide')
add('v7_only_variants','v7-only SNPs','Variants in Claude v7 but absent from C1-2000',
    [vbed(v7[k],'193,112,43') for k in set(v7)-keptvar],'hide')
add('claude_v7_tiles','Claude v7 tiles','498 actual 170-bp scanning inserts in Claude v7',
    [bed(o['start0'],o['end0'],o['oligo_id'],'193,112,43',strand=o['orientation']) for o in o7 if o['category'] in ['region_tile','locus_survey_tile']])
# Use the v7 ladder source coordinates, never its legacy 200-bp BED export.
v7coordkeys={}
for rows in load('data/processed/v7_control_ladder.json').values():
    for r in rows if isinstance(rows,list) else [rows]:
        if not r.get('window') or r.get('offset') is None:continue
        chrom,span=r['window'].split(':');a=int(span.split('-')[0])+int(r['offset'])
        v7coordkeys[(r.get('gene') or r.get('label'),r.get('piece'))]=(chrom,a,a+170)
v7coords={}
for r in load('data/processed/v7_intermediate.json')['control_ladder_rows']:
    coord=v7coordkeys.get((r.get('name') or r.get('label'),r.get('piece')))
    if coord:
        for oid in r.get('oligo_ids',[]):v7coords[oid]=coord
v7windows=[];v7controls=[]
for o in o7:
    coord=(o['chrom'],int(o['start0']),int(o['end0'])) if o['chrom'] and o['start0'] else v7coords.get(o['oligo_id'])
    if not coord:continue
    chrom,a,z=coord;assert z-a==170
    row=bed(a,z,o['oligo_id']+'|'+o['category'],'193,112,43',chrom,o['orientation'])
    v7windows.append(row)
    if 'control' in o['category']:v7controls.append(row)
add('claude_v7_windows','v7: 170bp windows','Actual v7 insert footprints; 56 ladder-control ends corrected from source to 170 bp',v7windows,'hide')
add('claude_v7_controls','v7 genomic controls','Genomic controls only; corrected 170-bp ladder footprints, shuffles and vectors unplaced',v7controls,'hide')
for cat,key,label,color in [('region_tile','evidence_tiles','232 evidence tiles','40,110,159'),('locus_survey_tile','survey_tiles','346 survey tiles','139,182,170')]:
    add(key,label,'Existing 200-bp insert genomic footprints; preserved in the proposal',
        [bed(o['start0'],o['end0'],o['oligo_id'],color,strand=o['orientation']) for o in oligos if o['category']==cat])
add('haplotype_existing','Haplotype REF+AB','Two existing states only; the two proposed A/B single mutants are NOT yet designed',
    [bed(o['start0'],o['end0'],o['oligo_id']+'|'+o['category'],'199,31,160',strand=o['orientation']) for o in oligos if o['category'].startswith('archaic_haplotype')],'pack')
# Recover the 56 genomic ladder controls from the already verified export.
coords={}
for line in (ROOT/'results/v6_6000/deliverables/igv/GPR82_v6_negative_and_activity_controls.bed').read_text().splitlines():
    if line.startswith(('track','browser','#')):continue
    f=line.split('\t');oid=re.search(r'GPR82v6_\d+',f[3]) if len(f)>3 else None
    if oid:coords[oid.group()]=(f[0],int(f[1]),int(f[2]))
mapped={}
for o in oligos:
    if o['chrom'] and o['start0']:mapped[o['oligo_id']]=(o['chrom'],int(o['start0']),int(o['end0']))
    elif o['oligo_id'] in coords:mapped[o['oligo_id']]=coords[o['oligo_id']]
for key,label,ids,color in [('windows_6000','6000 insert windows',{o['oligo_id'] for o in oligos},'85,114,151'),('windows_plan2000','C1-2000 windows',keepids,'33,133,121')]:
    rows=[]
    for o in oligos:
        oid=o['oligo_id']
        if oid in ids and oid in mapped:
            c,a,z=mapped[oid];rows.append(bed(a,z,oid+'|'+o['category'],color,c,o['orientation']))
    add(key,label,'200-bp genomic insert footprints, NOT 230-nt oligos; adapters have no genomic coordinates',rows,'hide')

# Validate coordinate bounds and BED structure before any publication.
sizes={a:int(z) for a,z in [x.split() for x in (OUT/'reference/hg38.chrom.sizes').read_text().splitlines()]}
for t in tracks:
    lines=(OUT/t['file']).read_text().splitlines()
    for line in lines:
        f=line.split('\t');assert len(f)==t['bed_fields']
        a,z=int(f[1]),int(f[2]);assert 0<=a<z<=sizes[f[0]]
        assert a<=int(f[6])<=int(f[7])<=z and len(f[3])<=255
        if t['bed_fields']==12:
            starts=[int(x) for x in f[11].strip(',').split(',')];lengths=[int(x) for x in f[10].strip(',').split(',')]
            assert len(starts)==len(lengths)==int(f[9]);assert all(0<=s and s+n<=z-a for s,n in zip(starts,lengths))
    if args.bedToBigBed:
        proc=subprocess.run([str(args.bedToBigBed),f'-type=bed{t["bed_fields"]}',str(OUT/t['file']),str(OUT/'reference/hg38.chrom.sizes'),str(OUT/'hub/hg38'/f'{t["id"]}.bb')],capture_output=True,text=True)
        if proc.returncode:raise RuntimeError(t['id']+': '+proc.stderr)

(OUT/'hub/hub.txt').write_text('hub GPR82_MPRA\nshortLabel GPR82 MPRA\nlongLabel GPR82 hg38: v6, Codex C1-2000 proposal, Claude v7\ngenomesFile genomes.txt\nemail guangshuaijia@gmail.com\n')
(OUT/'hub/genomes.txt').write_text('genome hg38\ntrackDb hg38/trackDb.txt\ndefaultPos chrX:41703000-41751000\n')
stanzas=[]
for i,t in enumerate(tracks):
    stanzas.append(f'track {t["id"]}\nshortLabel {t["title"]}\nlongLabel {t["description"]}\ntype bigBed {t["bed_fields"]}\nbigDataUrl {t["id"]}.bb\nvisibility {t["visibility"]}\nitemRgb on\npriority {i+1}\n')
(OUT/'hub/hg38/trackDb.txt').write_text('\n'.join(stanzas))
custom=['browser position chrX:41703000-41751000\n']
for t in tracks:
    custom.append(f'track name="{t["title"]}" description="{t["description"]}" type=bed {t["bed_fields"]} visibility={t["visibility"]} itemRgb="On"\n')
    custom.append((OUT/t['file']).read_text())
(OUT/'ucsc/GPR82_all_custom_tracks.bed').write_text(''.join(custom))

def session(path,remote=False,overview=False):
    root=ET.Element('Session',genome='hg38',locus='chrX:41435000-41755000' if overview else 'chrX:41703000-41751000',version='8')
    resources=ET.SubElement(root,'Resources');panel=ET.SubElement(root,'Panel',name='DataPanel')
    for t in tracks:
        if t['visibility']=='hide':continue
        ref=RAW+'/'+t['file'] if remote else '../'+t['file']
        ET.SubElement(resources,'Resource',path=ref,name=t['title'])
        ET.SubElement(panel,'Track',id=ref,name=t['title'],displayMode='EXPANDED' if t['visibility']=='pack' else 'COLLAPSED',height='45' if t['visibility']=='pack' else '25',visible='true')
    ET.indent(root);ET.ElementTree(root).write(path,encoding='UTF-8',xml_declaration=True)
for remote in [False,True]:
    for overview in [False,True]:session(OUT/'igv'/f'GPR82_{"remote" if remote else "local"}_{"overview" if overview else "main"}.xml',remote,overview)
huburl=RAW+'/hub/hub.txt'
def link(pos):return 'https://genome.ucsc.edu/cgi-bin/hgTracks?'+urllib.parse.urlencode(dict(db='hg38',hubUrl=huburl,position=pos))
links={'github':REPO,'hub_url':huburl,'ucsc_main':link('chrX:41703000-41751000'),'ucsc_overview':link('chrX:41435000-41755000'),
       'igv_remote_main':RAW+'/igv/GPR82_remote_main.xml','ucsc_custom_fallback':'https://genome.ucsc.edu/cgi-bin/hgTracks?'+urllib.parse.urlencode(dict(db='hg38',position='chrX:41703000-41751000',hgct_customText=RAW+'/ucsc/GPR82_all_custom_tracks.bed'))}
links['regions']={r['region']:link(f'chrX:{r["start0"]-500}-{r["end0"]+500}') for r in e}
savejson('links.json',links);savejson('metadata/track_manifest.json',tracks)
savejson('metadata/export_summary.json',dict(assembly='hg38',track_count=len(tracks),current_variants=len(variants),planned_variants=len(keptvar),current_oligos=6000,current_genomically_mapped=len(mapped),current_unmapped=6000-len(mapped),proposed_existing=1998,proposed_existing_mapped=len(keepids&set(mapped)),reserved_undesigned=2))
for fn in ['region_evidence.tsv','region_comparison.tsv','allocation_comparison.tsv','proposal_2000_summary.json']:
    shutil.copyfile(ROOT/review/fn,OUT/'metadata'/fn)
for fn in ['01_region_distribution.png','01_region_distribution.svg','02_allocation_6000_vs_2000.png','03_region_evidence.png']:
    shutil.copyfile(ROOT/review/fn,OUT/'figures'/fn)
shutil.copyfile(ROOT/(review+'设计评价与2000条方案.txt'),OUT/'metadata/design_review_zh.txt')
savejson('metadata/input_hashes.json',{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'results/v6_6000/oligos.tsv',ROOT/'results/v6_6000/contrast_readings.tsv',ROOT/'results/v7_2000/oligos.tsv',ROOT/'results/v7_2000/contrast_readings.tsv',ROOT/(review+'evidence_tracks.json'),ROOT/'data/processed/locus.json']})
readme=f'''# GPR82 MPRA genome browser tracks (hg38)

包含 **Claude v6（6000 条）**、**Codex C1-2000（2000 条预算方案）**及 **Claude v7（2000 条实际设计）**。C1-2000 是此前称作 PLAN2000 的方案：1998 条现有序列加 2 个尚未设计的单倍型单突变名额，不是合成订单。v7 则是 2000 条、每条总长 200 nt、插入片段 170 bp 的实际文件。

C1 与 v7 的逐位点、逐区比较见 [比较报告](metadata/C1_vs_v7_comparison_zh.txt) 与 [比较数据](metadata/C1_vs_v7_summary.json)。

## 一键打开 UCSC

- [主区：GPR82 与附近证据]({links['ucsc_main']})
- [全景：四个远端区域 + 主区]({links['ucsc_overview']})
- [备用：加载 BED custom tracks]({links['ucsc_custom_fallback']})
- Hub URL: `{huburl}`

可在 UCSC 的 My Data → Track Hubs → My Hubs / Connected Hubs 输入 Hub URL。轨道下方可调整 dense/pack/full，灰色低优先级或密集窗口轨道默认隐藏，可按需打开。

## IGV

参考基因组选择 **Human hg38**。在线使用 File → Open Session，选择 URL 并粘贴：

`{links['igv_remote_main']}`

如果当前 IGV 版本只提供本地会话文件，下载 `igv/GPR82_remote_main.xml` 再打开，它会从 GitHub 加载轨道。
离线使用：下载整个仓库或 `GPR82_IGV_UCSC_bundle.zip` 并解压，再 File → Open Session 打开 `igv/GPR82_local_main.xml`；保持 igv/ 与 tracks/ 的相对目录。
`*_overview.xml` 从全景开始。基因组参考序列可能仍需要联网或已有 hg38 缓存。

## 图

![区域与条数](figures/01_region_distribution.png)
![区域注释与基因组轨道](figures/04_genome_annotation.png)
![6000与2000](figures/02_allocation_6000_vs_2000.png)

## 轨道和解释

{len(tracks)} 个轨道包含：五块设计范围、13 个证据区、邻近基因跨度、GPR82 外显子/CDS、启动子、rE2G 连接区间、cCRE、五份 DNase、巨噬 H3K27ac、单核细胞代理、五个 STAT3 条件、1574 个当前变异、38 个手工位点、26 个古人类 ALT 位点、6 个 Sheet3 重点、691 个拟保留变异、883 个未保留变异、232+346 条扫描、已有单倍型状态及两版基因组窗口。

- 所有 BED 为 **hg38，0-based 半开区间**；浏览器显示为 1-based。仅使用 GRCh38/hg38。
- 窗口表示 **200 bp insert**，不把 30 nt 接头伪装为基因组坐标。
- `gene_spans` 是基因跨度；`gpr82_transcripts` 才是真实缓存外显子模型。设计 TSS 为 41,724,175；canonical 转录本的缓存起点为 41,724,181。
- 证据轨道是已调用的 **峰区间**，不是测序信号强度；不同来源/条件不等于独立重复。
- R01/R02 的 cCRE 不在缓存查询范围内；未显示不是阴性。DNase/H3K27ac 来自复核后的原始峰求交。
- rE2G 是预测，MPRA 活性不能单独证明内源 GPR82 靶向。R09 启动子的证据不能外推到整个基因体。
- 古人类 ALT 是相对于 hg38 的样本非参考等位基因，不等于古人类特有或已证明渗入。
- 未输出没有基因组坐标的打乱和载体序列；具体数目见 metadata/export_summary.json。两个预留的新单突变没有冒充已完成轨道。
- 未将旧 motif “gain/loss”当作已验证差异结合证据展示。

## 区域跳转

'''
for r in e:readme+=f'- [{r["region"]}：{r["interval_1based"]}]({links["regions"][r["region"]]})\n'
readme+='''
## 来源与复现

本地项目快照日期 2026-10-02；数据来源：Ensembl、ENCODE、GSE80727、GSE120943/ReMap2022 和本项目设计。来源与限制见 metadata/，逐轨道条数见 metadata/track_manifest.json。
生成器 `scripts/build_browser_bundle.py` 需要原 MPRA_GPR82 项目（大体积原始文件及合成序列未包含在此浏览仓库）。使用 UCSC bedToBigBed 生成，pybigtools 独立读取验证二进制轨道；使用 scripts/plot_annotation.py 从仓库 BED 绘图。

UCSC 官方文档：[Track Hubs](https://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html) · [Custom Tracks](https://genome.ucsc.edu/goldenPath/help/customTrack.html)
'''
(OUT/'README.md').write_text(readme)
print(json.dumps({'tracks':len(tracks),'current_mapped':len(mapped),'planned_existing_mapped':len(keepids&set(mapped)),'links':links},ensure_ascii=False,indent=2))
