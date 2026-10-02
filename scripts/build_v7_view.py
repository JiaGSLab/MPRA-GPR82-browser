"""A dedicated Claude v7 browser view, retaining the shared evidence tracks."""
import argparse,csv,json,hashlib,subprocess,urllib.parse,xml.etree.ElementTree as ET,zipfile
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--project-root',type=Path,required=True);ap.add_argument('--bedToBigBed',required=True);args=ap.parse_args()
O=Path(__file__).resolve().parents[1];R=args.project_root
raw='https://raw.githubusercontent.com/JiaGSLab/MPRA-GPR82-browser/main'
snapshot=json.loads((O/'metadata/C1_vs_v7_summary.json').read_text())['v7_pool_sha256']
assert hashlib.sha256((R/'results/v7_2000/oligos.fasta').read_bytes()).hexdigest()==snapshot,'v7 changed; refresh comparison before exporting'
def read(path):return list(csv.DictReader(open(path),delimiter='\t'))
variants={v['variant_id']:v for v in read(R/'results/v7_2000/contrast_readings.tsv')}
oligos=read(R/'results/v7_2000/oligos.tsv')
tracks=json.loads((O/'metadata/track_manifest.json').read_text())
exclude={'design_scope','variants_6000','variants_plan2000','variants_not_retained','c1_only_variants','v7_only_variants','evidence_tiles','survey_tiles','haplotype_existing','windows_6000','windows_plan2000'}
tracks=[t for t in tracks if t['id'] not in exclude]
def add(key,title,desc,rows,vis='pack'):
    rows.sort(key=lambda r:(r[0],r[1],r[2],r[3]));bed=O/'tracks'/f'{key}.bed'
    bed.write_text(''.join('\t'.join(map(str,r))+'\n' for r in rows))
    subprocess.run([args.bedToBigBed,'-type=bed9',str(bed),str(O/'reference/hg38.chrom.sizes'),str(O/'hub/hg38'/f'{key}.bb')],check=True,capture_output=True)
    tracks.append(dict(id=key,title=title,description=desc,count=len(rows),visibility=vis,bed_fields=9,file=f'tracks/{key}.bed'))
regions=json.loads((R/'results/codex_recheck_20261002/region_breakdown.json').read_text())['broad_regions']
rows=[]
for i,r in enumerate(regions):
    a,b=r['start0'],r['end0'];n=sum(a<=int(v['pos0'])<b for v in variants.values())
    rows.append(['chrX',a,b,('R%02d'%(i+1) if i<4 else 'MAIN')+f'|v7_variants={n}',0,'.',a,b,'193,112,43'])
add('v7_design_scope','v7 design scope','Five design intervals, labeled with actual v7 variant counts',rows)
rows=[]
for o in oligos:
    if not o['category'].startswith('archaic_haplotype'):continue
    a,b=int(o['start0']),int(o['end0']);rows.append(['chrX',a,b,o['oligo_id']+'|'+o['category'],0,o['orientation'],a,b,'199,31,160'])
add('v7_haplotype_existing','v7 hap: REF + AB','Actual 170-bp haplotype REF and AB only; no matching A or B singleton in this snapshot',rows)
tracks.sort(key=lambda t:({'v7_design_scope':-2,'evidence_regions':-1}.get(t['id'],0)))
(O/'hub_v7').mkdir(exist_ok=True)
(O/'hub_v7/hub.txt').write_text('hub GPR82_Claude_v7\nshortLabel GPR82 Claude v7\nlongLabel Claude v7 2000 oligos with macrophage ENCODE and archaic-ALT annotations\ngenomesFile genomes.txt\nemail guangshuaijia@gmail.com\n')
(O/'hub_v7/genomes.txt').write_text('genome hg38\ntrackDb trackDb.txt\ndefaultPos chrX:41703000-41751000\n')
stanzas=[]
for i,t in enumerate(tracks):
    stanzas.append(f'track {t["id"]}\nshortLabel {t["title"]}\nlongLabel {t["description"]}\ntype bigBed {t["bed_fields"]}\nbigDataUrl ../hub/hg38/{t["id"]}.bb\nvisibility {t["visibility"]}\nitemRgb on\npriority {i+1}\n')
(O/'hub_v7/trackDb.txt').write_text('\n'.join(stanzas))
for remote in [False,True]:
    root=ET.Element('Session',genome='hg38',locus='chrX:41703000-41751000',version='8');res=ET.SubElement(root,'Resources');panel=ET.SubElement(root,'Panel',name='DataPanel')
    for t in tracks:
        if t['visibility']=='hide':continue
        path=raw+'/'+t['file'] if remote else '../'+t['file']
        ET.SubElement(res,'Resource',path=path,name=t['title']);ET.SubElement(panel,'Track',id=path,name=t['title'],displayMode='EXPANDED' if t['visibility']=='pack' else 'COLLAPSED',height='40',visible='true')
        if not remote:assert (O/'igv'/path).resolve().is_file()
    ET.indent(root);ET.ElementTree(root).write(O/'igv'/f'GPR82_v7_{"remote" if remote else "local"}.xml',encoding='UTF-8',xml_declaration=True)
custom=['browser position chrX:41703000-41751000\n']
for t in tracks:
    custom.append(f'track name="{t["title"]}" description="{t["description"]}" type=bed visibility={t["visibility"]} itemRgb="On"\n'+(O/t['file']).read_text())
(O/'ucsc/GPR82_v7_with_evidence.bed').write_text(''.join(custom))
(O/'metadata/v7_track_manifest.json').write_text(json.dumps(tracks,indent=2)+'\n')
links=json.loads((O/'links.json').read_text());hub=raw+'/hub_v7/hub.txt'
links['v7_hub_url']=hub
for name,pos in [('main','chrX:41703000-41751000'),('overview','chrX:41435000-41755000')]:
    links['v7_ucsc_'+name]='https://genome.ucsc.edu/cgi-bin/hgTracks?'+urllib.parse.urlencode(dict(db='hg38',hubUrl=hub,position=pos))
links['v7_igv_remote']=raw+'/igv/GPR82_v7_remote.xml'
(O/'links.json').write_text(json.dumps(links,indent=2)+'\n')
instructions=f'''Claude v7 专用浏览包：hg38 / 170 bp insert / 2000 oligos

IGV: 解压后 File > Open Session 打开 igv/GPR82_v7_local.xml。请保留 tracks/ 相对目录。
在线会话: {links['v7_igv_remote']}
UCSC 主区: {links['v7_ucsc_main']}
UCSC 全景: {links['v7_ucsc_overview']}
Hub URL: {hub}

自动带入项目选定的 cCRE、5 份巨噬细胞 DNase、巨噬细胞 H3K27ac、单核细胞 H3K27ac 代理、5 个条件 STAT3、GPR82 基因模型与 rE2G 预测连接。
古人类轨道是 Altai/Vindija/Denisova/Chagyrskaya 中通过本项目筛选的 26 个 ALT 位点注释，非原始古 DNA reads/BAM，也不是全基因组古人类轨道。
ENCODE 证据为项目缓存的局部峰/注释，不是自动加载 ENCODE 全部数据，不会随数据库更新自动刷新。峰条带不是测序信号强度曲线。
默认隐藏密集 v7 oligo 窗口和基因组对照；可在 UCSC 调整显示。IGV 可另加载 tracks/claude_v7_windows.bed 和 tracks/claude_v7_controls.bed。
原 v7 的 56 条 ladder 控制区间长度错误在本浏览包中按来源修正到 170 bp。原合成文件没有改动。
首次访问 UCSC 可能需要手动完成人机验证。
'''
(O/'V7_VIEW_README.txt').write_text(instructions)
files=[O/'V7_VIEW_README.txt',O/'igv/GPR82_v7_local.xml',O/'igv/GPR82_v7_remote.xml',O/'ucsc/GPR82_v7_with_evidence.bed',O/'metadata/v7_track_manifest.json']+[O/t['file'] for t in tracks]
files += [p for stem in ['05_v7_distribution','06_v7_annotation'] for p in (O/'figures').glob(stem+'.*')]
files += [O/'metadata/v7_distribution_counts.json',O/'metadata/v7_distribution_counts.csv'] if (O/'metadata/v7_distribution_counts.json').exists() else []
with zipfile.ZipFile(O/'GPR82_v7_IGV_with_annotations.zip','w',zipfile.ZIP_DEFLATED) as z:
    for f in files:z.write(f,'GPR82_v7_browser/'+str(f.relative_to(O)))
print(json.dumps({'tracks':len(tracks),'visible':sum(t['visibility']!='hide' for t in tracks),'v7_ucsc_main':links['v7_ucsc_main']},indent=2))
