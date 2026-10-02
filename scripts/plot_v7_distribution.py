"""Plot actual Claude v7 primary oligo counts and exported hg38 coordinates."""
from pathlib import Path
from collections import Counter
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
O=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':['PingFang SC','Arial Unicode MS','DejaVu Sans'],'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,'svg.fonttype':'path'})
def read(k):return [l.split('\t') for l in (O/'tracks'/f'{k}.bed').read_text().splitlines()]
regions=read('v7_design_scope');windows=read('claude_v7_windows');summary=json.loads((O/'metadata/C1_vs_v7_summary.json').read_text())
colors=['#277F8E','#71B6AF','#BC5666','#8969A5','#DBA248','#7294B6','#CDD9E3']
names=['证据区扫描','背景区域扫描','手工＋古人类 SNV','其他 SNV','单倍型','基因组对照','无基因组坐标对照']
def group(cat):
 if cat=='region_tile':return 0
 if cat=='locus_survey_tile':return 1
 if cat.startswith(('manual_SNV','archaic_SNV')):return 2
 if 'haplotype' in cat:return 4
 if 'control' in cat:return 5
 return 3
counts=np.zeros((6,7),dtype=int);assigned={};seen=set()
for r in windows:
 oid,cat=r[3].split('|')[:2];assert oid not in seen;seen.add(oid)
 if 'control' in cat:counts[5,5]+=1;continue
 mid=(int(r[1])+int(r[2]))/2
 hits=[i for i,s in enumerate(regions) if r[0]==s[0] and int(s[1])<=mid<int(s[2])]
 assert len(hits)==1,(oid,hits)
 i=hits[0];counts[i,group(cat)]+=1;assigned[oid]=i
counts[5,6]=summary['v7_controls']-counts[5,5]
assert counts.sum()==2000 and counts[:5].sum()==1788
assert counts[:,0].sum()==250 and counts[:,1].sum()==248
labels=['R01','R02','R03','R04','主区 MAIN','独立对照']
metadata={'version':'Claude-v7-2000-170bp','pool_sha256':summary['v7_pool_sha256'],'unit':'unique oligo_id; primary category; controls counted separately','region_assignment':'midpoint of exported genomic window; every non-control maps to exactly one broad design interval','categories':names,'regions':[dict(region=labels[i],oligos=int(counts[i].sum()),category_counts=dict(zip(names,map(int,counts[i])))) for i in range(6)],'total':int(counts.sum())}
(O/'metadata/v7_distribution_counts.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n')
with (O/'metadata/v7_distribution_counts.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['region','total']+names)
 for i in range(6):w.writerow([labels[i],counts[i].sum()]+list(counts[i]))
ink='#263D50';grey='#62798A'
fig=plt.figure(figsize=(15,12));gs=fig.add_gridspec(3,2,height_ratios=[1,1.9,1.55],hspace=.55,wspace=.36)
fig.suptitle('Claude v7｜2000 条 oligo 的区域与功能分布',x=.07,y=.98,ha='left',fontsize=23,color=ink)
fig.text(.07,.925,'实际设计快照  ·  hg38  ·  170 bp 插入片段 + 30 nt 接头  ·  按唯一 oligo ID 计数',fontsize=11,color=grey)
ax=fig.add_subplot(gs[0,:]);ax.set_title('A  全景位置与区域总数（对照另外统计）',loc='left',fontweight='bold',pad=12)
for i,r in enumerate(regions):
 a,b=int(r[1])/1e6,int(r[2])/1e6
 ax.broken_barh([(a,b-a)],(.2,.15),facecolors=colors[0])
 ax.annotate(f'{labels[i]}\n{counts[i].sum():,} 条',((a+b)/2,.36),((a+b)/2,.67 if i!=3 else .86),ha='center',fontsize=10,arrowprops={'arrowstyle':'-','lw':.7,'color':grey})
ax.set_xlim(41.435,41.755);ax.set_ylim(.05,1.1);ax.set_yticks([]);ax.set_xlabel('chrX 坐标（Mb）');ax.spines['left'].set_visible(False)
ax=fig.add_subplot(gs[1,0]);ax.set_title('B  各区域的 oligo 分配',loc='left',fontweight='bold',pad=18)
left=np.zeros(6)
for j in range(7):ax.barh(np.arange(6),counts[:,j],left=left,color=colors[j],height=.66,label=names[j]);left+=counts[:,j]
for i,n in enumerate(counts.sum(axis=1)):ax.text(n+19,i,f'{n:,}',va='center',fontweight='bold',color=ink)
ax.set_yticks(np.arange(6),labels);ax.invert_yaxis();ax.set_xlim(0,1795);ax.set_xlabel('唯一 oligo 数');ax.grid(axis='x',alpha=.12);ax.set_axisbelow(True)
ax=fig.add_subplot(gs[1,1]);ax.set_title('C  全库功能构成（合计 2000 条）',loc='left',fontweight='bold',pad=18)
y=np.arange(7);tot=counts.sum(axis=0);ax.barh(y,tot,color=colors,height=.65)
for i,n in enumerate(tot):ax.text(n+12,i,f'{n:,}  ({n/20:.1f}%)',va='center',fontsize=10)
ax.set_yticks(y,names);ax.invert_yaxis();ax.set_xlim(0,max(tot)*1.29);ax.set_xlabel('唯一 oligo 数');ax.grid(axis='x',alpha=.12);ax.set_axisbelow(True)
ax=fig.add_subplot(gs[2,:]);ax.set_title('D  主区位置密度（806 个变异中 728 个位于主区）',loc='left',fontweight='bold',pad=16)
bins=np.arange(41703000,41751001,500)
for ids,label,col in [([r for r in windows if r[3].split('|')[0] in assigned and assigned[r[3].split('|')[0]]==4],'oligo 窗口中心',colors[3]),(read('claude_v7_variants'),'变异位置',colors[0])]:
 pos=[(int(r[1])+int(r[2]))/2 for r in ids];h,_=np.histogram(pos,bins)
 ax.stairs(h,bins/1e6,label=label,lw=1.5,color=col)
ax.axvspan(41.724175,41.730130,color=colors[0],alpha=.08,label='GPR82 基因跨度');ax.axvline(41.724175,color=grey,ls=':',lw=1)
ax.set_xlim(41.703,41.751);ax.set_ylabel('数量 / 500 bp');ax.set_xlabel('chrX 坐标（Mb）');ax.legend(frameon=False,fontsize=9,loc='upper right');ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x:.3f}'));ax.grid(axis='y',alpha=.12)
fig.subplots_adjust(left=.09,right=.96,top=.875,bottom=.10)
fig.text(.07,.027,'条数 ≠ 变异数：v7 共 806 个变异、498 条扫描、212 条对照；REF 可被多个比较共享。\n区域按窗口中心唯一归属；对照单列（116 条有基因组坐标，96 条无坐标）。本图不表示 MPRA 实测活性。',fontsize=10,color=grey)
for ext in ['png','svg']:fig.savefig(O/'figures'/f'05_v7_distribution.{ext}',dpi=180,bbox_inches='tight')
plt.close(fig)
# Dedicated main-locus annotation view.
fig,ax=plt.subplots(figsize=(16,9));fig.suptitle('Claude v7｜主区设计与功能注释',x=.05,y=.98,ha='left',fontsize=23,color=ink)
fig.text(.05,.908,'hg38 / chrX  ·  v7 实际 170 bp 窗口  ·  全库：806 个变异、498 条扫描、38 个手工位点、26 个古人类 ALT 位点',fontsize=11,color=grey)
def bars(rows,y,color,h=.35):
 intervals=[(int(r[1]),int(r[2])) for r in rows if r[0]=='chrX'];ax.broken_barh([(a/1e6,(b-a)/1e6) for a,b in intervals],(y-h/2,h),facecolors=color)
for r in read('evidence_regions'):
 if int(r[1])<41703000:continue
 bars([r],12,colors[0]);num=int(r[3].split('|')[0][1:]);ax.text((int(r[1])+int(r[2]))/2e6,12.4 if num not in [6,12] else 12.8,f'R{num:02}',ha='center',fontsize=9)
r=next(x for x in read('gpr82_transcripts') if x[3].startswith('ENST00000302548'));a=int(r[1]);ax.plot([a/1e6,int(r[2])/1e6],[11,11],color=colors[0],lw=1)
for offset,length in zip(r[11].strip(',').split(','),r[10].strip(',').split(',')):
 ax.broken_barh([((a+int(offset))/1e6,int(length)/1e6)],(10.83,.34),facecolors=colors[0])
ax.text(41.7308,11,'GPR82-201 →',va='center',color=colors[0])
layers=[('gpr82_promoter',10,colors[4]),('re2g',9,colors[3]),('ccre',8,colors[5]),('h3_macrophage_H3K27ac',6,colors[4]),('h3_ENCFF560PQQ',5,'#DEC381'),('claude_v7_tiles',3,colors[0])]
for key,y,c in layers:bars(read(key),y,c)
for key in ['dnase_ENCFF580ICE','dnase_ENCFF174FFX','dnase_ENCFF619TQP','dnase_ENCFF586HMB','dnase_ENCFF948HMV']:bars(read(key),7,'#368963')
for key in ['stat3_IFNg','stat3_IFNg-LPS','stat3_IL10','stat3_LPS','stat3_resting']:bars(read(key),4,'#986195')
for key,y,c in [('claude_v7_variants',2,colors[3]),('manual_38',1,colors[2]),('archaic_26',0,'#8156AA'),('sheet3_6',-1,colors[4])]:
 rows=read(key);ax.scatter([int(r[1])/1e6 for r in rows],[y]*len(rows),marker='|',s=36,color=c)
ax.set_yticks(list(range(12,-2,-1)),['证据区域','GPR82 外显子','设计启动子','rE2G 预测连接','ENCODE cCRE','巨噬 DNase 合并峰','巨噬 H3K27ac','单核 H3K27ac 代理','STAT3 合并峰','v7 扫描窗口','v7 变异位置','手工重点位点','古人类 ALT 位点','Sheet3 重点'])
ax.set_ylim(-1.7,13.4);ax.set_xlim(41.703,41.751);ax.axvline(41.724175,color=grey,lw=.8,ls=':');ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x:.3f}'));ax.set_xlabel('chrX 坐标（Mb）');ax.tick_params(axis='y',length=0);ax.spines['left'].set_visible(False);ax.grid(axis='x',alpha=.13);ax.set_axisbelow(True)
fig.subplots_adjust(left=.17,right=.97,top=.845,bottom=.15)
fig.text(.05,.035,'DNase / H3K27ac / STAT3 显示缓存峰区间，不是信号强度；合并多个条件不代表独立重复。\n古人类轨道为筛选出的 ALT 位点注释，不是古 DNA 原始 reads。rE2G 是预测连接；主区图未显示四个远端区域。',fontsize=10,color=grey)
for ext in ['png','svg']:fig.savefig(O/'figures'/f'06_v7_annotation.{ext}',dpi=180,bbox_inches='tight')
print(json.dumps(metadata,ensure_ascii=False,indent=2))
