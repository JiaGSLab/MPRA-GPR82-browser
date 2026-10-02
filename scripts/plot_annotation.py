"""Plot published BED tracks; no private project inputs required."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
O=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':['PingFang SC','Arial Unicode MS','DejaVu Sans'],'font.size':10,
 'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,'svg.fonttype':'path'})
def read(key):return [l.split('\t') for l in (O/'tracks'/f'{key}.bed').read_text().splitlines()]
def merged(keys):
    xs=sorted((int(r[1]),int(r[2])) for k in keys for r in read(k));out=[]
    for a,b in xs:
        if out and a<=out[-1][1]:out[-1][1]=max(out[-1][1],b)
        else:out.append([a,b])
    return out
def bars(ax,intervals,y,color,height=.36):
    ax.broken_barh([(a/1e6,(b-a)/1e6) for a,b in intervals],(y-height/2,height),facecolors=color)
blue='#286E9F';teal='#218579';orange='#C1702B';purple='#8563A6';grey='#9AA8B5';ink='#23394D'
fig=plt.figure(figsize=(16,13));gs=fig.add_gridspec(3,1,height_ratios=[1.2,3.6,1.35],hspace=.33)
fig.suptitle('GPR82｜区域注释与设计位置',x=.06,ha='left',fontsize=23,color=ink,y=.98)
fig.text(.06,.947,'hg38 / chrX  ·  Claude v6、Codex C1-2000 与 Claude v7  ·  峰区间及设计位置，不是实测 MPRA 活性',fontsize=11,color='#65788A')
ax=fig.add_subplot(gs[0]);ax.set_title('A  全景：四个远端证据区与主设计区',loc='left',fontsize=13,fontweight='bold',pad=12)
for i,r in enumerate(read('design_scope')):
    a,b=int(r[1]),int(r[2]);bars(ax,[(a,b)],2,blue if i<4 else teal,.35)
    lab=r[3].split('|')[0];txt=f'{lab}\n{[75,30,146,40,5349][i]} 条（v6）'
    ax.annotate(txt,((a+b)/2e6,2.2),((a+b)/2e6,2.62 if i!=3 else 3.03),ha='center',fontsize=9,
                arrowprops=dict(arrowstyle='-',color=grey,lw=.5))
for name,y,col in [('CASK',.55,grey),('GPR34',1.05,purple),('GPR82',1.05,teal),('NYX',1.05,grey)]:
    r=next(x for x in read('gene_spans') if x[3]==name);a,b=int(r[1]),int(r[2]);bars(ax,[(a,b)],y,col,.13)
    ax.text((max(a,41435000)+min(b,41755000))/2e6,y-.16,name+(' ←' if r[5]=='-' else ' →'),ha='center',va='top',fontsize=9,color=col)
ax.set_xlim(41.435,41.755);ax.set_ylim(0,3.6);ax.set_yticks([]);ax.set_xlabel('chrX 坐标（Mb）')
ax.axvspan(41.703,41.751,color=teal,alpha=.06)

ax=fig.add_subplot(gs[1]);ax.set_title('B  主区注释：证据峰、扫描及核心位点',loc='left',fontsize=13,fontweight='bold',pad=18)
for j,r in enumerate(read('evidence_regions')):
    a,b=int(r[1]),int(r[2]);
    if b<41703000:continue
    bars(ax,[(a,b)],11,blue)
    ax.text((a+b)/2e6,11.38 if j not in [5,11] else 11.72,r[3].split('|')[0],ha='center',fontsize=9)
r=next(x for x in read('gpr82_transcripts') if x[3].startswith('ENST00000302548'))
a,b=int(r[1]),int(r[2]);ax.plot([a/1e6,b/1e6],[10,10],color=teal,lw=.8)
starts=[int(x) for x in r[11].strip(',').split(',')];lens=[int(x) for x in r[10].strip(',').split(',')]
bars(ax,[(a+s,a+s+n) for s,n in zip(starts,lens)],10,teal,.32)
bars(ax,[(int(r[6]),int(r[7]))],10,teal,.55)
ax.text(41.7308,10,'GPR82-201 →',va='center',fontsize=9,color=teal)
layers=[(['gpr82_promoter'],9,orange),(['re2g'],8,purple),
    (['dnase_ENCFF580ICE','dnase_ENCFF174FFX','dnase_ENCFF619TQP','dnase_ENCFF586HMB','dnase_ENCFF948HMV'],7,'#2B8F63'),
    (['h3_macrophage_H3K27ac'],6,'#CC7333'),(['h3_ENCFF560PQQ'],5,'#DAA959'),
    (['stat3_IFNg','stat3_IFNg-LPS','stat3_IL10','stat3_LPS','stat3_resting'],4,'#9C578E'),
    (['evidence_tiles','survey_tiles'],3,teal),(['claude_v7_tiles'],2,orange)]
for keys,y,color in layers:bars(ax,merged(keys),y,color)
for key,y,color in [('manual_38',1,'#D62728'),('archaic_26',0,purple),('sheet3_6',-1,'#E68919')]:
    xx=[int(r[1])/1e6 for r in read(key)];ax.scatter(xx,[y]*len(xx),s=24,marker='|',color=color)
ax.set_yticks(list(range(11,-2,-1)),['13 个证据区','GPR82 外显子 / CDS','设计启动子','rE2G 预测连接区','DNase 合并峰','巨噬 H3K27ac 合并峰','单核 H3K27ac 代理','STAT3 合并峰','C1 扫描覆盖','v7 扫描覆盖','手工重点','古人类 ALT','Sheet3 重点'])
ax.set_ylim(-1.7,12);ax.set_xlim(41.703,41.751);ax.tick_params(axis='y',length=0);ax.spines['left'].set_visible(False)
ax.grid(axis='x',color='#E6ECF1',lw=.6);ax.set_axisbelow(True);ax.axvline(41.724175,color=ink,ls=':',lw=1)
ax.set_xlabel('chrX 坐标（Mb）');ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x:.3f}'))

ax=fig.add_subplot(gs[2]);ax.set_title('C  主区变异密度（每 500 bp；按变异计数，不按 oligo 条数）',loc='left',fontsize=13,fontweight='bold',pad=12)
bins=np.arange(41703000,41751001,500)
for key,label,col in [('variants_6000','Claude v6：1574 个（全范围）',grey),('variants_plan2000','C1-2000：691 个',teal),('claude_v7_variants','Claude v7：806 个',orange)]:
    xs=[int(r[1]) for r in read(key)];h,_=np.histogram(xs,bins=bins);ax.stairs(h,bins/1e6,color=col,label=label,lw=1.4)
ax.set_xlim(41.703,41.751);ax.set_ylabel('变异数 / 500 bp');ax.set_xlabel('chrX 坐标（Mb）');ax.legend(loc='upper right',frameon=False,fontsize=9)
ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x:.3f}'));ax.grid(axis='y',color='#E6ECF1',lw=.6)
fig.subplots_adjust(left=.19,right=.975,top=.9,bottom=.095)
fig.text(.06,.025,'C1-2000：200 bp insert，1998 条现有序列 + 2 个待设计名额；Claude v7：170 bp insert，2000 条现有序列。\n'
    'B 的 DNase / STAT3 为来源合并峰，不能据条带高度判断信号强弱。R09 局部启动子强证据不能外推到整个区间。',fontsize=10,color='#65788A')
for ext in ['png','svg']:fig.savefig(O/'figures'/('04_genome_annotation.'+ext),dpi=180,bbox_inches='tight')
print('Saved annotation PNG/SVG')
