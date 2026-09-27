"""다층 골조의 표·그림·층별 집계. multi_frame.py가 호출합니다."""
import numpy as np
from frame_analysis import save_csv, local_displacement


def number(value,digits=6):
    # 반올림하면 0인 미소 음수는 표시에서만 -0을 제거. 원시 CSV는 유지.
    return f'{0.0 if abs(value)<.5*10**(-digits) else value:.{digits}f}'


def aggregate(rows,label):
    keys=[k for k in rows[0] if k.endswith(('_m','_percent','_J'))]
    return dict(member=label,**{k:sum(r[k] for r in rows) for k in keys})


def table_text(rows):
    lines=['| 부재 | 축 (mm) | 전단 (mm) | 휨 (mm) | 합계 (mm) | 참여율 (%) |',
           '|---|---:|---:|---:|---:|---:|']
    for r in rows:
        vals=[r[k+'_contribution_m']*1000 for k in ['axial','shear','bending','total']]
        lines.append('| '+r['member']+' | '+' | '.join(number(v) for v in vals)
                     +f' | {r["participation_percent"]:.4f} |')
    return lines


def save_outputs(result,out):
    p=result['p']
    rows=sorted(result['rows'],key=lambda r:(r['floor'],r['kind']!='column',r['position']))
    total=aggregate(rows,'TOTAL')
    numeric=[{k:r[k] for k in total} for r in rows]+[total]
    save_csv(out/'contributions.csv',numeric)
    portions=[]
    for r in rows:
        item=dict(member=r['member'])
        for key in ['axial','shear','bending']:
            item[key+'_global_percent']=r[key+'_percent']
            item[key+'_within_member_percent']=(100*r[key+'_contribution_m']/r['total_contribution_m']
                                                if abs(r['total_contribution_m'])>1e-20 else float('nan'))
        portions.append(item)
    save_csv(out/'component_portions.csv',portions)
    save_csv(out/'member_map.csv',[{k:r[k] for k in ['member','kind','floor','position','start_node','end_node']} for r in rows])
    floors=[aggregate([r for r in rows if r['floor']==f],f'{f}F') for f in range(1,p.stories+1)]
    save_csv(out/'floor_contributions.csv',floors+[total])
    drift=[]
    nline=p.bays+1
    for floor in range(1,p.stories+1):
        for line in range(nline):
            upper=floor*nline+line
            lower=upper-nline
            d=float(result['u'][3*upper]-result['u'][3*lower])
            h=result['nodes'][upper,1]-result['nodes'][lower,1]
            drift.append(dict(floor=floor,column_line=line+1,drift_m=d,drift_ratio=d/h))
    save_csv(out/'story_drifts.csv',drift)
    save_csv(out/'nodal_results.csv',[
        dict(node=i,floor=i//nline,column_line=i%nline+1,x_m=float(x),y_m=float(y),
             ux_m=float(result['u'][3*i]),uy_m=float(result['u'][3*i+1]),rotation_rad=float(result['u'][3*i+2]),
             reaction_x_N=float(result['reactions'][3*i]),reaction_y_N=float(result['reactions'][3*i+1]),
             reaction_moment_Nm=float(result['reactions'][3*i+2]))
        for i,(x,y) in enumerate(result['nodes'])])
    maximum=max(r['participation_percent'] for r in rows)
    winners=[r['member'] for r in rows if np.isclose(r['participation_percent'],maximum,rtol=1e-7,atol=1e-8)]
    lines=[f'# {p.stories}층 {p.bays}경간 결과','',
           f'절점 {len(result["nodes"])}개 · 부재 {len(rows)}개 · 자유도 {len(result["u"])}개', '',
           f'최상층 모든 절점 ux={p.delta*1000:g} mm. 중간층 수평변위 및 모든 상부 절점의 uy·회전은 해석으로 결정합니다.', '',
           f'수평 구동력 합: **{result["Q"]/1000:.6f} kN**. 최대 참여 부재: **{", ".join(winners)}**, 각각 **{maximum:.4f}%**.', '',
           '## 부재별 축·전단·휨 기여도','']+table_text(numeric)
    lines+=['','## 층별 기여도','',
            '해당 층 기둥과 해당 층 상단 보를 묶은 집계입니다. 층간변위 자체의 기여도는 아닙니다.','']+table_text(floors+[total])
    lines+=['','## 성분별 portion (%)','',
            '전체 기준: 세 성분을 더하면 해당 부재 참여율입니다. 부재 내부 기준: 세 성분을 더하면 100%입니다.', '',
            '| 부재 | 전체 기준 축 | 전체 기준 전단 | 전체 기준 휨 | 부재 내부 축 | 부재 내부 전단 | 부재 내부 휨 |',
            '|---|---:|---:|---:|---:|---:|---:|']
    for r in portions:
        keys=[k+'_'+suffix for suffix in ['global_percent','within_member_percent'] for k in ['axial','shear','bending']]
        lines.append('| '+r['member']+' | '+' | '.join(number(r[k],4) for k in keys)+' |')
    lines+=['','비율은 최상층 지정 변위 기준입니다. 각 C/G 번호의 위치는 그림과 member_map.csv를 참조하세요.', '',
            '![골조와 부재별 기여도](overview.png)','',
            '![상세 결과 표](table_01.png)','',
            '부재가 24개보다 많으면 상세 표 그림을 24개씩 나누어 저장합니다.', '',
            '전단면적 As/A 기본값은 5/6인 예제 가정입니다. 해석 가정과 수정 방법은 MULTI_FRAME_GUIDE.md를 참조하세요.']
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    try:
        plot_result(result,numeric,out)
    except ImportError:
        print('matplotlib 미설치: CSV 및 Markdown 결과만 저장했습니다.')


def plot_result(result,numeric,out):
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.patches import Rectangle
    fonts={f.name for f in font_manager.fontManager.ttflist}
    for name in ['Malgun Gothic','NanumGothic','Noto Sans CJK KR']:
        if name in fonts:
            plt.rcParams['font.family']=name
            break
    plt.rcParams['axes.unicode_minus']=False
    p=result['p']
    nodes=result['nodes']
    xs,ys=np.unique(nodes[:,0]),np.unique(nodes[:,1])
    W,H=xs[-1],ys[-1]
    bw,hh=min(np.diff(xs)),min(np.diff(ys))
    samples=[]
    for el in result['elements']:
        x=np.linspace(0,el['L'],65)
        ue=el['T']@result['u'][el['dofs']]
        samples.append((el,x,local_displacement(el,ue,x)))
    max_disp=max(float(np.max(np.linalg.norm(d,axis=0))) for _,_,d in samples)
    scale=.12*min(bw,hh)/max_disp
    ranked=sorted(result['rows'],key=lambda r:r['participation_percent'],reverse=True)
    # 큰 모델에서도 읽을 수 있도록 개요 그래프는 상위 24개, 모든 수치는 표에 제공.
    selected=ranked[:24]
    fig,(ax,bar,portion_ax)=plt.subplots(1,3,figsize=(19,max(7,len(selected)*.28+2)),
                            gridspec_kw={'width_ratios':[1.3,1,.75]},layout='constrained')
    fig.suptitle(f'{p.stories}층 {p.bays}경간 라멘골조 | 축 · 전단 · 휨 기여도',fontsize=18,fontweight='bold')
    for index,(el,x,d) in enumerate(samples):
        start,end=nodes[el['i']],nodes[el['j']]
        ax.plot(*np.array([start,end]).T,color='#263c56',lw=3,label='원래 골조' if index==0 else None)
        xy=start[:,None]+el['T'][:2,:2].T@(np.array([x,np.zeros_like(x)])+scale*d)
        ax.plot(*xy,color='#df624a',lw=1.6,ls='--',label=f'변형 ×{scale:.1f}' if index==0 else None)
        is_column=el['name'].startswith('C')
        ax.annotate(el['name'],(start+end)/2,xytext=(-15,0) if is_column else (0,-13),
                    textcoords='offset points',ha='center',va='center',fontsize=9,fontweight='bold',
                    bbox=dict(fc='white',ec='none',alpha=.85,pad=1))
    for x in xs:
        ax.add_patch(Rectangle((x-.05*bw,-.08*hh),.1*bw,.08*hh,hatch='////',fc='#dce4ed',ec='#263c56'))
        ax.plot([x,x],[-.37*hh,0],lw=.7,color='#aebaca')
    for a,b in zip(xs[:-1],xs[1:]):
        ax.annotate('',(a,-.3*hh),(b,-.3*hh),arrowprops=dict(arrowstyle='<->',color='#667788'))
        ax.text((a+b)/2,-.39*hh,f'{b-a:g} m',ha='center',va='top',fontsize=9)
    for y in ys:
        ax.plot([-.3*bw,0],[y,y],lw=.7,color='#aebaca')
    for a,b in zip(ys[:-1],ys[1:]):
        ax.annotate('',(-.25*bw,a),(-.25*bw,b),arrowprops=dict(arrowstyle='<->',color='#667788'))
        ax.text(-.31*bw,(a+b)/2,f'{b-a:g} m',rotation=90,ha='right',va='center',fontsize=9)
    ax.text(W/2,-.85*hh,f'총 경간 {W:g} m  |  총 높이 {H:g} m',ha='center',fontsize=10)
    ax.text(W/2,H+.55*hh,f'최상층 ux = {p.delta*1000:g} mm',ha='center',color='#b33d2b')
    for x in xs:
        ax.annotate('',(x+np.sign(p.delta)*.12*bw,H+.18*hh),(x,H+.18*hh),
                    arrowprops=dict(arrowstyle='->',color='#df624a',lw=1.5))
    ax.set(xlim=(-.65*bw,W+.45*bw),ylim=(-1.1*hh,H+.9*hh))
    ax.set_aspect('equal')
    ax.axis('off')
    ax.legend(loc='lower center',bbox_to_anchor=(.5,-.04),ncol=2,frameon=False)
    left=np.zeros(len(selected))
    for key,label,color in [('axial','축','#4878bf'),('shear','전단','#e9a442'),('bending','휨','#2b9b8f')]:
        vals=np.array([r[key+'_percent'] for r in selected])
        bar.barh(np.arange(len(selected)),vals,left=left,color=color,label=label,height=.68)
        left+=vals
    for i,v in enumerate(left):
        bar.text(v+max(left)*.015,i,f'{v:.2f}%',va='center',fontsize=9)
    bar.set_yticks(np.arange(len(selected)),labels=[r['member'] for r in selected])
    bar.invert_yaxis()
    bar.set_xlim(0,max(left)*1.2)
    bar.set_xlabel('최상층 지정 변위 대비 기여율 (%)')
    bar.set_title('부재별 참여율 순위'+(' (상위 24개)' if len(ranked)>24 else ''),pad=32)
    bar.legend(loc='lower right',bbox_to_anchor=(1,1),frameon=False,ncol=3)
    bar.spines[['top','right']].set_visible(False)
    bar.set_axisbelow(True)
    bar.grid(axis='x',alpha=.15)
    # 작은 성분도 읽을 수 있도록 각 막대와 같은 높이에 숫자를 표시.
    portion_ax.set_ylim(bar.get_ylim())
    portion_ax.set_xlim(0,4)
    portion_ax.axis('off')
    portion_ax.set_title('성분별 portion (%)\n전체 최상층 변위 기준',pad=25,fontsize=11)
    headers=[('축','#4878bf'),('전단','#b67a20'),('휨','#2b9b8f')]
    for col,(label,color) in enumerate(headers):
        portion_ax.text(.65+col*1.1,1.015,label,transform=portion_ax.get_xaxis_transform(),
                        ha='center',color=color,fontweight='bold',fontsize=10)
    for i,r in enumerate(selected):
        if i%2==0:
            portion_ax.axhspan(i-.42,i+.42,color='#f0f3f7',zorder=0)
        for col,key in enumerate(['axial','shear','bending']):
            portion_ax.text(.65+col*1.1,i,number(r[key+'_percent'],4),ha='center',va='center',fontsize=9)
    fig.savefig(out/'overview.png',dpi=180)
    fig.savefig(out/'overview.svg')
    plt.close(fig)

    # 표는 24부재 단위로 나누고 마지막 페이지에 전체 합계를 표시.
    data=numeric[:-1]
    for offset in range(0,len(data),24):
        chunk=data[offset:offset+24]
        if offset+24>=len(data):
            chunk=chunk+[numeric[-1]]
        fig,ax=plt.subplots(figsize=(12,2+.38*(len(chunk)+1)),layout='constrained')
        ax.axis('off')
        ax.set_title(f'{p.stories}층 {p.bays}경간 | 부재별 성분 기여량',fontsize=16,pad=20)
        cells=[[r['member']]+[number(r[k+'_contribution_m']*1000) for k in ['axial','shear','bending','total']]
               +[f"{r['participation_percent']:.4f}"] for r in chunk]
        table=ax.table(cellText=cells,colLabels=['부재','축 (mm)','전단 (mm)','휨 (mm)','합계 (mm)','참여율 (%)'],
                       cellLoc='center',bbox=[0,.13,1,.85])
        table.auto_set_font_size(False)
        table.set_fontsize(11)
        for (r,c),cell in table.get_celld().items():
            cell.set_edgecolor('white')
            if r==0:
                cell.set_facecolor('#263c56'); cell.set_text_props(color='white',weight='bold')
            elif chunk[r-1]['member']=='TOTAL':
                cell.set_facecolor('#dbeae9'); cell.set_text_props(weight='bold')
            else:
                cell.set_facecolor('#edf1f6' if r%2 else 'white')
        ax.text(0,.04,'mm는 가상일로 환산한 최상층 횡변위 기여량입니다. 비율의 분모는 전체 지정 변위입니다.',fontsize=10)
        name=f'table_{offset//24+1:02d}'
        fig.savefig(out/f'{name}.png',dpi=180)
        fig.savefig(out/f'{name}.svg')
        plt.close(fig)
