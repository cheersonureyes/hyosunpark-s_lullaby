"""1층 1경간 평면골조: 단계별 직접강성법 / 가상일 분석.

실행: python frame_analysis.py
필수: numpy. 그림 저장에만 matplotlib 사용.
단위: N, m, Pa, rad.
"""
from dataclasses import dataclass, replace
from pathlib import Path
import csv
import json
import numpy as np


# STEP 1 ─ 입력: 이 부분을 수정해서 골조를 바꿉니다.
@dataclass(frozen=True)
class Parameters:
    span: float = 6.0
    height: float = 3.0
    E: float = 200e9
    nu: float = 0.3
    column_shear_factor: float = 5/6
    beam_shear_factor: float = 5/6
    column_A: float = 0.02
    column_I: float = 8e-5
    beam_A: float = 0.02
    beam_I: float = 8e-5
    delta: float = 0.01
    left_I_factor: float = 1.0
    right_I_factor: float = 1.0


def make_model(p):
    positive = [p.span, p.height, p.E, p.column_A, p.column_I,
                p.beam_A, p.beam_I, p.left_I_factor, p.right_I_factor,
                p.column_shear_factor, p.beam_shear_factor]
    if not np.isfinite(p.nu) or not -1 < p.nu < 0.5:
        raise ValueError('등방성 재료의 포아송비는 -1 < nu < 0.5여야 합니다.')
    if not np.all(np.isfinite(positive)) or min(positive) <= 0:
        raise ValueError('길이, 탄성계수, 단면값, 강성 배율은 양수여야 합니다.')
    if not np.isfinite(p.delta) or p.delta == 0:
        raise ValueError('참여율 계산에는 0이 아닌 유한한 변위가 필요합니다.')
    # A=0, B=1: 기초 / C=2, D=3: 보 양 끝
    nodes = np.array([[0., 0.], [p.span, 0.], [0., p.height], [p.span, p.height]])
    members = [
        ('C1', 0, 2, p.column_A, p.column_I*p.left_I_factor),
        ('G1', 2, 3, p.beam_A, p.beam_I),
        ('C2', 1, 3, p.column_A, p.column_I*p.right_I_factor),
    ]
    return nodes, members


# STEP 2 ─ 부재 강성행렬과 좌표변환, 전체 강성행렬 조립.
def local_stiffness(E, A, I, L, GA=np.inf):
    # 등단면 Timoshenko 요소의 닫힌형 강성. GA→∞이면 Euler–Bernoulli.
    phi = 12*E*I/(GA*L**2)
    a = E*A/L
    b, c = 12*E*I/(L**3*(1+phi)), 6*E*I/(L**2*(1+phi))
    d, e = (4+phi)*E*I/(L*(1+phi)), (2-phi)*E*I/(L*(1+phi))
    return np.array([
        [a,0,0,-a,0,0], [0,b,c,0,-b,c], [0,c,d,0,-c,e],
        [-a,0,0,a,0,0], [0,-b,-c,0,b,-c], [0,c,e,0,-c,d],
    ], dtype=float)


def assemble(p, nodes, members):
    K = np.zeros((3*len(nodes), 3*len(nodes)))
    elements = []
    for name, i, j, A, I in members:
        dx, dy = nodes[j]-nodes[i]
        L = np.hypot(dx, dy)
        c, s = dx/L, dy/L
        R = np.array([[c,s,0], [-s,c,0], [0,0,1]])
        T = np.zeros((6,6))
        T[:3,:3] = T[3:,3:] = R
        G = p.E/(2*(1+p.nu))
        As = A*(p.beam_shear_factor if name.startswith('G') else p.column_shear_factor)
        k = local_stiffness(p.E, A, I, L, G*As)
        dofs = np.array([3*i,3*i+1,3*i+2,3*j,3*j+1,3*j+2])
        K[np.ix_(dofs,dofs)] += T.T @ k @ T
        elements.append(dict(name=name, i=i, j=j, A=A, I=I, E=p.E,
                             G=G, As=As, GA=G*As, L=L, T=T, k=k, dofs=dofs))
    return K, elements


# STEP 3 ─ 지정 변위를 포함한 연립방정식. 구속 반력은 K u - f.
def solve_with_prescribed(K, loads, prescribed):
    fixed = np.array(sorted(prescribed), dtype=int)
    free = np.setdiff1d(np.arange(len(K)), fixed)
    u = np.zeros(len(K))
    u[fixed] = [prescribed[int(i)] for i in fixed]
    u[free] = np.linalg.solve(K[np.ix_(free,free)],
                            loads[free] - K[np.ix_(free,fixed)] @ u[fixed])
    return u, K @ u - loads


# STEP 4 ─ 실제 반력 패턴을 합력 1 N으로 정규화한 가상 해석.
def solve_states(K, delta):
    base = {i: 0.0 for i in range(6)}
    u, reactions = solve_with_prescribed(K, np.zeros(12), {**base, 6:delta, 9:delta})
    Q = reactions[6] + reactions[9]
    virtual_load = np.zeros(12)
    virtual_load[[6,9]] = reactions[[6,9]]/Q  # 수치상 합력 1 N
    z, _ = solve_with_prescribed(K, virtual_load, base)
    return u, reactions, Q, virtual_load, z


# STEP 5 ─ 단부력 평형으로 N,V,M 복원. 부재 중간 하중이 없는 경우.
def section_forces(el, local_u, x):
    q = el['k'] @ local_u
    return np.array([-q[0], -q[1], -q[2]+q[1]*x])


def local_displacement(el, ue, x):
    q = el['k'] @ ue
    ux = ue[0]+(ue[3]-ue[0])*x/el['L']
    # theta'=M/EI, v'=theta+V/GA를 적분. 전단변형도 그림에 반영.
    vy = ue[1]+ue[2]*x+(-q[2]*x*x/2+q[1]*x**3/6)/(el['E']*el['I'])-q[1]*x/el['GA']
    return np.array([ux,vy])


def member_contributions(elements, u, z, delta):
    rows = []
    # 등단면 부재에서 모멘트 곱은 2차식: 2점 Gauss 적분으로 정확.
    points, weights = np.polynomial.legendre.leggauss(2)
    for el in elements:
        ue = el['T'] @ u[el['dofs']]
        ze = el['T'] @ z[el['dofs']]
        contribution = np.zeros(3)
        energies = np.zeros(3)
        rigidity = np.array([el['E']*el['A'], el['GA'], el['E']*el['I']])
        for point, weight in zip(points, weights):
            x = (point+1)*el['L']/2
            actual = section_forces(el, ue, x)
            virtual = section_forces(el, ze, x)
            w = weight*el['L']/2
            contribution += actual*virtual/rigidity*w
            energies += 0.5*actual*actual/rigidity*w
        energy = float(0.5*ue @ el['k'] @ ue)
        np.testing.assert_allclose(sum(energies), energy, rtol=1e-9, atol=1e-9)
        np.testing.assert_allclose(local_displacement(el, ue, el['L']), ue[3:5], atol=1e-12)
        cn, cv, cm = contribution
        un, uv, um = energies
        rows.append(dict(member=el['name'], axial_contribution_m=float(cn),
                         shear_contribution_m=float(cv), bending_contribution_m=float(cm),
                         total_contribution_m=float(cn+cv+cm),
                         axial_percent=float(100*cn/delta), shear_percent=float(100*cv/delta),
                         bending_percent=float(100*cm/delta),
                         participation_percent=float(100*(cn+cv+cm)/delta),
                         axial_energy_J=float(un), shear_energy_J=float(uv),
                         bending_energy_J=float(um), energy_J=energy))
    return rows


# STEP 6 ─ 모든 단계를 연결하고 핵심 항등식 검증.
def analyze(p=Parameters()):
    nodes, members = make_model(p)
    K, elements = assemble(p, nodes, members)
    u, reactions, Q, virtual_load, z = solve_states(K, p.delta)
    rows = member_contributions(elements, u, z, p.delta)
    C = sum(r['total_contribution_m'] for r in rows)
    U = sum(r['energy_J'] for r in rows)
    np.testing.assert_allclose(C, p.delta, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(U, 0.5*Q*p.delta, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(z, u/Q, rtol=1e-9, atol=1e-12)
    for row in rows:
        np.testing.assert_allclose(row['participation_percent'], 100*row['energy_J']/U,
                                   rtol=1e-9, atol=1e-8)
    return dict(p=p, nodes=nodes, elements=elements, K=K, u=u, reactions=reactions,
                Q=Q, virtual_load=virtual_load, z=z, rows=rows,
                virtual_work_error_m=C-p.delta, energy_error_J=U-0.5*Q*p.delta)


# STEP 7 ─ 대칭성, 변위 배율, 비대칭성, 강성 민감도 검증.
def verify(p, result):
    rows = result['rows']
    # 독립 해석해: 끝단 횡력을 받는 Timoshenko 외팔보.
    E, A, I, L, GA, P = 200e9, .02, 8e-5, 3., 1e8, 1000.
    kc = local_stiffness(E,A,I,L,GA)
    uc, _ = solve_with_prescribed(kc, np.array([0.,0,0,0,P,0]), {0:0.,1:0.,2:0.})
    np.testing.assert_allclose(uc[4], P*L**3/(3*E*I)+P*L/GA, rtol=1e-10)
    np.testing.assert_allclose(uc[5], P*L**2/(2*E*I), rtol=1e-10)
    np.testing.assert_allclose(local_stiffness(E,A,I,L,1e30), local_stiffness(E,A,I,L), rtol=1e-12)
    symmetric = analyze(replace(p, left_I_factor=1.0, right_I_factor=1.0))
    np.testing.assert_allclose(symmetric['rows'][0]['energy_J'], symmetric['rows'][2]['energy_J'], rtol=1e-9)
    twice = analyze(replace(p, delta=2*p.delta))
    np.testing.assert_allclose(twice['Q'], 2*result['Q'], rtol=1e-9)
    for a,b in zip(rows, twice['rows']):
        np.testing.assert_allclose(b['energy_J'], 4*a['energy_J'], rtol=1e-9)
        np.testing.assert_allclose(b['participation_percent'], a['participation_percent'], rtol=1e-9)
    # 고정 횡하중에서 보 I 민감도를 중앙차분과 가상 해석으로 교차검증.
    base = {i:0.0 for i in range(6)}
    load = result['virtual_load']*result['Q']
    c = result['virtual_load']  # 기준 반력 패턴의 가중 평균 수평변위
    eps = 1e-4
    values = []
    for factor in [1-eps, 1+eps]:
        pp = replace(p, beam_I=p.beam_I*factor)
        nn, mm = make_model(pp)
        kk, _ = assemble(pp, nn, mm)
        uu, _ = solve_with_prescribed(kk, load, base)
        values.append(c@uu)
    fd = (values[1]-values[0])/(2*eps)
    predicted = -rows[1]['bending_contribution_m']
    np.testing.assert_allclose(fd, predicted, rtol=1e-6, atol=1e-12)
    asymmetric = analyze(replace(p, left_I_factor=0.5, right_I_factor=1.0))
    assert not np.isclose(asymmetric['rows'][0]['energy_J'], asymmetric['rows'][2]['energy_J'])
    return dict(status='PASS', virtual_work_error_m=result['virtual_work_error_m'],
                energy_error_J=result['energy_error_J'],
                beam_I_sensitivity_finite_difference_m=float(fd),
                beam_I_sensitivity_adjoint_m=float(predicted),
                checks=['Timoshenko cantilever exact solution', 'Euler-Bernoulli limit',
                        'virtual work', 'energy', 'energy integral vs matrix', 'deformed endpoints',
                        'symmetry', 'displacement scaling', 'asymmetry', 'stiffness sensitivity'])


# STEP 8 ─ CSV, JSON, 그림 저장.
def save_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def display_rows(result):
    by_name = {r['member']: r for r in result['rows']}
    rows = [by_name[n] for n in ['C1','C2','G1']]
    total = {k: sum(r[k] for r in rows) for k in rows[0] if k != 'member'}
    return rows+[dict(member='TOTAL', **total)]


def save_report(result, out):
    rows = display_rows(result)
    lines = ['# 축·전단·휨 변형 기여도', '',
             '전체 지정 수평변위에 대한 가상일 환산 기여량입니다. 비율의 분모는 전체 변위입니다.', '',
             '| 부재 | 축 (mm) | 전단 (mm) | 휨 (mm) | 합계 (mm) | 참여율 (%) |',
             '|---|---:|---:|---:|---:|---:|']
    for r in rows:
        values = [r[k+'_contribution_m']*1000 for k in ['axial','shear','bending','total']]
        lines.append('| '+r['member']+' | '+' | '.join(f'{v:.6f}' for v in values)
                     +f" | {r['participation_percent']:.4f} |")
    lines += ['', '## 전체 변위 대비 성분별 비율', '',
              '| 부재 | 축 (%) | 전단 (%) | 휨 (%) | 합계 (%) |', '|---|---:|---:|---:|---:|']
    for r in rows:
        lines.append('| '+r['member']+' | '+' | '.join(f"{r[k]:.4f}" for k in
                     ['axial_percent','shear_percent','bending_percent','participation_percent'])+' |')
    lines += ['', f"경간 {result['p'].span:g} m · 높이 {result['p'].height:g} m · "
              f"지정 변위 {result['p'].delta*1000:g} mm · 수평 반력 합 {result['Q']/1000:.4f} kN", '',
              'C1: 왼쪽 기둥, C2: 오른쪽 기둥, G1: 보.', '',
              '![골조와 성분별 기여도](participation.png)', '',
              '![보 휨강성 비교](beam_stiffness_sweep.png)']
    text = '\n'.join(lines)+'\n'
    (out/'results_summary.md').write_text(text, encoding='utf-8')
    print(text.split('##')[0])


def save_plot(result, sweep, path):
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.patches import Rectangle
    fonts = {f.name for f in font_manager.fontManager.ttflist}
    for candidate in ['Malgun Gothic', 'NanumGothic', 'Noto Sans CJK KR']:
        if candidate in fonts:
            plt.rcParams['font.family'] = candidate
            break
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['mathtext.fontset'] = 'dejavusans'
    plt.rcParams['font.size'] = 10
    p = result['p']
    rows = display_rows(result)
    fig = plt.figure(figsize=(13,10), layout='constrained', facecolor='#f5f7fb')
    gs = fig.add_gridspec(2,2, height_ratios=[2.2,1])
    ax = fig.add_subplot(gs[0,0])
    chart = fig.add_subplot(gs[0,1])
    table_ax = fig.add_subplot(gs[1,:])
    fig.suptitle('1층 1경간 라멘골조 | 축 · 전단 · 휨 기여도', fontsize=19, fontweight='bold')
    scale = 0.07*min(p.span,p.height)/abs(p.delta)
    for n, el in enumerate(result['elements']):
        start, end = result['nodes'][el['i']], result['nodes'][el['j']]
        ax.plot(*np.array([start,end]).T, color='#263c56', lw=5,
                label='원래 골조' if n==0 else None)
        ue = el['T']@result['u'][el['dofs']]
        x = np.linspace(0,el['L'],100)
        disp = local_displacement(el,ue,x)
        xy = start[:,None]+el['T'][:2,:2].T@(np.array([x,np.zeros_like(x)])+scale*disp)
        ax.plot(*xy, color='#e05c46', ls='--', lw=2.2,
                label=f'변형 형상 (×{scale:.1f})' if n==0 else None)
        mid = (start+end)/2
        offsets = {'C1':(-30,0), 'C2':(32,0), 'G1':(0,-24)}
        ax.annotate(el['name'], mid, xytext=offsets[el['name']], textcoords='offset points',
                    ha='center', va='center', fontsize=13, fontweight='bold',
                    bbox=dict(fc='white',ec='none',alpha=.85))
    for xx in [0,p.span]:
        ax.add_patch(Rectangle((xx-.035*p.span,-.06*p.height), .07*p.span,.06*p.height,
                              hatch='////', facecolor='#dbe2eb', edgecolor='#263c56'))
    yd = -.27*p.height
    ax.annotate('', (0,yd), (p.span,yd), arrowprops=dict(arrowstyle='<->',color='#607080'))
    ax.text(p.span/2,yd-.07*p.height,f'경간 L = {p.span:g} m',ha='center',va='top')
    xd = -.24*p.span
    ax.annotate('', (xd,0),(xd,p.height),arrowprops=dict(arrowstyle='<->',color='#607080'))
    ax.text(xd-.035*p.span,p.height/2,f'높이 H = {p.height:g} m',rotation=90,ha='right',va='center')
    for xx in [0,p.span]:
        ax.plot([xx,xx],[yd-.02*p.height,0],color='#bac4d0',lw=.8)
    for yy in [0,p.height]:
        ax.plot([xd,0],[yy,yy],color='#bac4d0',lw=.8)
    arrow_y = 1.16*p.height
    ax.annotate('',(.7*p.span+np.sign(p.delta)*.16*p.span,arrow_y),(.7*p.span,arrow_y),
                arrowprops=dict(arrowstyle='->',color='#e05c46',lw=2))
    ax.text(.5*p.span,1.3*p.height,f'보 양 끝 수평변위 Δ = {p.delta*1000:g} mm',
            ha='center',color='#b33d2b')
    ax.set(xlim=(-.4*p.span,1.25*p.span),ylim=(-.5*p.height,1.5*p.height))
    ax.set_aspect('equal')
    ax.axis('off')
    ax.legend(loc='lower center',bbox_to_anchor=(.5,-.04),frameon=False,ncol=2)
    modes = [('axial','축','#4878bf'),('shear','전단','#e9a442'),('bending','휨','#2b9b8f')]
    bottom = np.zeros(3)
    for key,label,color in modes:
        vals = np.array([r[key+'_percent'] for r in rows[:3]])
        chart.bar(['C1','C2','G1'], vals, bottom=bottom, width=.55, color=color,label=label)
        bottom += vals
    for i,val in enumerate(bottom):
        chart.text(i,val+.65,f'{val:.2f}%',ha='center',fontweight='bold')
    chart.set_ylim(0,max(bottom)*1.22)
    chart.set_ylabel('전체 수평변위 대비 기여율 (%)')
    chart.set_xticks(range(3), labels=[r['member']+'\n'+
                     '\n'.join(f'{label} {r[key+"_percent"]:.4f}%' for key,label in
                               [('axial','축'),('shear','전단'),('bending','휨')]) for r in rows[:3]], fontsize=9)
    chart.set_title('부재별 참여율과 변형 성분',pad=18)
    chart.spines[['top','right']].set_visible(False)
    chart.set_axisbelow(True)
    chart.grid(axis='y',alpha=.15)
    chart.legend(frameon=False,ncol=3,loc='upper right')
    table_ax.axis('off')
    columns = ['부재','축 (mm)','전단 (mm)','휨 (mm)','합계 (mm)','참여율 (%)']
    cells=[]
    for r in rows:
        cells.append([r['member']]+[f"{r[k+'_contribution_m']*1000:.6f}" for k in
                                   ['axial','shear','bending','total']]+[f"{r['participation_percent']:.4f}"])
    table=table_ax.table(cellText=cells,colLabels=columns,cellLoc='center',bbox=[0,.22,1,.72])
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    for (r,c),cell in table.get_celld().items():
        cell.set_edgecolor('white')
        if r==0:
            cell.set_facecolor('#263c56')
            cell.set_text_props(color='white',weight='bold')
        elif r==4:
            cell.set_facecolor('#dbeae9')
            cell.set_text_props(weight='bold')
        else:
            cell.set_facecolor('#edf1f6' if r%2 else '#ffffff')
    table_ax.text(0,.10,f'Timoshenko 요소  |  반력 합 Q = {result["Q"]/1000:.3f} kN  |  '
                  f'G = {p.E/(2*(1+p.nu))/1e9:.3f} GPa  |  '
                  f'As/A: 기둥 {p.column_shear_factor:.4f}, 보 {p.beam_shear_factor:.4f}',fontsize=10)
    table_ax.text(0,.01,'표의 mm는 가상일로 환산한 횡변위 기여량입니다. 부재 자체의 늘어남과는 다릅니다.',
                  fontsize=10,color='#596777')
    fig.savefig(path,dpi=200)
    fig.savefig(path.with_suffix('.svg'))
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8,4.7),layout='constrained')
    for name in ['C1','C2','G1']:
        data=[r for r in sweep if r['member']==name]
        ax.plot([r['beam_I_factor'] for r in data],[r['participation_percent'] for r in data],
                marker='o',ls='--' if name=='C2' else '-',label=name)
    ax.set(xscale='log',xlabel='보 I / 기본 보 I',ylabel='전체 변위 대비 참여율 (%)',
           title='보 휨강성 변화에 따른 참여율 (대칭 기둥 C1·C2는 겹침)')
    ax.set_xticks([0.1,1,10,100], labels=['0.1','1','10','100'])
    ax.legend()
    ax.grid(alpha=.2)
    fig.savefig(path.parent/'beam_stiffness_sweep.png',dpi=180)
    plt.close(fig)


def main():
    out = Path(__file__).resolve().parent/'results'
    out.mkdir(exist_ok=True)
    p = Parameters()
    result = analyze(p)
    checks = verify(p, result)
    save_csv(out/'baseline.csv', display_rows(result))
    save_report(result, out)
    sweep = []
    for factor in [0.1,0.25,0.5,1,2,4,10,100]:
        r = analyze(replace(p, beam_I=p.beam_I*factor))
        for row in r['rows']:
            sweep.append(dict(beam_I_factor=factor, horizontal_reaction_N=float(r['Q']), **row))
    save_csv(out/'beam_stiffness_sweep.csv', sweep)
    asym = analyze(replace(p, left_I_factor=0.5))
    save_csv(out/'asymmetric.csv', asym['rows'])
    save_csv(out/'nodal_results.csv', [dict(node=name, ux_m=float(result['u'][3*i]),
                 uy_m=float(result['u'][3*i+1]), rotation_rad=float(result['u'][3*i+2]),
                 reaction_x_N=float(result['reactions'][3*i]),
                 reaction_y_N=float(result['reactions'][3*i+1]),
                 reaction_moment_Nm=float(result['reactions'][3*i+2]))
                 for i,name in enumerate(['A','B','C','D'])])
    summary = dict(parameters=p.__dict__, total_horizontal_reaction_N=float(result['Q']),
                   lateral_stiffness_N_per_m=float(result['Q']/p.delta),
                   total_energy_J=sum(r['energy_J'] for r in result['rows']), verification=checks)
    (out/'verification.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding='utf-8')
    try:
        save_plot(result, sweep, out/'participation.png')
    except ImportError:
        print('matplotlib 미설치: 수치 결과만 저장했습니다.')
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    for row in result['rows']:
        print(f"{row['member']:14s}: {row['total_contribution_m']*1000:.6f} mm, "
              f"{row['participation_percent']:.6f}%")


if __name__ == '__main__':
    main()
