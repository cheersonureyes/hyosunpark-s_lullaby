"""실행: python test_multi_frame.py. 수치·회귀·대칭·분할·민감도 검증."""
from dataclasses import replace
from pathlib import Path
import json
import numpy as np
from multi_frame import Parameters,analyze,make_model,boundary_dofs
from frame_analysis import (analyze as single_analyze, Parameters as SingleParameters,
                            assemble,solve_with_prescribed,member_contributions,local_stiffness)


def check_case(p):
    r=analyze(p)
    assert len(r['nodes'])==(p.stories+1)*(p.bays+1)
    assert len(r['rows'])==p.stories*(2*p.bays+1)
    twice=analyze(replace(p,delta=2*p.delta))
    np.testing.assert_allclose(twice['Q'],2*r['Q'],rtol=1e-8)
    for a,b in zip(r['rows'],twice['rows']):
        np.testing.assert_allclose(a['participation_percent'],b['participation_percent'],rtol=1e-8,atol=1e-8)
        np.testing.assert_allclose(b['energy_J'],4*a['energy_J'],rtol=1e-8,atol=1e-8)
    # 대칭 모델의 각 층에서 거울 위치 부재끼리 비교.
    if not p.bay_widths and p.left_I_factor==p.right_I_factor:
        for row in r['rows']:
            n=p.bays+1 if row['kind']=='column' else p.bays
            mirror=next(x for x in r['rows'] if x['floor']==row['floor'] and
                        x['kind']==row['kind'] and x['position']==n+1-row['position'])
            np.testing.assert_allclose(row['participation_percent'],mirror['participation_percent'],atol=1e-7)
    # 모든 부재를 두 요소로 나누어 기존 절점 응답 및 물리적 부재의 기여도 보존 검증.
    nodes,members,_=make_model(p)
    refined_nodes=nodes.tolist()
    refined_members=[]
    for name,i,j,A,I in members:
        mid=len(refined_nodes)
        refined_nodes.append(((nodes[i]+nodes[j])/2).tolist())
        refined_members.extend([(name+'_a',i,mid,A,I),(name+'_b',mid,j,A,I)])
    K,els=assemble(p,np.array(refined_nodes),refined_members)
    base,roof=boundary_dofs(p)
    u,rx=solve_with_prescribed(K,np.zeros(len(K)),{**base,**{int(i):p.delta for i in roof}})
    np.testing.assert_allclose(u[:len(r['u'])],r['u'],rtol=1e-7,atol=1e-11)
    Q=sum(rx[roof])
    vf=np.zeros(len(K)); vf[roof]=rx[roof]/Q
    z,_=solve_with_prescribed(K,vf,base)
    rr=member_contributions(els,u,z,p.delta)
    for i,row in enumerate(r['rows']):
        for key in ['axial','shear','bending','total']:
            field=key+'_contribution_m'
            np.testing.assert_allclose(rr[2*i][field]+rr[2*i+1][field],row[field],rtol=1e-7,atol=1e-11)
    # 모든 보 I를 같은 비율로 변경. 고정 하중 응답의 민감도는 보 휨 기여량 합의 음수.
    eps=1e-4
    vals=[]
    for factor in [1-eps,1+eps]:
        pp=replace(p,beam_I=p.beam_I*factor)
        nn,mm,_=make_model(pp)
        kk,_=assemble(pp,nn,mm)
        uu,_=solve_with_prescribed(kk,r['virtual_load']*r['Q'],base)
        vals.append(r['virtual_load']@uu)
    fd=(vals[1]-vals[0])/(2*eps)
    adj=-sum(x['bending_contribution_m'] for x in r['rows'] if x['kind']=='beam')
    np.testing.assert_allclose(fd,adj,rtol=2e-5,atol=1e-10)
    return dict(stories=p.stories,bays=p.bays,node_count=len(nodes),member_count=len(members),
                Q_N=r['Q'],status='PASS',refinement_Q_error_N=float(Q-r['Q']),
                sensitivity_difference_m=float(fd-adj),**{k:v for k,v in r['verification'].items() if k!='status'})


def main():
    original=single_analyze(SingleParameters())
    generalized=analyze(Parameters(stories=1,bays=1))
    np.testing.assert_allclose(original['u'],generalized['u'],rtol=1e-9,atol=1e-12)
    old={r['member']:r for r in original['rows']}
    for r in generalized['rows']:
        np.testing.assert_allclose(r['participation_percent'],old[r['member']]['participation_percent'],rtol=1e-9)
    # 독립 외팔보 해석해와 요소 강성 확인.
    E,A,I,L,GA,P=200e9,.02,8e-5,3.,1e8,1000.
    K=local_stiffness(E,A,I,L,GA)
    u,_=solve_with_prescribed(K,np.array([0.,0,0,0,P,0]),{0:0.,1:0.,2:0.})
    np.testing.assert_allclose(u[4],P*L**3/(3*E*I)+P*L/GA,rtol=1e-10)
    np.testing.assert_allclose(u[5],P*L**2/(2*E*I),rtol=1e-10)
    cases=[Parameters(stories=s,bays=b) for s,b in [(1,1),(2,2),(3,3),(2,3),(3,2),(5,5)]]
    cases.append(Parameters(stories=3,bays=3,bay_widths=(4.,5.,6.),story_heights=(3.5,3.,2.8),
                            delta=-.008,left_I_factor=.7,beam_shear_factor=.55))
    results=[check_case(p) for p in cases]
    for bad in [Parameters(stories=0),Parameters(bays=2,bay_widths=(4.,)),Parameters(delta=0)]:
        try:
            analyze(bad)
        except ValueError:
            pass
        else:
            raise AssertionError('잘못된 입력이 거부되지 않았습니다.')
    payload=dict(status='PASS',checks=['single-frame regression','cantilever exact solution',
                 'member/node counts','virtual work','energy','force/moment equilibrium','symmetry',
                 'displacement scaling','two-element refinement','beam stiffness sensitivity','invalid inputs'],cases=results)
    out=Path(__file__).resolve().parent/'results'/'multi_verification.json'
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print(json.dumps(payload,indent=2))


if __name__=='__main__':
    main()
