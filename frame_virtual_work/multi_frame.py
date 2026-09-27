"""다층·다경간 해석. 실행: python multi_frame.py --stories 3 --bays 3
예제 전체: python multi_frame.py --examples
기존 frame_analysis.py의 검증된 요소·적분 함수를 재사용합니다.
"""
from dataclasses import dataclass, replace, asdict
from pathlib import Path
import argparse
import json
import numpy as np
from frame_analysis import (Parameters as SingleParameters, assemble, solve_with_prescribed,
                            member_contributions, local_stiffness, save_csv)


# STEP 1: 층수·경간수와 필요 시 서로 다른 경간·층고 입력.
@dataclass(frozen=True)
class Parameters(SingleParameters):
    stories: int = 2
    bays: int = 2
    bay_widths: tuple = ()
    story_heights: tuple = ()


def geometry(p):
    for value in [p.stories, p.bays]:
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError('층수와 경간수는 1 이상의 정수여야 합니다.')
    widths = np.asarray(p.bay_widths if len(p.bay_widths) else [p.span]*p.bays, dtype=float)
    heights = np.asarray(p.story_heights if len(p.story_heights) else [p.height]*p.stories, dtype=float)
    if widths.shape != (p.bays,) or heights.shape != (p.stories,):
        raise ValueError('bay_widths와 story_heights의 길이가 각각 경간수·층수와 같아야 합니다.')
    positive = [*widths, *heights, p.E,p.column_A,p.column_I,p.beam_A,p.beam_I,
                p.column_shear_factor,p.beam_shear_factor,p.left_I_factor,p.right_I_factor]
    if not np.all(np.isfinite(positive)) or min(positive) <= 0:
        raise ValueError('치수·재료·단면·강성 배율은 유한한 양수여야 합니다.')
    if not np.isfinite(p.nu) or not -1 < p.nu < .5:
        raise ValueError('포아송비 범위: -1 < nu < 0.5')
    if not np.isfinite(p.delta) or p.delta == 0:
        raise ValueError('지정 변위는 0이 아닌 유한한 값이어야 합니다.')
    return np.r_[0,np.cumsum(widths)], np.r_[0,np.cumsum(heights)]


# STEP 2: 층별로 절점과 부재 생성. C/G 번호는 아래층부터 왼쪽→오른쪽.
def make_model(p):
    xs, ys = geometry(p)
    nline = p.bays+1
    nodes = np.array([(x,y) for y in ys for x in xs])
    members, metadata = [], []
    for floor in range(1,p.stories+1):
        for line in range(nline):
            name = f'C{(floor-1)*nline+line+1}'
            factor = p.left_I_factor if line==0 else p.right_I_factor if line==p.bays else 1.
            i, j = (floor-1)*nline+line, floor*nline+line
            members.append((name,i,j,p.column_A,p.column_I*factor))
            metadata.append(dict(member=name, kind='column', floor=floor, position=line+1,
                                 start_node=i, end_node=j))
        for bay in range(p.bays):
            name = f'G{(floor-1)*p.bays+bay+1}'
            i, j = floor*nline+bay, floor*nline+bay+1
            members.append((name,i,j,p.beam_A,p.beam_I))
            metadata.append(dict(member=name,kind='beam',floor=floor,position=bay+1,
                                 start_node=i,end_node=j))
    return nodes, members, metadata


def boundary_dofs(p):
    nline=p.bays+1
    base={i:0.0 for i in range(3*nline)}
    roof=np.arange(p.stories*nline,(p.stories+1)*nline)*3
    return base, roof


# STEP 3: 최상층의 모든 절점에 같은 ux=delta 지정. 중간층은 자유.
def analyze(p=Parameters()):
    nodes,members,metadata=make_model(p)
    K,elements=assemble(p,nodes,members)
    base,roof=boundary_dofs(p)
    prescribed={**base,**{int(i):p.delta for i in roof}}
    u,reactions=solve_with_prescribed(K,np.zeros(len(K)),prescribed)
    Q=float(sum(reactions[roof]))
    if Q*p.delta <= 0:
        raise ValueError('양의 횡강성을 얻지 못했습니다. 모델을 확인하세요.')

    # STEP 4: 구동 반력 패턴 / Q = 합력 1 N인 가상하중.
    virtual=np.zeros(len(K))
    virtual[roof]=reactions[roof]/Q
    z,_=solve_with_prescribed(K,virtual,base)

    # STEP 5: 각 부재의 축·전단·휨 적분. 기존 검증된 함수를 재사용.
    rows=member_contributions(elements,u,z,p.delta)
    for row,meta in zip(rows,metadata):
        row.update(meta)

    # STEP 6: 평형·가상일·에너지 검증. 반력 모멘트도 포함.
    C=sum(r['total_contribution_m'] for r in rows)
    U=sum(r['energy_J'] for r in rows)
    np.testing.assert_allclose(C,p.delta,rtol=1e-8,atol=1e-12)
    np.testing.assert_allclose(U,.5*Q*p.delta,rtol=1e-8,atol=1e-8)
    np.testing.assert_allclose(z,u/Q,rtol=1e-8,atol=1e-12)
    free=np.setdiff1d(np.arange(len(K)),list(prescribed))
    free_res=float(np.max(np.abs(reactions[free])))
    np.testing.assert_allclose(reactions[free],0,atol=max(abs(Q),1)*1e-8)
    force=reactions.reshape(-1,3)
    np.testing.assert_allclose(force[:,:2].sum(axis=0),0,atol=max(abs(Q),1)*1e-8)
    moment=float(np.sum(nodes[:,0]*force[:,1]-nodes[:,1]*force[:,0]+force[:,2]))
    np.testing.assert_allclose(moment,0,atol=max(abs(Q)*np.max(nodes),1)*1e-8)
    for r in rows:
        np.testing.assert_allclose(r['participation_percent'],100*r['energy_J']/U,atol=1e-7)
    checks=dict(status='PASS',virtual_work_error_m=float(C-p.delta),
                energy_error_J=float(U-.5*Q*p.delta),free_dof_residual_max=free_res,
                moment_equilibrium_error_Nm=moment)
    return dict(p=p,nodes=nodes,elements=elements,rows=rows,K=K,u=u,z=z,Q=Q,
                reactions=reactions,virtual_load=virtual,verification=checks)


# STEP 7: 층별·부재별 결과와 그림 저장.
def export(result,out):
    from multi_frame_reports import save_outputs
    out=Path(out)
    out.mkdir(parents=True,exist_ok=True)
    save_outputs(result,out)
    summary=dict(parameters=asdict(result['p']),node_count=len(result['nodes']),
                 member_count=len(result['rows']),dof_count=len(result['u']),
                 horizontal_reaction_N=result['Q'],lateral_stiffness_N_m=result['Q']/result['p'].delta,
                 verification=result['verification'])
    (out/'verification.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    return summary


# STEP 8: 명령행 실행. 기존 1층 1경간은 별도 스크립트 그대로 지원.
def main():
    parser=argparse.ArgumentParser(description='다층·다경간 Timoshenko 라멘골조')
    parser.add_argument('--stories',type=int,default=2)
    parser.add_argument('--bays',type=int,default=2)
    parser.add_argument('--span',type=float,default=6.)
    parser.add_argument('--height',type=float,default=3.)
    parser.add_argument('--delta-mm',type=float,default=10.)
    parser.add_argument('--examples',action='store_true')
    args=parser.parse_args()
    sizes=[(1,1),(2,2),(3,3)] if args.examples else [(args.stories,args.bays)]
    for stories,bays in sizes:
        p=Parameters(stories=stories,bays=bays,span=args.span,height=args.height,delta=args.delta_mm/1000)
        r=analyze(p)
        out=Path(__file__).resolve().parent/'results'/f'{stories}F_{bays}B'
        export(r,out)
        top=max(r['rows'],key=lambda x:x['participation_percent'])
        print(f'{stories}F {bays}B: {len(r["nodes"])} nodes, {len(r["rows"])} members; '
              f'Q={r["Q"]/1000:.4f} kN; max={top["member"]} {top["participation_percent"]:.4f}%; PASS')


if __name__=='__main__':
    main()
