# hyosunpark-s_lullaby
hyosunpark's highrise building class

## 라멘골조 가상일 기여도 분석

Python과 NumPy로 2D 라멘골조의 축·전단·휨 변형이 최상층 수평변위에 기여하는 양을 계산합니다. OpenSees 없이 Timoshenko 요소의 직접강성법과 단위 가상하중 해석을 수행합니다.

- [다층·다경간 실행 설명서](frame_virtual_work/MULTI_FRAME_GUIDE.md)
- [1층 1경간의 이론 및 단계별 설명](frame_virtual_work/README.md)
- [2층 2경간 결과](frame_virtual_work/results/2F_2B/report.md)
- [3층 3경간 결과](frame_virtual_work/results/3F_3B/report.md)

### 실행

Python 3.9 이상에서:

```bash
cd frame_virtual_work
python -m pip install -r requirements.txt
python multi_frame.py --examples
python test_multi_frame.py
```

직접 층수·경간수를 지정하려면:

```bash
python multi_frame.py --stories 3 --bays 3 --delta-mm 10
```

기본 조건은 경간마다 6 m, 층고마다 3 m, 최상층 모든 절점에 같은 수평변위 10 mm입니다. 중간층 변위는 해석으로 구합니다. 차트에 전체 최상층 변위 대비 축·전단·휨 portion %를 표시하며, 부재 내부 성분 비율은 보고서에 별도 제공합니다.

![3층 3경간 골조 및 기여도](frame_virtual_work/results/3F_3B/overview.png)

선형탄성·소변형 교육용 모델입니다. 단면과 유효전단면적 등 입력 가정은 설명서를 참조하세요.
