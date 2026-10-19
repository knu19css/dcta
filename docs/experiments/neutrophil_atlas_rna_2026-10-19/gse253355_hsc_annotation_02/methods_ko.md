# 방법 — 군집 잠정 세포 유형 후보 (R75)

1. 마커: node_markers_sign_review_2026-10-18.csv 원본 사본 (marker_input/). use_in_current_panel = true, is_annotation_candidate = true, marker_sign = ±1 행만 쓴다.
2. 모델: bm_dag_broad 와 bm_dag_slide2 를 따로 쓰고, 각각에 t_cell_broad_categories 를 더한다. node · gene · sign 중복은 지우고, 같은 node 의 양 · 음 충돌은 보류한다. DAG edge 는 쓰지 않는다.
3. 발현: log1p(1e4 · count / 모든 유전자 library). 유전자 전체 세포 평균으로 24 개 같은 개수 bin 을 만든다(동률은 symbol 순서).
4. 대조: 마커 유전자마다 같은 bin 의 비마커 유전자 50 개를 뽑는다(seed 20261020, Mersenne-Twister / Rejection, 마커 symbol 순서, 모자라면 전부). 목록은 control_genes.csv.
5. 유전자 보정값 = log 발현 − 대조 평균. 양성 · 음성 패널은 측정 유전자 동일 가중 평균이고, 없는 마커는 missing_markers 에 적는다.
6. 군집 점수 = 군집 평균 양성 − max(0, 군집 평균 음성). 음성 패널이 없으면 0 을 빼고, 양성 측정 유전자가 없으면 NA.
   선형인 부분은 양성 · 음성 패널 평균의 군집 집계까지다: 그래서 유전자별 군집 평균으로 계산했고 세포 단위 재계산과 같다. 음성 벌점 max(0, ·) 는 비선형이며 세포마다가 아니라 군집 평균 음성에 적용한다.
7. donor 동일 가중: donor x 군집 >= 5 세포 donor 의 양성 · 음성 평균을 같은 무게로 평균한 뒤 같은 식을 쓴다. 해당 donor 가 없으면 미판정이다.
8. 후보: 1 위가 0 이하이면 미판정이다. 2 위 = 다음 유형과 미판정 (0) 중 큰 것이다. 동률은 이름 하나로 깨지 않는다.
9. 합의: 두 모델 1 위가 같은 단일 유형이고 점수 > 0 → 잠정 후보. 둘 다 미판정 → 미판정. 그 밖 → 불일치.
10. 성숙 단계: 기존 6 단계 점수는 그대로다. Neutrophil 잠정 군집에만 기존 최고 단계 후보를 보조로 적고 unvalidated 로 둔다.
