# 단일 경로–분기 경로 비교 구현 계획

Goal: 같은 세포의 RNA count를 단일 곡선과 분기별 곡선으로 설명하고, 학습에 쓰지 않은 세포/공여자에서 예측을 비교하는 실행 가능한 R 시제품.
Architecture: 고정된 Route 후보(세포 ID, 경로상 위치, 공통 구간과 분기 배정)를 받아 유전자별 count 모형을 학습한다. 후보 생성, 생물학적 라벨, DAG 매칭과 분리된 비교 함수이며, 현재 단계에서 실제 분화를 입증하지 않는다.
Tech stack: R, stats, splines, MASS; testthat가 가용하면 테스트에 사용. 패키지 설치 없이 기존 dcta 환경을 우선한다.
Spec: 이 문서는 공개 가능한 구현 계약이다. 내부 조율 기록과 실자료 분석은 게시하지 않는다.

## 승인·소유권

- 사용자는 별도 branch의 구현 및 GitHub 게시를 명시 승인했다. 결과를 직접 검토하기 전 merge하지 않는다.
- Claude Desktop Route 세션이 아래 R·테스트·예제·결과 요약을 구현한다. Codex는 명세/검증/commit/push 담당이다. 다른 작업자가 있으므로 아래 범위 밖 수정/되돌림은 금지한다.
- 컴퓨터 유즈 없이 전달 전용 relay와 bridge 파일로 소통한다. 요청 모델은 Opus 5.5·고속 모드이며 실제 메타데이터와 제한은 보고한다. 이전 '검토만' 메시지의 쓰기 제한은 이번 구현 범위에서 해제된다.
- 원본 Route/Flow/DAG와 원자료를 수정하지 않는다. 운영 annotation으로 자동 채택하거나 merge하지 않는다.

## 파일과 인터페이스

- `R/branch_candidates/branch_curve_comparison.R`: 입력 검증, split 생성, 모형 학습·예측 비교. 외부 원자료 경로에 의존하지 않는 독립 함수.
- `R/branch_candidates/run_branch_curve_demo.R`: 재현 가능한 작은 모의자료 실행, TSV와 PNG 생성.
- `tests/R/branch_candidates/test_branch_curve_comparison.R`: 아래 중요한 동작 검증.
- `tests/R/branch_candidates/run_tests.R`: 저장소 root에서 실행 가능한 테스트 진입점.
- `docs/branch_candidates/README_ko.md`: 사용법·목적함수 설명·제한. 수식을 사용하면 모든 기호와 값의 해석을 설명한다.
- `docs/branch_candidates/examples/`: 작은 합성자료 결과/그림만. 실제 세포별 원자료는 게시하지 않는다.
- `results/branch_candidates_2026-10-16/implementation_report.md`: 실행 환경/명령/결과/한계/변경 파일; 로컬 bridge.

주 함수 `compare_branch_curves(counts, cell_data, split_point, spline_df=3L, folds=NULL, seed=1L, geometry_source="fixed_full_graph")`:
- counts: 유전자×세포의 비음수 정수 count, 고유 행/열 이름. log-normalized 값에 음이항 모형을 적용하지 않는다.
- cell_data: cell_id, donor_id, position(0~1), arm, library_size(양수). optional replicate_id는 실제 별도 채취/라이브러리 구분이며 가상 공여자가 아니다.
- arm: 공통 구간은 trunk, 분기 뒤는 2개 이상의 후보 arm. position/arm은 후보 탐색 산출물이며 시간/진짜 cell type이라는 주장은 하지 않는다.
- split_point: 사전 고정 후보의 분기 위치. 평가 자료의 예측 오차로 같은 단계에서 골라 최종 성능처럼 보고하지 않는다.
- 출력: 후보/유전자/fold별 동일 평가 세포 수, 예측 손실, 차이, 실패 이유, 포함/제외 목록, 공여자 수, 평가 층, 재현 seed/매개변수. 원래 세포 순서를 보존한다.

## 모형과 비교 목적

- 단일 모형: 각 유전자 log 평균 count = log(library_size의 학습 기준으로 정규화한 크기) + 하나의 자연 cubic spline(position).
- 분기 모형: 단일 모형에 각 arm의 분기 이후 변화 항을 더하되, 같은 차수의 arm별 계수 합이 0이 되도록 중심화한 대비를 사용한다. 처음 구현은 h=max(position−split_point,0)의 1·2차 항과 sum-to-zero arm 대비의 곱을 사용한다. 이는 각 갈래가 공통 곡선에서 얼마나 벗어나는지를 나타내며, 특정 arm만 공통 곡선으로 특혜를 주지 않는다. 분기점에서 추가 항은 0이므로 공통 구간을 인위적으로 갈라놓지 않는다. arm 재명명·정렬 변경에도 모형 공간과 예측 손실이 같아야 한다. 모든 arm에 제약 없이 항을 추가하여 공통 굴곡의 자유도까지 늘리는 방식은 사용하지 않는다.
- 두 모형 모두 음이항 분포. spline 기준/knots, 정규화 기준, 분산 등 학습 정보는 학습 fold로만 정한다. 평가 자료에서 모형을 다시 적합하지 않는다. 학습 범위 밖 사영은 명시한다.
- 비교값: 동일 평가 세포의 negative log predictive probability를 유전자 내 평균, 유전자 간 평균, 마지막으로 평가 donor 간 동일 가중치 평균한다. 작을수록 예측이 낫다. 단일 손실−분기 손실이 양수면 분기 모형이 더 낮은 예측 손실을 보인다.
- 전체 입력의 독립 donor가 1명일 때 해당 donor 내 cell/replicate 비교로 분류한다. 평가 fold에 donor가 1명인 것만으로 이렇게 분류하지 않는다. 전체 donor 수와 학습·평가 donor의 중복 여부를 확인한다. 평가 세포를 복제해 표본 크기를 부풀리거나 donor 간 신뢰구간으로 변환하지 않는다.
- 작은 학습 집합, 학습에 없는 arm, rank deficiency, 전부 0인 유전자, 비수렴, 무한 예측을 숨기지 않는다. 두 모형의 비교 분모를 맞추고 제외/실패 수를 함께 반환한다. 한 모형의 실패를 다른 모형 승리로 바꾸지 않는다.
- 결과는 branch_preferred_exploratory / single_preferred_exploratory / unresolved 정도의 기술적 비교다. delta=0 또는 비교 불능은 unresolved. 임의 p값/확률/정확도는 생성하지 않는다. 작은 양수 차이를 검증된 생물학적 분기로 해석하지 않는다.

## 공여자 수에 따른 평가

1. 공여자 2명 이상: leave-one-donor-out. 학습/평가 donor ID 교집합은 0. arm이 donor와 완전히 겹쳐 학습에 없는 arm을 예측해야 하면 unavailable로 남긴다.
2. 공여자 1명·독립 replicate_id 2개 이상: replicate holdout. 같은 사람의 기술/채취 반복 평가이며 사람 간 일반화가 아니다.
3. 공여자 1명·별도 반복 없음: 고정 seed로 arm/position 범위를 고려한 cell split(기본 5-fold; 표본이 부족하면 명시적 실패 또는 가능한 더 적은 fold를 기록). `within_donor_cell_split`, `generalization_verified=FALSE`, `exploratory=TRUE`를 반드시 출력한다.
4. 전체 자료에서 만든 graph/position/arm을 공유하면 공여자를 나눠도 `fixed_graph_*_stability`라고 표시한다. 외부 독립 검증으로 승격하지 않는다.
5. 경로 추정 유전자와 평가 유전자가 같은지 provenance에 표시한다. 미지정은 unknown이며 독립 검증을 주장하지 않는다. 이번 시제품은 입력된 고정 후보를 비교하는 도구이고 원자료에서 독립 경로를 새로 만드는 도구는 아니다.

## 구현·검증 순서

- [ ] 실패 테스트: count/metadata ID 불일치, 음수·소수 count, donor/arm 결측, library_size 0, 분기점 범위 오류.
- [ ] 최소 함수: 유효 입력과 공통 구간/분기 후 기준을 검증하고 두 모형을 비교한다.
- [ ] 모의자료: 매끄럽게 휘어진 단일 경로, 발현이 갈라지는 두 arm, arm과 donor가 겹치는 경우, 공여자 1명, 희소/전부 0 유전자.
- [ ] 실제 기대 동작: 단일 모형을 직선으로 제한하지 않음; 분기 데이터에서 보류 예측 차이 확인; 한 공여자 결과를 독립 donor 검증으로 표시하지 않음; 결측/비수렴을 분기 증거로 쓰지 않음; 같은 seed 재현.
- [ ] arm 이름·정렬을 바꿔도 같은 split에서 예측과 손실이 수치 허용 오차 내에서 동일함. 두 모형의 추가 자유도 차이가 공통 곡선의 유연성 증가가 아니라 갈래별 차이에 해당함.
- [ ] 작은 TSV/PNG 예제와 한국어 설명. 원자료 없이 실행 가능한 명령을 제공한다.
- [ ] Codex 실제 코드·테스트·그림 검증 후 필요한 수정만 반영한다. 합의 자체는 검증 근거가 아니다.

실행 예: 고정된 의존성 환경에서 `Rscript tests/R/branch_candidates/run_tests.R`, 이어서 `Rscript R/branch_candidates/run_branch_curve_demo.R`. 공개 재현 기록에는 R/MASS·splines·testthat 버전을 포함하되 개인 컴퓨터의 절대경로는 포함하지 않는다.

현재 생물학적 실자료 성능은 미검증이다. module 패널이나 MPP annotation을 이 구현에서 임의 변경하지 않는다.
