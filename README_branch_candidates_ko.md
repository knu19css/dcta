# Route–Flow 분기 후보 비교 — 개발 branch

이 branch는 RNA 발현이 하나의 변화 경로로 설명되는지, 여러 갈래가 필요한지 비교하는 R 기능을 개발한다. **운영 모델 반영과 merge는 사용자 검토 전 보류**한다.

- [구현 계획](docs/IMPLEMENTATION_BRANCH_COMPARISON_2026-10-16.md)
- 공개 문서는 일반화한 설계와 재현 방법만 포함하며 내부 조율 기록은 포함하지 않는다.
- 상태(2026-10-18 파일 확인): 테스트 명세와 작은 합성 count 자료 생성기를 작성한 단계다.
- 있는 파일: `R/branch_candidates/simulate_branch_data.R`, `tests/R/branch_candidates/test_branch_curve_comparison.R`, `tests/R/branch_candidates/run_tests.R`.
- 남은 파일: 핵심 `branch_curve_comparison.R`, 실행 예제, 비교 결과와 그림. 핵심 함수가 없어 현재 테스트 전체를 실행할 수 없다. 테스트 통과나 분기 검출 성능을 주장하지 않는다.
- 구현 계약은 공통 곡선에 계수 합이 0인 갈래별 변화 항을 더하는 비교다. 갈래 이름·정렬에 따른 결과 불변성과 공여자별 학습·평가 분리를 테스트 명세에 포함했다.
- 공여자가 한 명이면 동일 공여자 내 재현성만 평가하며, 독립 공여자 검증으로 보고하지 않는다.
- 실제 원자료와 대형 분석 결과는 이 branch 게시 범위에 포함하지 않는다.

RNA 중심 DAG 개정은 별도 `codex/rna-dag-review-2026-10-16` branch에서 검토한다.
