# Neutrophil Atlas RNA Implementation Plan

> **For agentic workers:** use executing-plans inline; no subagents. Pair ownership below is already approved.

**Goal:** Atlas RNA 기반 패널을 별도 브랜치에서 적용하고 기준선과 비교한다.
**Architecture:** source 전용 패널 산출을 먼저 고정한 뒤 target 관측/연결성에 전달한다. 기존 R05–62는 바이트 보존, 새 R63 이상만 작성한다.
**Tech Stack:** R4.5.2, run-r dcta, Matrix/hdf5r/igraph/digest/jsonlite/testthat; 새 설치 없음.
**Spec:** docs/experiments/neutrophil_atlas_rna_2026-10-19/design.md

## Global Constraints
- AGENTS.md/Spec의 읽기전용 source, 263159분모, 순서6단계, source-only selection, 원그래프, 독립성/정확도 제한 유지.
- Claude owns 새 R·테스트·분석·그림·claude_* bridge. Codex owns design/plan/input records·codex_* 검증·Git.

## Review Focus
- 중복 technical columns와 donor-stage missingness는 paired 분모를 바꾼다.
- ENSG와 승인 symbol의 비일대일 대응/없는 target 유전자를 조용히 제거하지 않는다.
- interior stage에 단독 peak가 없으면 빈 panel이 정답일 수 있다.
- 파라미터·유전자 순위는 target RNA/라벨/연결을 보기 전에 고정한다.
- 빈 panel/부족한 target genes/더 큰 panel의 검출 증가를 정확도 증가로 오해하지 않는다.

### Task 1: source 전용 RNA 대비와 panel
**Files:** R/branch_biology/63_atlas_rna_reference.R, tests/R/test_atlas_rna_reference.R.
**Interfaces:** atlas_read_source(data_dir) → expression/metadata/map; atlas_build_reference(source) → all_contrasts, selected_transitions, panels, exclusions, nominal-group summary.
- [ ] 작은 손계산 예제에 중복 column, missing stage, TPM0, mixed sign, ambiguous HGNC mapping, empty interior panel을 포함하여 실패 테스트를 쓴다.
- [ ] spec 고정 cutoff/정렬을 구현하고 테스트 통과 확인. 실자료 metadata와 대비·panel을 저장·hash 고정한다.
- [ ] 4인/5코드 관계 근거가 없으면 독립 CV를 생략하고 미확인으로 기록한다.

### Task 2: 고정 panel target 적용과 원 edge 연결
**Files:** R/branch_biology/64_run_atlas_rna.R, 65_atlas_rna_view.R, tests/R/test_atlas_rna_application.R; README_neutrophil_atlas_rna_ko.md.
**Interfaces:** 별도 runner 입력 source_dir/M/O/Nbaseline/out; atlas 고정 reference → cell/cluster program evidence + R62형 connectivity, PNG.
- [ ] missing vs observed zero, 빈 panel, reference/hash 불일치, 263159/order 보존, 원 edge만 사용, 기존 output 덮어쓰기 거절을 테스트한다.
- [ ] 기존 불변 primitive를 재사용하고 새 소스에서만 조립한다. 필요하면 번호66을 추가한다. 기존 NSO_PINNED/NCA_PIN을 완화하지 않는다.
- [ ] source-only freeze 후 raw target을 순차 chunk로 읽어 원 count 기반 근거·allcells 분모를 계산한다.
- [ ] baseline manifest 비교표와 panel 크기·연결성·원 UMAP 그림을 생성한다. 원저자 라벨은 posthoc만.

### Task 3: 검증과 되돌리기
- [ ] Claude 테스트·실행 command/environment/source/input hash·25줄 이내 report와 진행 bridge를 남긴다.
- [ ] Codex가 source 대비 일부를 원표로 검산하고 cell-count/edge/path/분모 및 baseline file hashes를 확인한다.
- [ ] Claude가 실제 결과와 Codex 해석을 한정 교차 검토한다. 실행 가능한 오류만 수정·재검증한다.
- [ ] Codex가 이 브랜치의 선별 commit과 baseline commit을 기록한다. 기존 branch와 결과는 보존, merge 없음.

Self-review: spec→각 task 연결, 고정 cutoff/함수 경계/입출력·실패 사례 확인. 사용자 이미 구현·분리를 승인했으므로 추가 승인 대기 없이 지정 앱 구현으로 진행.
