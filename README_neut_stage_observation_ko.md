# 호중구 성숙 단계 관측 — 첫 구현 묶음 (2026-10-19)

2026-10-19 사용자는 DAG 수정을 잠시 멈췄다. 그리고 DCTA 의 중심을 **성숙 단계에 따라 세포군을 구분하고 단계 이름을 자동 주석하는 알고리즘**으로 되돌렸다(spec §14). 첫 대상은 기존 형태 순서의 호중구 성숙이다.

이번 묶음은 그 첫 단계인 관측 범위 점검과 기존 UMAP 표시다. 호중구 형태 6 단계(Myeloblast → PM → MY → MM → BN → SN)에 대해 원 RNA 근거와 실제 관측 연결을 빠짐없이 정리했다. **단계 분류기 완성이나 주석 정확도 개선이 아니다.** 예측 stage_label 은 만들지 않았다.

- 계약: `results/claude_bridge/neut_stage_observation_2026-10-19/codex_instruction.md`
- 계획: `/Users/edemhw/Honglab/CS/04_Claude/14_dcta/docs/superpowers/plans/2026-10-19-neutrophil-stage-observation.md`
- 목표: `.../14_dcta/docs/SESSION_STAGE_ANNOTATION_DIRECTION_2026-10-18.md` §14
- 진행 기록과 Ruling: 같은 bridge 폴더의 `claude_progress.md`. 보고서: `claude_report.md`
- 세션: Claude Desktop 「DCTA 호중구 성숙 단계 구현」, claude-opus-5-5, effort max. 사용자 결정으로 Ultracode 와 고속 모드는 쓰지 않았다.

## 1. 계산

- **입력 고정**
  - 다음 넷의 sha256 이 codex_input_baseline 및 계획 값과 같아야 실행한다. 참고표 `02_nodes.csv` / `03_edges.csv`, `graph_base.npz`(`e63392de…`), `hay_hca_source.h5ad`(`4149b357…`).
  - graph · graph_ids · h5ad 는 O 의 NRV structure 단계 기록(`input_sha256.tsv`)과도 같아야 한다.
  - O 의 NRV p92 · p51 `stage_structure`(cells_lossless · edges · pairs)는 stage manifest hash 로 확인하고 그대로 썼다.
- **reference** (`58_nso_inputs.R`)
  - 6 단계 순서를 고정한다. 단계가 빠지거나, 순서가 틀리거나, 다음 단계로 가지 않는 단계 사이 edge 가 있으면 멈춘다.
  - marker 는 "RNA 프로그램 요약"으로 둔다. 여러 단계가 같이 쓰는 marker 는 지우지 않고 표시한다(LTF · LCN2 = MY · MM, MMP9 = MM · BN, FCGR3B = BN · SN).
  - 단계 사이 변화는 E05–E09 를 core / conditional × up / down 으로 나눈다. 빈칸 = 근거 미지정이다.
  - SN → PB_PMN(E10)은 골수 → 혈액 이동이라 별도 표에만 둔다. 11 node 와 10 edge 는 모두 scope 표에 역할과 함께 남긴다.
- **원 RNA 근거** (`59_nso_evidence.R`)
  - raw/X count 를 chunk 로 한 번 읽는다. x = log1p(10^4 · count / 세포 전체 유전자 library) 이다(R25 와 같다).
  - 참고 유전자는 18 개, 세포는 263,159 개 전부다.
  - 군집 x block 의 n · 평균 logRNA · 검출률을 낸다. p92 가 주결과, p51 이 민감도다. block 은 R24 donor block 42 개다.
- **프로그램**
  - 단계 marker 6 개와 단계 사이 변화 11 개다. 값은 관측 유전자 logRNA 의 동일 가중 평균(서술값)이다.
  - 입력에 유전자가 없으면 unavailable(0 점으로 바꾸지 않음), 있으나 count 가 모두 0 이면 observed_zero 다.
  - 평균이 높은 프로그램을 단계 정답으로 바꾸지 않는다.
- **관측 연결**
  - R24 군집 쌍 표를 다시 계산했고 NRV 저장본과 같다(p92 1,053 쌍 · p51 460 쌍).
  - 경계 세포 = 상대 군집과 실제 edge 가 있는 세포다. 전체 경계 대비와, 같은 block 대비(R26, 한쪽 세포 5 개 이상, R23 기준)를 냈다.
  - 대비 = cluster_b 경계 − cluster_a 경계다. 반대 방향은 부호 반전일 뿐 시간 · 분화 방향 검증이 아니다.
- **B 내부 연속성**: 모든 관측 삼중쌍 A–B–C 에 대해 B 의 내부 edge 연결 성분으로 R24 `nrv_internal_coordinate` 와 같은 상태를 판정한다(작은 그래프 테스트에서 R24 결과와 같음).
- **참조 열**
  - R57 군집 주석(predicted_type · status · warnings)을 썼다. 그 파일의 원저자 열은 읽지 않았다.
  - R54 M 중 결과를 보기 전에 정한 경쟁 후보 7 개(GMP · Neutrophil · Eosinophil · Basophil · Monocyte · CD14/CD16 Monocyte)를 썼다. 선택에는 쓰지 않았다.
- **고정 뒤 표시 · 사후 표**: 계산 표 18 개를 `calc_manifest.json` 으로 고정한 뒤에야 원 UMAP 과 원저자 라벨을 읽었다. 원저자 라벨은 `posthoc_` 표에만 썼다.

## 2. 파일 · 실행

- 코드: `R/branch_biology/58_nso_inputs.R`(reference · 입력 검사 · 정보 분리), `59_nso_evidence.R`(RNA 근거 · 연결 · 내부 연속성), `60_run_nso.R`(실행기), `61_nso_view.R`(그림). 기존 helper(R05 · R17 · R18 · R23–R27 · R34 · R36 · R37 · R46 · R50)는 고치지 않았다.
- 테스트: `tests/R/test_neut_stage_observation.R` 는 17 tests / 244 expectations 이다. 기존 7 suite(106 / 1,385)도 모두 통과했다.
- 테스트 실행: `run-r dcta results/claude_bridge/neut_stage_observation_2026-10-19/claude_run_tests.R`
- 실행:
  ```
  run-r dcta R/branch_biology/60_run_nso.R --data-root /Users/edemhw/Honglab/CS/04_Claude/14_dcta \
    --o-root /Users/edemhw/Honglab/CS/04_Claude/14_dcta/.claude/worktrees/expression-edge-2026-10-04 \
    --reference-dir /Users/edemhw/.codex/visualizations/2026/10/03/01a0ffdf-e36a-70b0-b52b-3389d008cd25/neutrophil_expression/cytoscape \
    --out-dir <N>/results/neut_stage_observation_2026-10-19/attempt_02
  ```
  실제 전체 자료 실행은 약 1 분이다.
- **결과는 `attempt_02`** 다.
  - `attempt_01` 은 처음 실행본이다. 계산 파일 18 개는 attempt_02 와 바이트 단위로 같다.
  - 그림 1 의 옅은 색 구간 때문에 드문 단계가 잘 보이지 않아, 그림만 고쳐 다시 실행했다. attempt_01 은 그대로 보존했다.
- 주요 산출물 (`attempt_02/`)
  - `stage_reference_availability.tsv`: 단계 marker · 변화 유전자의 입력 존재 · 전체 세포 검출 수 · 공유. 관측 가능성을 확정한 표가 아니다.
  - `cell_stage_evidence.tsv.gz`: 263,159 세포의 원 cell_id 순서와 군집. 단계 marker 6 개(검출 유전자 수 · logRNA 평균), 유전자별 logRNA, unavailable 이유, R57 군집 주석 참조를 담는다.
  - `cluster_block_programs.tsv.gz` · `cluster_block_genes.tsv.gz` · `cluster_programs.tsv`: 전체 분모를 보존한 군집 x block · 군집 요약.
  - `observed_pairs.tsv` · `observed_edge_evidence.tsv.gz`: 관측 쌍과 경계 대비(유전자 · 프로그램, 상태 포함).
  - `internal_continuity.tsv.gz` · `internal_cluster_summary.tsv`: 삼중쌍 상태와 군집 내부 성분.
  - `stage_mapping_limitations.tsv`: 공유 marker · 단일 유전자 프로그램 · 조건부 변화 · 범위 · 연결 해석의 한계.
  - `spot_check_cells.tsv`: p92 군집마다 첫 세포의 원 count · library · logRNA.
  - `figures/fig01_neut_stage_evidence_umap.png`: 같은 원 UMAP 의 6 단계 marker 프로그램, 공통 척도.
  - `figures/fig02_observed_connections_umap.png`: 관측 쌍 요약선, 경계 세포, 내부 단절.
  - `posthoc_author_*.tsv`, `checks.tsv`, `input_sha256.tsv`, `source_snapshot.tsv`, `environment.tsv`, `lossless.rds`, `calc_manifest.json`, `manifest.json`, `README_ko.md`

## 3. 결과 (attempt_02, 서술값)

| 단계 | marker (관측 / 목록) | 1 개 이상 검출 세포 | p92 군집 평균 최고 (서술) |
|---|---|---|---|
| Myeloblast | MPO · ELANE · PRTN3 (3/3) | 60,177 (22.9%) | 31 (3.95), 25 (3.30), 47 (2.94) |
| PM | DEFA4 · BPI · DEFA1 · DEFA3 (4/4) | 1,112 (0.4%) | 최고 0.02 |
| MY | LTF · CAMP · LCN2 (3/3) | 5,525 (2.1%) | 84 (0.10), 57 (0.09) |
| MM | LTF · LCN2 · MMP9 (3/3) | 5,077 (1.9%) | 84 (0.10), 57 (0.09) |
| BN | MMP9 · FCGR3B (2/2) | 752 (0.3%) | 최고 0.01 |
| SN | FCGR3B · FPR1 · CXCR2 (3/3) | 13,200 (5.0%) | 27 (0.26), 37 (0.25), 74 (0.23) |

- **유전자별 검출 세포 수**
  - 초기 과립 marker 는 풍부하다: MPO 49,464, ELANE 24,203, PRTN3 21,081.
  - PM–SN marker 는 LCN2 와 FPR1 을 빼면 모두 전체 세포의 0.4% 미만이다: DEFA4 105, BPI 730, **DEFA1 1**, DEFA3 286, LTF 135, CAMP 989, MMP9 532, FCGR3B 220, CXCR2 752.
  - 예외: LCN2 4,426 (1.68%), FPR1 12,318 (4.68%).
  - 변화 전용 유전자: MKI67 37,540, TOP2A 38,055, MYC 64,975, CEBPE 2,935.
- **연결**: 관측 쌍은 p92 1,053 · p51 460 이고 NRV 저장본과 같다.
- **B 내부 연속성**
  - 모든 군집의 내부 연결 성분은 1 개다.
  - p92 삼중쌍 26,839: span_ok 24,635 · no_internal_span 2,204 · internal_disconnected 0.
  - p51 삼중쌍 9,040: span_ok 8,191 · no_internal_span 849 · internal_disconnected 0.
- **같은 block 대비 가용성** (p92 유전자 행 18,954): ok 5,688 · single_block 1,656 · no_eligible_block 6,390 · no_same_block_boundary 5,220.
- **사후** (계산 고정 뒤, 선택에 쓰지 않음, 정확도 아님)
  - Myeloblast 프로그램이 가장 높은 원저자 라벨: GMP-Neut 3.91, GMP-Cycle 3.13, GMP-Mono 2.80.
  - SN 프로그램: CD14 Mono 0.20, CD16 Mono 0.11. 대부분 FPR1 이다.
  - MY · MM 프로그램: Megakaryocyte 0.48 (54 세포), Megakaryocyte Precursor 0.09.
  - 원저자 라벨 중 "Neut" 이 들어간 것은 GMP-Neut 하나다.
- **검증**
  - 실행기 확인 15/15.
  - Claude 독립 검산 14/14: 원 h5 를 직접 읽어 spot 세포 92 x 유전자 18 의 count · library · logRNA 를 대조했다(차이 4.4e-16). 경계 평균은 두 쌍을 직접 다시 계산했다(차이 1.8e-15). 전체 분모 263,159 와 manifest hash 도 확인했다.
  - Codex 의 직접 검산은 아직이다.

## 4. 해석 주의

- 값은 서술값이다. 단계 정답 확률이나 예측 단계가 아니다. 그림의 어떤 패널도 원저자 단계 정답이 아니다.
- **이 자료에서 단계 근거가 고르지 않다.** 초기 과립 프로그램(Myeloblast)은 넓고 강하게 관측된다. 그러나 PM–SN marker 는 LCN2(1.68%)와 FPR1(4.68%)을 빼면 모두 전체 세포의 0.4% 미만에서만 검출된다.
  - 후기 marker 신호 일부는 사후 비교에서 단핵구(FPR1)와 거핵구 계열(LCN2 · LTF) 쪽에 있다.
  - 검출이 적거나 0 인 것은 생물학적 부재의 증거가 아니다. 범위 문헌 검토에도 BoneMarrowMap 의 선별 · 균형 표집에 따른 후기 관측 한계가 기록돼 있다.
- **DEFA1 은 1 세포에서만 검출됐다.** 입력에는 참고표에 없는 DEFA1B 가 따로 있다. 참고표 유전자 · 부호를 바꾸지 않는다는 계획에 따라 넣지 않았다. 바꿀지는 사용자 · Codex 판단 사항이다.
- 공유 marker(MY · MM 의 LTF · LCN2, MM · BN 의 MMP9, BN · SN 의 FCGR3B) 때문에 단계 프로그램은 정의상 겹친다.
- 군집 내부 단절이 0 이므로, 이 분할에서는 내부 연속성이 구별 정보를 주지 않는다. 그래프 연결과 내부 연속은 따로 남겼다.
- 같은 block 대비는 많은 쌍에서 적격 block 이 없다. 그래프 · 분할은 전체 자료로 만들어졌으므로 block 대비는 조건부 서술이다.
- **다음 묶음 전에 정할 것**(사용자 · Codex)
  - 형태 · 분리 실험 단계 주석이나 성숙 호중구가 충분한 평가 자료가 필요한지 (spec §14 1 단계 결정 사항)
  - 참고표 유전자 기호(DEFA1 / DEFA1B)를 다시 볼지
  - 단계 직접 판정 기준모형의 설계
- 이번 묶음에서 하지 않은 것: 단계 확정 · 강제 6 군집 분할 · 새 stage 모형 / 가중치 · DAG 수정 · 12 개 비교 도구 설치 / 실행.
