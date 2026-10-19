# 호중구 atlas RNA 패널 실험 (2026-10-19)

branch `codex/neutrophil-atlas-rna-2026-10-19` 전용이다. 기존 기본모델 · DAG · 결과는 바꾸지 않았고 병합하지 않는다.

## 요약

1. GSE294330 bulk RNA atlas 의 골수 6 단계(Myeloblast → PM → MY → MM → BN → SN)에서 두 가지를 source 전용으로 먼저 고정했다.
   - 인접 단계 대비
   - 단계 양성 패널 (단계마다 10 유전자)
2. 고정한 패널을 Hay 골수 263,159 세포에 적용했다. 패널 · 변화 유전자 177 개가 모두 target 에 있다.
3. atlas 패널은 호중구 계열 안에서 단계끼리 비교해 고른 유전자다. 그래서 골수 전체 단일세포에서는 다른 계열에서도 넓게 검출된다(전체 세포의 7–63%).
   - 단핵구 · T 세포 · 적혈구계 · 전구세포
4. 기준선 NSO 패널(단계당 2–4 유전자)보다 검출 · 연결 세포가 크게 늘었다. 그러나 이것은 단계 정확도 향상이 아니다.
5. 예측 단계(stage_label)는 만들지 않았다.

## 계약 · 소유

- 지시서와 계획
  - 지시서: `results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_instruction.md`
  - 설계: `docs/experiments/neutrophil_atlas_rna_2026-10-19/design.md`
  - 계획: `docs/superpowers/plans/2026-10-19-neutrophil-atlas-rna.md`
- Claude 가 새로 작성한 것
  - `R/branch_biology/63_atlas_rna_reference.R` · `64_run_atlas_rna.R` · `65_atlas_rna_view.R`
  - `tests/R/test_atlas_rna_reference.R` · `test_atlas_rna_application.R`
  - 이 문서, `claude_*` bridge 파일, `results/neutrophil_atlas_rna_2026-10-19/`
- 건드리지 않은 것: R05–R62 와 기존 테스트, Codex 파일, Git.
  - R05–R62 와 기존 테스트 21 개는 `baseline_sources.json` 의 hash 와 같다.

## 방법

1. **source (R63).** GEO TPM 40 열을 full SOFT 로 GSM · 단계 · code · 세부 분획 · 병합 표기에 대응했다.
   - 같은 code x 단계의 재시퀀싱 · 세부 분획은 TPM 중앙값 하나로 묶었다. 그 뒤 log2(TPM + 1) 로 바꿨다.
   - code A–E 는 nominal collection group 이다. 논문 4 인 + 반복 1 과의 동일인 관계는 확인되지 않았다. 그래서 공여자 CV · p 값 · 정확도 · CI 는 내지 않는다.
2. **인접 대비 5 개 (뒤 − 앞).** 두 단계가 모두 있는 code 만 쓴다.
   - 후보 조건: 쌍 ≥ 3, 높은 쪽 TPM ≥ 1 인 code ≥ 3, |중앙값 차이| ≥ 0.5, 부호 일치율 ≥ 0.75.
   - 부호별로 top 20 을 고른다. 전체 대비 표(유전자 x 대비)도 그대로 둔다.
   - 이것은 서술적 효과이며 DE 검정이 아니다.
3. **단계 양성 패널.** 이웃한 모든 단계보다 높은 유전자다. 최소 인접 효과 순으로 top 10 을 고른다.
   - 빈 패널을 허용한다. 실제로 6 단계 모두 10 개가 찼다.
   - DCTA 실험 규칙이며 저자 공식 패널이 아니다.
4. **유전자 대응.** HGNC 의 정확한 ENSG 일대일 대응이면서 protein-coding 인 유전자만 쓴다. 나머지는 `source_exclusions.csv` 에 이유와 함께 남긴다.
   - 혈액 PMN 은 별도 표로 두고 대비 · 패널에는 쓰지 않는다.
5. **고정 순서 (R64).**
   1. source hash 를 `source_downloads.json` 과 코드 고정값에 대조한다.
   2. `source/` CSV 를 쓰고 `source_freeze.json` 을 남긴다.
   3. 그 뒤에야 target 을 연다: h5ad · graph · NRV 세포표, 기준선 manifest 고정값.
   4. raw count 를 1 회 읽고(R59: logRNA = log1p(1e4 · count / library)) 연결 · 비교를 계산한다.
   5. `calc_manifest.json` 으로 고정한다.
   6. 그 뒤에만 원 UMAP 그림과 `posthoc_` 저자 표를 만든다.
6. **target 유전자.** h5ad `raw/var/_index` 의 ENSG 로 정확히 대응한다.
   - 없는 유전자는 missing_gene, count 가 모두 0 이면 observed_zero 로 구분한다. 0 으로 바꾸지 않는다.
   - symbol 이 target 과 다른 두 유전자도 ENSG 로 잡았다: CCL3L3 (target 표기 CCL3L1), ADISSP (target 표기 C20orf27).
7. **연결 (R62 규칙 그대로).** 패널 유전자 ≥1 / ≥2 검출 집합을 쓴다. 원 edge 만 쓰고, zero-hop 전진은 없다. 모든 원 edge 와 같은 block edge 를 따로 계산했다.
   - 판정불가 단계(빈 패널 · target 유전자 없음 · 기준보다 적은 관측 유전자)는 0 이 아니라 결측이다.
   - 판정불가 단계 뒤의 도달은 평가하지 않는다.
8. **그림 (R65).** 두 장을 만들었다.
   - fig01: 같은 원 UMAP 위 atlas 6 패널
   - fig02: 기준선 vs atlas 의 패널 크기 · ≥1 검출 · ge1 순서 도달

## 실행 (W 에서, dcta 환경만)

```bash
env -u ANTHROPIC_API_KEY /Users/edemhw/Documents/Codex/2026-10-07/x20/outputs/agent-workflow/scripts/run-r dcta results/claude_bridge/neutrophil_atlas_rna_2026-10-19/claude_run_tests.R test_neut_stage_observation.R test_nso_connectivity_audit.R test_atlas_rna_reference.R test_atlas_rna_application.R
```

```bash
env -u ANTHROPIC_API_KEY /Users/edemhw/Documents/Codex/2026-10-07/x20/outputs/agent-workflow/scripts/run-r dcta R/branch_biology/64_run_atlas_rna.R --data-root /Users/edemhw/Honglab/CS/04_Claude/14_dcta --o-root /Users/edemhw/Honglab/CS/04_Claude/14_dcta/.claude/worktrees/expression-edge-2026-10-04 --baseline-dir /Users/edemhw/Honglab/CS/outputs/dcta_neighbor_state_route_2026-10-16/results/neut_stage_observation_2026-10-19/attempt_02 --baseline-connectivity-dir /Users/edemhw/Honglab/CS/outputs/dcta_neighbor_state_route_2026-10-16/results/neut_stage_observation_2026-10-19/connectivity_audit_02 --out-dir results/neutrophil_atlas_rna_2026-10-19/attempt_02 --chunk 5000
```

exit code: 0 정상 · 2 인자 · 4 불일치 / 검사 실패 · 6 결과 폴더가 비어 있지 않음(덮어쓰지 않는다). 실행 시간은 약 2 분이다.

## 결과 (최종 `attempt_02`)

attempt_01 과 계산 결과가 같다. 다른 것은 시각이 들어간 freeze hash 한 항목뿐이다.

| 단계 | atlas 패널 (target 10/10) | 기준선 크기 | ≥1 검출 기준선 → atlas | ge1 도달 · 모든 edge | ge1 도달 · 같은 block | GMP 라벨 비율* |
|---|---|---|---|---|---|---|
| Myeloblast | ID1 AVP DNTT CD34 ID3 FCER1A GATA2 CRHBP RGS1 CAVIN1 | 3 | 60,177 → 145,503 (55.3%) | 60,177 → 145,503 | 60,177 → 145,503 | 32.2% → 9.8% |
| PM | G0S2 MMP9 PRTN3 DEFB1 CTSG FPR1 AQP9 FCGR3B CCL4L2 ELANE | 4 | 1,112 → 81,926 (31.1%) | 1,065 → 81,923 | 556 → 68,290 | 4.6% → 21.6% |
| MY | DEFA1B DEFA1 DEFA3 TACSTD2 SLPI SCD PRG3 PRG2 ZNF367 BPI | 3 | 5,525 → 19,876 (7.6%) | 4,725 → 19,792 | 825 → 14,867 | 1.2% → 14.1% |
| MM | LTF CAMP MMP8 LCN2 PTX3 CD101 TBC1D8 OLR1 CEACAM1 CDKN1C | 3 | 5,077 → 18,872 (7.2%) | 4,151 → 18,030 | 759 → 8,118 | 0.8% → 3.6% |
| BN | MMP9 CTSH DOK2 ARHGEF17 DMWD TTC38 MAP2K6 SLC24A3 UCP2 VNN1 | 2 | 752 → 165,761 (63.0%) | 101 → 165,761 | 12 → 105,831 | 0.1% → 6.5% |
| SN | CCL4L2 ALPL SECTM1 ZFP36L1 RAB20 CCL4 CCL3L3 IL1B IL13RA1 MME | 3 | 13,200 → 138,811 (52.7%) | 12,165 → 138,811 | 2 → 95,819 | 1.1% → 6.1% |

\* GMP 라벨 비율 = ≥1 검출 세포 중 저자 라벨에 "GMP" 가 들어간 세포의 비율. 해당 라벨은 Early GMP · GMP-Cycle · GMP-Mono · GMP-Neut (21,500 세포)이다.
- 계산 고정 뒤의 사후 서술이다. 이 문자열 규칙은 결과를 본 뒤 정했다.
- 출처: 이번 실행과 기준선 attempt_02 의 `posthoc_author_program_summary.tsv`.

- **순서 도달.** ge1 은 두 그래프 모두 첫 끊김이 없다(기준선도 없다). ge2 도 atlas 는 끊김이 없다(기준선은 PM>MY 에서 끊긴다).
- **사후: 검출 세포의 상위 저자 라벨.**
  - MB: HSC · Early GMP · 순환 전구세포
  - PM: CD14 Mono 21%
  - MY: 적혈구계
  - MM: CD16 · CD14 Mono
  - BN: CD14 Mono · CD4 T
  - SN: CD14 Mono · CFU-E · HSC
- **단일 유전자 검출률이 높다.** 전체 세포 대비 UCP2 47.1% (BN), ZFP36L1 28.0% (SN), IL1B 19.9% (SN), CD34 19.0% (MB).
- **GMP 라벨 세포 중 검출된 비율 (기준선 → atlas).** PRTN3 · CTSG · ELANE 이 atlas PM 패널에 있다.
  - PM: 0.2% → 82.3%
  - BN: 0.0% → 50.2%
  - SN: 0.7% → 39.3%

## 검증

- **테스트.** 38 tests / 435 expectations, 실패 0 (`claude_tests.log`).
  - Codex 기준선: NSO 17/244 + connectivity 8/48
  - atlas: reference 7/61 + application 6/82
  - 합성 실행기 테스트가 다루는 것: missing vs observed zero, 빈 패널 = 판정불가, hash · source 불일치 중지, 덮어쓰기 거절, 263,159 형 분모 · 순서 보존, 원 edge 만, source 고정 전 target 읽기 금지, 계산 고정 전 UMAP · 저자 읽기 금지.
- **돌연변이 확인.** 아래 4 개를 일부러 망가뜨리면 테스트가 모두 실패했다.
  - 판정불가 뒤 평가
  - missing → 0
  - source 고정 전 target 읽기
  - 기준선 hash 생략
- **실행 검사.** 19 개가 모두 TRUE 다(`checks.tsv`).
- **독립 spot-check.** 고정 규칙으로 고른 세포 8 개의 패널 60 유전자를 비교했다.
  - h5ad raw/X 를 hdf5r 로 직접 읽어 480 값을 대조했다. library · logRNA 차이는 0 이다.
  - 단계별 검출 유전자 수 48 개도 같다.
- **Codex 의 Task 1 검산.** 패널 60 행 · 변화 200 행이 원 TPM 과 일치했다(`codex_reference_numeric_checks.*`).

## 해석 · 한계

- **검출 증가는 정확도가 아니다.** atlas 패널은 호중구 계열 안의 인접 단계 대비다. 골수 전체에서는 계열 특이성이 없어, ≥1 검출과 "끊김 없는 도달"이 다른 계열 세포로 채워진다. 기준선 NSO 패널도 PM–SN 의 GMP 라벨 비율이 낮았다(0.1–4.6%).
- **target 의 한계.** 성숙 호중구 · 형태 단계 정답 라벨이 없다. 그래서 단계 정확도를 계산할 수 없고, 계산하지 않았다.
- **PM 패널.** 후기 유전자(MMP9 · FCGR3B · FPR1 · AQP9 · CCL4L2)가 들어 있다. 분류한 PM 분획에 다른 단계가 섞였을 가능성이 있으나 검증하지 않았다.
- **패널 크기.** 기준선 2–4 개, atlas 10 개로 다르다. 직접 비교는 서술용이다.
- **선택 변화 프로그램.** 군집 표(`cluster_programs.tsv`)에만 있다. 해석 · 그림에는 쓰지 않았다.
- **후속 후보** (이번 범위 밖, 설계 결정 필요): 계열 특이성 조건을 패널에 더할지 여부. 이는 "manifold 에 계열 유전자 가중치를 주지 않는다"는 설계 원칙과 함께 검토해야 한다.

## 판단 기록 (Ruling, 자세한 이유는 `claude_progress.md`)

- A1: 선택 변화도 HGNC 일대일 protein-coding 만 쓴다.
- A2: target 은 ENSG 정확 대응, 표기는 HGNC symbol 이다.
- A3: source 표는 CSV(design), target 표는 TSV(NSO 관례)다.
- A4: 판정불가 뒤 도달은 결측이고, 끊김 원인을 구분한다.
- A5: ≥2 기준에서 관측 유전자가 2 개 미만이면 판정불가다.
- A6: 변화 프로그램 표기는 `atlas` 다.
- A7: 유전자별 세포 값은 lossless.rds 에만 둔다.
- A8: 혈액 PMN 은 별도 표 2 개다.
- A9: 전체 그래프 성분표는 다시 만들지 않는다.
- A10: fig02 축 제목 겹침을 고치려고 attempt_01 을 보존하고 attempt_02 를 만들었다.

## 주요 파일 (`results/neutrophil_atlas_rna_2026-10-19/attempt_02/`)

- **source 고정**
  - `source/`: metadata · groups · 전체 대비(gz) · transitions · panels · panel_status · exclusions · HGNC map / summary · PMN 표 2 개 · params · files · RDS
  - `source_freeze.json`
- **target 계산**
  - `target_gene_availability.tsv`, `cell_atlas_evidence.tsv.gz` (263,159 행, 예측 단계 없음)
  - `cluster_block_programs.tsv.gz` · `cluster_block_genes.tsv.gz` · `cluster_programs.tsv` · `program_definitions.tsv`
  - `detection_sets.tsv` · `stage_components.tsv` · `adjacent_pair_edges.tsv` · `ordered_reachability.tsv`
  - `baseline_comparison.tsv` · `reach_comparison.tsv` · `cluster_comparison.tsv` (p92 · p51)
  - `lossless.rds`, `calc_manifest.json`
- **표시 · 사후**
  - `figures/fig01_atlas_panel_umap.png` · `figures/fig02_baseline_vs_atlas.png` · `figure_text.tsv`
  - `posthoc_author_cluster_composition.tsv` · `posthoc_author_program_summary.tsv`
- **기록**
  - `checks.tsv` · `input_sha256.tsv` · `environment.tsv` · `source_snapshot.tsv` · `manifest.json` · `run.log`
