# RNA 중심 조혈 DAG 검토안 (사용자 검토 대기)

문헌 기반 조혈 DAG v1.1.0 의 76 노드와 89 연결을 **하나도 빼지 않고** 대응시킨 검토안이다. 두 가지를 나눴다.

- RNA 상태 후보
- 실험적 분류: assay gate, 기능 · 배양 정의, taxonomy

HSC 와 MPP 는 별도 후보로 표시했다.

- **승인된 RNA ontology 가 아니다.** 모든 항목의 `user_approval` 은 `pending` 이다.
- 사용자 승인 전에는 다른 branch 에 병합하지 않는다. 운영 모델에도 채택하지 않는다.
- RNA 로 상태가 구별되는지는 **평가하지 않았다.** 해당 열은 `not_assessed` 이고, HSC · MPP 는 `unresolved` 이다.
  - 표의 근거 ID 는 원본의 방향 · 계통 관계 근거다. RNA 구별 검증 근거가 아니다.

## 그림 (`figures/`, PNG · PDF)

| 파일 | 내용 |
|---|---|
| `rna_dag_original_vs_proposed` | 원본 node_kind → 제안 역할, 원본 support_class → 제안 처리. HSPC_continuum 하나에서 HSC · MPP 별도 후보로 바뀐 모습 |
| `rna_dag_core` | 제안 view 전체. 계통 묶음 안에 RNA 상태 후보 · 산물 · 보류 후보가 있다. assay 정의와 taxonomy 는 아래 참고 panel 에 있다 |
| `rna_dag_lineage_<group>` | 계통별 상세 6 장. 노드에는 ID · 원본 이름 · 배지가 있다. 왼쪽 회색 열 = assay 정의, 오른쪽 청록 열 = taxonomy, 오른쪽 표 = 그림 속 원본 연결 |

계통별 그림은 다음 6 장이다.

- `mk_ery`
- `granulocyte`
- `mono_dc`
- `b`
- `t`
- `nk_ilc`

**읽는 법**

- 세로 층은 원본 연결의 위상 순서일 뿐이다 (앞선 상태가 위). **분화 시간이나 속도가 아니다.** 가로 위치는 상자가 겹치지 않게 놓은 배치다.
- 선 종류
  - 검은 실선: 유지 제안 연결. 원본 방향 근거는 그대로이고 승인 대기다.
  - 파란 점선: 출발점 미해결.
  - 빨간 점선: 원본 보류.
  - 회색 선: assay 정의 끝점 참고 연결.
  - 회색 점선: 산물 생성, 또는 gate 포함 관계.
- 파란 점선은 HSC · MPP 묶음 상자의 **아래 변 한 점**에서 나온다. 어느 한 후보 상자에서 나오지 않는다. 즉 출발점을 정하지 않았다는 뜻이다.
- 꺾인 선은 다른 상자를 피해 돌아간 경로다. 의미는 직선과 같다.
- 주황 표시는 원본 정의가 단백질 표면항원 · 기능/배양 측정에 의존하는 후보다.

## 표 (`data/priors/rna_dag_review/proposal/`)

| 표 | 행 | 내용 |
|---|---|---|
| `source_node_map.tsv` | 76 | 원본 노드 → 처리 · 제안 ID · 이유 · RNA 구별 상태 · 정의 / 연결 근거 ID |
| `proposed_nodes.tsv` | 77 | 제안 항목 (`view_role`), 원본 ID, 불확실성, 배지 |
| `source_edge_map.tsv` | 89 | 원본 연결 → 처리 · 표시 panel · 제안 연결 ID · 추론 사용 여부 · 이유. 원본 layer · adult_support 등은 그대로 |
| `proposed_edges.tsv` | 62 | RNA 상태 후보 view 에 그리는 연결만 (유지 제안 49 + 출발점 미해결 13). 방향 근거와 RNA 구별 근거를 다른 열에 둔다 |

처리 기준과 개수는 다음과 같다.

- **노드**
  - `keep` 61: RNA 상태 후보로 유지한다.
  - `split` 1: HSPC_continuum 을 HSC_candidate · MPP_candidate 로 나눈다.
  - `reference_only` 12
    - assay 정의 7: 6 gate 와 EryP (BFU-E / CFU-E)
    - 산물 1: Platelet
    - taxonomy 4
  - `held` 2
  - 제안 77 항목의 역할 구성: RNA 상태 후보 63 (HSC · MPP 포함) · assay 정의 7 · 산물 1 · 보류 2 · taxonomy 4.
  - HSPC_continuum 은 77 항목에 넣지 않는다. 원본 대응표에 provenance 로 남긴다.
- **연결**
  - `keep` 49: 두 끝이 모두 RNA 상태 후보다.
  - `unresolved` 19: 원본 HSPC_continuum 출발 coarse 관계다.
    - E002–E015 와 E201–E205 가 해당한다.
    - 두 후보에 복제하지 않는다.
    - 13 개는 RNA view 에 그린다. 6 개는 도착점이 assay 정의라 참고 panel 에만 그린다.
  - `reference_only` 13: assay 정의 끝점 12 + Platelet 생산 1.
  - `held` 8: 원본 support_class 가 held 다.
  - 추론용 후보 (`candidate_pending_approval`) 는 유지 제안 49 개뿐이다. 나머지는 `not_for_inference` 다.
- HSC_candidate 와 MPP_candidate 사이 연결은 원본에 없다. 그래서 그리지 않았다. HSC → MPP 방향도 자동 승인하지 않는다.
- 배지 `protein_or_function_defined` 는 다음 8 개에 달았다.
  - EoP, eNeP, cMoP, MCp, NKP, ILCP, NK_CD56bright, NK_CD56dim
  - 원본 이름에 표면항원이 들어 있는 후보도 배지 대상인지 검토가 필요하다. 예: Mono_classical / intermediate / nonclassical 의 CD14 · CD16, cDC1 의 CD141, cDC2 의 CD1c.
  - 이번에는 새로 배지를 달지 않았다.
- 새 마커 · 생물학적 근거는 만들지 않았다.

## 사용자 확인 (`approval_checklist.tsv`)

| 항목 | 질문 | 제안 |
|---|---|---|
| A1 | 개정안 전체를 받아들일지 | 위 77 항목 · 89 연결 대응 |
| A2 | HSC/MPP 출발 연결 19 개의 출발점 | 모두 미해결 유지. 어느 후보에도 붙이지 않음 |
| A3 | assay 정의 · 산물 · taxonomy 분리 | 삭제하지 않고 참고 panel 로 분리 |

- 각 항목에는 (a) 그대로 / (b) 수정 / (c) 되돌림 선택지가 있다.
- 선택은 `user_choice` 에 적는다. 고칠 노드 · 연결 ID 와 내용은 `user_notes_ko` 에 적는다.
- 빌더는 이 칸을 채우지 않는다.
- 이미 적힌 검토 기록이 있으면 다시 실행해도 덮어쓰지 않는다. 아무것도 쓰지 않고 멈춘다. 이는 `user_choice` · `user_notes_ko` 가 비어 있지 않거나 표의 `user_approval` 이 `pending` 이 아닌 경우다. 새 결과는 다른 `--out` / `--checklist` 경로에 쓴다.

## 원본 복사본 (`data/priors/rna_dag_review/source_v1_1_0/`)

복사 방식은 `SHA256SUMS.tsv` 에 파일마다 적혀 있다.

- **바이트 그대로 복사**: nodes, edges, edge_evidence, taxonomy, gate_membership.
  - nodes · edges 의 SHA-256 은 v1.1.0 값과 같다. 빌더가 읽을 때 다시 확인한다.
- **`evidence_public_subset.tsv`**: 원본 evidence 표에서 검토자 기록 열(`reviewed_by`)을 뺐다.
  - `review_level` 은 문헌 검토 수준 값만 남겼다.
  - 나머지 열 · 행 · 값은 원본 그대로다.
- **`adult_evidence_audit.tsv` 는 공개 복사본에 넣지 않았다.** 원본 SHA 만 기록했다.
  - 자유 서술에 내부 작업 메모가 있다.
  - 이 검토의 빌더는 이 표를 쓰지 않는다. 연결의 adult_support 는 edges 에 있다.
- 실측 세포 자료 · 원 RNA 는 복사하지 않는다.

## 재현

R 4.5 와 igraph · digest · testthat 가 필요하다. 그림은 base graphics 로 그린다. 저장소 루트에서 실행한다.

```sh
Rscript R/rna_dag_review/copy_source_v1_1_0.R --source-dir <v1.1.0 표 폴더>   # 원본 복사 (이미 있으면 같은 내용인지만 확인)
Rscript R/rna_dag_review/build_rna_dag_review.R                               # 대응표 4 개 + approval_checklist.tsv
Rscript R/rna_dag_review/render_rna_dag_review.R                              # 그림 16 개 (--lang en 으로 영어)
Rscript tests/R/rna_dag_review/run_tests.R                                    # 테스트
```

렌더러는 저장된 표가 원본에서 새로 만든 표와 다르면 멈춘다.

macOS 에서는 시스템 한글 글꼴(Apple SD Gothic Neo)을 이름으로 등록해 쓴다. 다른 환경에서 한글이 깨지면 `--lang en` 을 쓴다.

## 한계

- 이 개정안은 **표시 · 분류 제안**이다. 원본 문헌 근거의 등급 · 연령 범위를 다시 평가하지 않았다.
- RNA 로 상태가 구별되는지 검증하지 않았다.
- 계통 묶음은 원본 module 을 그대로 썼다. 묶음 안 층 순서는 원본 연결에서 계산했다.
- 원본의 species · context · caveat 같은 열은 대응표에 복사하지 않았다. 원본 파일과 연결 ID 로 이어서 본다.
