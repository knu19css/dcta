# RNA 중심 DAG 개정안 구현 계획

Goal: 기존 조혈 DAG의 전사체 상태 후보와 실험적 분류를 분리하고 HSC/MPP 두 후보를 표시한, 사용자가 직접 검토할 수 있는 노드·연결표와 그림을 만든다.
Architecture: 원본 76노드/89연결을 동결 복사하고 명시적 대응표를 통해 제안 view를 만든다. 기존 생물학적 근거와 불확실성을 보존한다. 모델 채택은 별도 승인 단계다.
Tech stack: R, base graphics 또는 가용 igraph/ggplot2. 새 설치 없이 재현 가능한 렌더러.
Spec: 사용자 2026-10-16 '다른 브랜치에서 만들고 괜찮으면 병합', '나한테 확인받으라고' 지시; HSC·MPP 별도 후보 유지 결정.

## 범위·역할

- branch: `codex/rna-dag-review-2026-10-16`; base `a88524e7826486cce6702a1ea5c525284482d1c4`.
- Codex 명세/검증/commit/push, Claude Desktop R 구현/작은 데이터표/그림. 여러 작업자가 있으므로 지정된 파일만 수정한다.
- 첫 결과는 제안이지 승인된 RNA ontology가 아니다. 모든 새/변경 항목의 `user_approval`은 pending으로 둔다. HSC·MPP를 별도 후보로 유지하는 결정만 confirmed이며 연결·마커 세부 승인으로 확대하지 않는다.

## 입력 기준

- 원본: 로컬에 보존한 문헌 기반 v1.1.0의 nodes.tsv, edges.tsv, evidence.tsv, edge_evidence.tsv, adult_evidence_audit.tsv 및 필요한 작은 대응 근거표. 공개용 복사본은 개인 경로·내부 세션 기록·실측 세포 정보가 없는지 확인한 뒤 포함한다.
- nodes SHA-256 `438fd18aa47a5c201f993f3f09ee018eb980e3ee06f118438c567fa718a81eca`, edges `ccc67f11019d7107adfa1a1c967f1d6bca1a5b038cfe19253f55f43c82fc8004`.
- 원본 자료 자체는 읽기 전용. 새 branch의 `data/priors/rna_dag_review/source_v1_1_0/`에 위 작은 문헌 기반 TSV만 명시적으로 복사하고 해시 목록을 남긴다. 실측 세포/원 RNA는 복사하지 않는다.

## 출력 계약

- `source_node_map.tsv`: 원본 76 ID 전부, 원본 정의, proposed IDs, 처리(keep/split/reference_only/held), 이유, RNA 구별 검증 상태, 근거 ID, 사용자 승인 상태.
- `proposed_nodes.tsv`: 고유 ID/label/node_kind/view_role/source_node_ids, RNA 상태 근거와 불확실성, approval. HSC_candidate와 MPP_candidate는 별도 행이다. 원본 HSPC_continuum은 provenance/reference summary로 보존한다.
- `source_edge_map.tsv`: 원본 89 ID 전부, source/target, relation, layer/evidence, 새 대응/보류 이유. 기존 연결 ID를 잃지 않는다.
- `proposed_edges.tsv`: RNA 상태 후보 view에서 그릴 연결만. 생물학적 방향 근거와 RNA 구별 근거를 분리하고 source_edge_ids/review_status/approval을 남긴다. HSC/MPP 구별이 불확실한 기존 HSPC 연결은 어느 쪽에 붙일지 unresolved로 보존한다. 두 후보에 자동 복제하지 않는다.
- assay gate/CMP_gate와 기능·taxonomy는 RNA 후보와 시각적으로 분리한 참고 panel로 표시한다. 모든 assay-defined 노드를 생물학적으로 잘못된 노드라고 부르거나 삭제하지 않는다.
- 새 marker signature를 발명하지 않는다. 검토한 근거가 RNA상 독립 유형을 지지하지 않으면 not_assessed 또는 unresolved로 남긴다. 문헌의 표면항원 측정을 RNA 단조 변화 규칙으로 옮기지 않는다.

읽기 전용 원본 검토에 따른 첫 표시안(사용자 승인 대기): RNA 상태 후보 63, assay 정의 7(기존 6 gate와 EryP/BFU-E·CFU-E), 생산 산물 Platelet 1, 보류 2, taxonomy 4로 총 77개 제안 항목이다. HSPC_continuum은 source/provenance 표에 남기며 77개 본체에 또 추가하지 않는다. 분류는 검증된 전사체 유형 목록이라는 뜻이 아니다.

HSPC 출발 19개 원본 관계(E002–E015, E201–E205)는 coarse_ancestry/coarse_model_inference다. 어느 후보가 출발점인지 unresolved로 두고 active edge로 자동 채택하지 않는다. EoP/eNeP/cMoP/MCp/NKP/ILCP/CD56bright·dim 등 단백질·기능 정의 의존 후보에는 배지를 단다. 기존 layer/adult_support를 덮어쓰지 않는다.

## 파일

- `R/rna_dag_review/build_rna_dag_review.R`: 입력 검증, 대응표 생성, 후보/참고 view 구성.
- `R/rna_dag_review/render_rna_dag_review.R`: 원본 대 제안 그림 및 핵심/계통별 그림. 작은 PNG/PDF와 재현 명령.
- `tests/R/rna_dag_review/test_rna_dag_review.R` 및 `run_tests.R`.
- `docs/rna_dag_review/README_ko.md`, `figures/`, `approval_checklist.tsv`.
- `results/rna_dag_review_2026-10-16/implementation_report.md`: 로컬 회신.

## 완료 기준

- [ ] 원본 노드 76개와 연결 89개가 각각 대응표에서 누락 없이 추적된다. 새 ID·참조 무결성·중복·cycle을 검증한다. reference-only 연결을 추론용 edge로 섞지 않는다.
- [ ] HSC/MPP 별도 후보가 존재하고 사용자 미승인 연결은 명확하게 표시된다. 검토 대기 항목이 'validated'로 나가지 않는 테스트를 넣는다.
- [ ] 그림은 한눈에 볼 핵심안과 읽을 수 있는 계통별 상세도를 제공한다. 원본에 있던 것, 변경 제안, 보류를 구분한다. 상자 위치가 분화 속도/시간이라는 오해를 피한다.
- [ ] 사용자 확인 항목을 3개 이내의 구체적 선택으로 정리한다. 전체 개정안 승인 여부와 수정할 노드/연결을 적을 수 있게 한다. 사용자 승인 칸을 자동으로 채우지 않는다.
- [ ] Codex가 실제 생성 표·그림·테스트를 검증한 뒤 GitHub branch에 게시하고 사용자에게 그림과 변경점을 보여준다. 이때 확인 요청을 한다. 승인 전 merge하지 않는다.

예상 명령: 고정된 의존성 환경에서 `Rscript R/rna_dag_review/build_rna_dag_review.R`, `Rscript R/rna_dag_review/render_rna_dag_review.R`, `Rscript tests/R/rna_dag_review/run_tests.R`.
