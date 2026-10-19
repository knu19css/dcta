# Atlas RNA 패널 실험 설계

- 사용자 승인: 첫 논문 RNA 기반 재정리/적용, 별도 브랜치, 기존 모델 보존. 기본 DAG 변경과 병합 없음.
- 논문: A molecular atlas of human granulopoiesis, Research Square v1, DOI10.21203/rs.3.rs-6184761/v1. GSE294330 공개 counts/TPM + GEO sample SOFT; 임의 표면단백질 마커 목록을 RNA 근거로 대체하지 않는다.
- 목적: RNA 인접단계 대비 패널을 재현 가능하게 만들고 현재 세포에서 관측 및 연결성을 시험. 배타 stage 자동주석/단일세포 정확도 검증은 이 실험 결과만으로 주장하지 않는다.
- M/O/N 원 checkout 읽기 전용. 새 checkout 내 baseline 복제본도 hash 보존. 세포263159·p92/p51·원 edge·원 UMAP 고정.
- raw 다운로드 data/external/는 Git에서 제외하고 URL·SHA256만 기록한다. HGNC 공개 표의 정확한 Ensembl ID→승인 symbol 일대일 대응만 사용; 다대일/미대응은 표에 남기고 패널에서 제외. 표적 RNA 값/원저자 라벨/그래프/UMAP으로 유전자를 선택하지 않는다.
- 저자 원고는 4인+생물학적 반복1, GEO에는 A–E 코드5. 어떤 코드가 같은 사람인지 확인되지 않으면 nominal_collection_group으로 표시하고 공여자 독립 검증·p값·정확도·독립 bootstrap CI를 산출하지 않는다. 같은 코드의 재시퀀싱/세부 분획은 같은 그룹에 묶는다.
- SOFT에서 각 matrix column의 GSM·stage·collection group·merge note를 추적. 40mRNA columns만 사용. PMN은 골수6단계와 분리한다. Myeloblast/PM/MY/MM/BN/SN 순서 유지(MC→MY,B→BN,S→SN).
- source expression: 같은 코드×stage의 TPM 중앙값 후 log2(TPM+1). 이를 독립세포로 취급하지 않는다. TPM 차이는 기술적·구성적 영향이 있을 수 있는 서술적 효과이며 DE 검정 결과로 부르지 않는다.
- 인접 대비5개: pair가 있는 코드만 사용. 각 유전자에 뒤-stage minus 앞-stage log2(TPM+1) 차이, 짝수, 중앙값, 부호 일치율, 양 stage TPM을 저장. 양쪽 ≥3코드, 높은 쪽 TPM≥1인 코드≥3, |중앙값 차이|≥0.5, 해당 부호 일치율≥0.75가 후보 기준. 각 부호 effect 크기 순 top20(동률 ENSG 정렬); 전체 대비도 보존.
- node 양성 패널: 각 stage가 존재하는 모든 인접 stage보다 높음(각 대비 동일 기준)을 만족하는 후보를 교집합으로 고른다. endpoint는 한 이웃. 최소 인접 효과 내림차순/ENSG 동률 순 top10. HGNC protein-coding gene만 사용. 후보가 없으면 빈 패널/판정불가로 남기고 목록을 억지로 채우지 않는다. 이 기준은 DCTA 실험 규칙이며 논문 저자의 공식 패널이 아니다.
- 증감은 edge별 RNA 대비로 CSV에 보존. node 양성 마커와 transition up/down을 구별한다. 이웃 양쪽 대비가 다르면 하나의 전역 음성마커로 합치지 않는다.
- 새 panel/reference를 source-only 단계에서 고정한 뒤에만 target H5AD를 읽는다. target gene availability는 선택 후 보고만 한다. missing과 observed zero를 구별. raw count/정규화/증강 없음 정책은 기존 NSO 동일(log1p(1e4*count/all-gene library)).
- 기존 NSO 원표와 비교: stage별 패널 크기·원 유전자 검출·≥1/≥2세포 수·겹침·p92/p51 군집 요약·원 UMAP 표시. 패널 크기가 달라 검출/연결 증가를 정확도 개선이라 부르지 않는다.
- connectivity: 원 cell graph, 같은 block 조건/전체 조건, 기존 R62 순서도달 규칙 그대로. 원 edge만, no zero-hop advance, 전체 분모·불연결·중복 stage 집합 유지. 새 신호 없는 stage는 0과 unavailable 구분.
- posthoc 저자 구성표는 계산/패널 hash 고정 후에만 작성. 형태 단계 정답이 없으므로 classification accuracy 금지.
- 완료: source metadata/contrast/panel CSV + source hash; target 적용표·연결성·그림; 의미 있는 테스트; Codex 검산·Claude 한정 검토; 25줄 이내 보고와 rollback 기준 commit.
