# DCTA atlas RNA panel experiment (2026-10-19)

- 최신 사용자 승인: 첫 논문 A molecular atlas of human granulopoiesis의 RNA를 재정리·적용하되 별도 브랜치에서만 진행한다. 기존 기본모델·DAG·결과는 수정/병합하지 않는다.
- 이 worktree: codex/neutrophil-atlas-rna-2026-10-19. 이전 M/O/N checkout은 읽기 전용이다. 공개 RNA 다운로드와 이 브랜치만의 선별 commit은 이번 승인 범위다. 원자료/대형결과/인증정보는 commit하지 않는다.
- Pair 유지: Codex 설계·검증·기록, Claude Desktop Code의 DCTA 호중구 성숙 단계 구현 / df42d737-325f-4c9b-bfe9-dd5f4ff6707d가 R 구현·분석·테스트·그림 담당. 전달 전용 relay만 사용한다. 앱 사용자의 최신 설정 Opus5.5/max/일반 속도, Ultracode·subagent 사용 안 함.
- 기존 복제 R05–62 및 테스트는 baseline_sources.json의 바이트를 유지한다. 새 R63 이상·새 테스트·별도 panel/reference/output만 작성한다. 기존 소스의 hash pin을 풀지 않는다.
- run-r dcta 환경만 사용하며 R실행파일/.libPaths()/직접 의존성과 버전을 기록한다. 설치와 다른 R/Conda 혼합 금지.
- 2026-10-19 사용자 승인 예외: GSE253355 실제 UMAP 두 그림을 위해 dcta 환경에 uwot 및 필수 신규 의존성을 추가한다. 기존 설치 패키지 버전은 유지하고 다른 환경을 섞지 않는다. 이 예외는 마커·DAG·기존 분석 결과 변경 승인이 아니다.
- 원 세포263159, p92/p51, 원 graph와 UMAP을 보존한다. 저자 라벨/UMAP/연결성은 사후 평가·표시용이며 source gene 선택이나 cutoff 결정에 쓰지 않는다. 인위 연결·결측 RNA 보충·원 count 수정 금지.
- 골수6단계 Myeloblast/PM/MY/MM/BN/SN 유지. 혈액PMN은 제외·별도 표. bulk 형태확인과 같은 단일세포 사진/RNA 정답을 구별한다.
- 공여자·재측정·재시퀀싱을 구별하며 독립성 확인 안 된 공여자 검증/정확도는 주장하지 않는다.
- 결과 보고 MD는 25줄 이내. 변경·검증·한계만 간결하게.
- 이전 CLAUDE.md/STATE.md의 일반 설정보다 이 실험 계약과 새 계획이 우선한다.
