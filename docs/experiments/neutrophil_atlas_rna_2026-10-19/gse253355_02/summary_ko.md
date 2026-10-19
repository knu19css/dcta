# GSE253355 골수 — 고정 atlas 6 단계 패널의 군집 단위 재적용 (요약)

- 입력: S2_steady_GSE253355.local_copy.rds (sha256 83fd17576c48…, Codex local copy, GEO 원본과 바이트 동일 아님). 세포 82,742 · 유전자 29,452 · donor 12 (donor).
- 패널 유전자 측정: 59 / 60 (정확한 symbol). 없음: CCL3L3 (점수는 고정 분모 10, 관측 유전자 분모는 진단).
- 그래프 · 군집 · 좌표: new_expression_snn (not original) · new Louvain clusters (not author clusters) · 새 DrL 레이아웃(UMAP 아님).
- 군집 26 개 1 위 후보: Myeloblast 11 · PM 4 · MY 0 · MM 2 · BN 2 · SN 7 · 동률 0 · 0 0 (작은 군집 경고 0).
- 단계별 최고 군집 (사후 저자 최다 라벨): Myeloblast 8(HSC) · PM 6(Neutrophil) · MY 12(Early Myeloid Progenitor) · MM 15(Late Myeloid) · BN 25(Ba/Eo/Ma) · SN 1(Adipo-MSC).
- 저자 라벨: Neutrophil 3,857 세포 · donor 12 · 주 군집 6:3855(1.00) 12:1(0.00) 17:1(0.00); Late Myeloid 3,294 세포 · donor 12 · 주 군집 15:2042(0.62) 12:952(0.29) 6:295(0.09) 13:3(0.00) 1:1(0.00) (성숙 호중구 = Neutrophil 라벨).
- 그 세포들 군집의 1 위 후보 비율: Neutrophil [Myeloblast 0.00 · PM 1.00 · MY 0.00 · MM 0.00 · BN 0.00 · SN 0.00]; Late Myeloid [Myeloblast 0.00 · PM 0.09 · MY 0.00 · MM 0.91 · BN 0.00 · SN 0.00].
- atlas PM 패널은 attempt_02 에서 기록한 후기 유전자(MMP9 · FCGR3B · FPR1 · AQP9 · CCL4L2)를 포함한다; Neutrophil 군집의 PM 1 위는 이 한계와 함께 읽는다.
- 연속 단계 최고 군집 연결 (전체 그래프): Myeloblast>PM cell_path · PM>MY cell_path · MY>MM cell_path · MM>BN cell_path · BN>SN cell_path.
- 같은 연결 (과립구 라벨 유도 그래프, 사후): Myeloblast>PM cell_path · PM>MY cell_path · MY>MM cell_path · MM>BN cell_path · BN>SN cell_path.
- 유도 그래프 안 최고 군집 세포 (과립구 라벨/전체): Myeloblast 8 3,381/3,442 · PM 6 4,150/4,152 · MY 12 2,339/2,339 · MM 15 2,043/2,043 · BN 25 8/629 · SN 1 24/9,812. 소수 라벨 세포가 잇는 연결은 성숙 근거가 아니다.
- cluster_path_only = 군집 요약에는 경로가 있으나 실제 세포 그래프에서는 끊김. 연결은 분화 증거가 아니고, 끊김만으로 생물학적 부재를 말하지 않는다.
- 후보는 확률 · 확정 cell type · 분화 판정이 아니다. 저자 라벨은 계산 고정 뒤 사후 맥락으로만 썼다. 패널 유전자는 호중구 밖 계열에서도 발현된다.
- 한계: 새 그래프 · 군집 · 좌표는 원 저자 그래프의 재현이 아니다. donor 동일 가중은 donor x 군집 >= 5 세포만(없으면 결측). 정확도 주장 없음.
- 표: cluster_candidates · cluster_scores · stage_peak_clusters · panel_gene_map · peak_paths · posthoc_*.csv. 그림: figures/fig1 · fig2 (PNG · PDF).
