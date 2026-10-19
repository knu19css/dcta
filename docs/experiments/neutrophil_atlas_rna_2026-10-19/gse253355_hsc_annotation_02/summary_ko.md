# GSE253355 — HSC 부터 전체 세포 유형 군집 잠정 후보 (요약)

- 사용: run-r dcta R/branch_biology/75_gse253355_hsc_annotation.R --out-dir <새 폴더>. 입력 sha256: 마커 55d8f2056ffc… · 좌표 1f91bc928f06… · lossless c423ecfe33be… · 객체 83fd17576c48….
- 마커: 검토 후보 CSV (approved_for_model = false) 의 사용 행만. broad 30 유형 (양성 202 · 음성 7), slide2 30 유형 (양성 205 · 음성 74). 제외 75 행은 이유와 함께 panel_markers_excluded.csv.
- 점수: 배경 보정 module score (log1p CP10k, 24 bin, 마커당 대조 50 유전자, seed 20261020), 군집 평균 양성 − max(0, 음성). 0 이하 = 미판정. 확률 · 검증된 분류가 아니다.
- 결과 26 군집: 잠정 26 · 불일치 0 · 미판정 0. 잠정: 0 Plasma cell · 1 Megakaryocyte · 2 NKT cell · 3 Plasma cell · 4 pDC · 5 CD4 T cell · 6 Classical Monocyte · 7 GMP · 8 CDP · 9 CLP · 10 Pre B cell · 11 Erythroblast · 12 Neutrophil · 13 Erythroblast · 14 CLP · 15 Classical Monocyte · 16 CDP · 17 Neutrophil · 18 Plasma cell · 19 HSC · 20 B cell · 21 CDP · 22 Pre B cell · 23 Plasma cell · 24 Pre B cell · 25 Mast cell.
- HSC: HSC — 1 위 군집 broad 19 · slide2 19; 점수 최고 군집 broad 8 (0.13, 저자 HSC) · slide2 19 (0.11, 저자 Fibro-MSC).
- MPP: MPP — 1 위 군집 broad 없음 · slide2 없음; 점수 최고 군집 broad 5 (0.36, 저자 CD4+ T-Cell) · slide2 5 (0.36, 저자 CD4+ T-Cell).
- 그림 2 주요 군집: 6 Classical Monocyte (broad Classical Monocyte · slide2 Classical Monocyte) · 7 GMP (broad GMP · slide2 GMP) · 8 CDP (broad CDP · slide2 CDP) · 12 Neutrophil (broad Neutrophil · slide2 Neutrophil) · 15 Classical Monocyte (broad Classical Monocyte · slide2 Classical Monocyte) · 16 CDP (broad CDP · slide2 CDP) · 17 Neutrophil (broad Neutrophil · slide2 Neutrophil).
- PM: 기존 6 단계 점수 · 후보는 바꾸지 않았다. 군집 12 기존 최고 단계 후보 MM (unvalidated) · 군집 17 기존 최고 단계 후보 PM (unvalidated). PM 문제는 여전히 해결되지 않았고 성숙 단계는 확정하지 않는다.
- 한계: MSC · 내피 · 골모세포 등은 패널이 없어 이 군집들도 혈액 유형 중 가장 높은 것이 나오거나 미판정이 된다. 미판정이 이 계통을 배제한다는 뜻이 아니다.
- 한계: 마커 목록 (검토 후보) · 군집 해상도 · 대조 표본에 의존한다. bootstrap 없음. 저자 라벨은 계산 고정 뒤 사후 대조에만 썼다.
- 표: cluster_provisional_annotation · cluster_type_scores · type_top_clusters · panel_markers_used / excluded · control_genes · gene_bins. 그림: figures/fig1 · fig2.
- 검사 12 / 12 통과 (checks.csv).
