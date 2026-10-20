# GSE253355 — 최초 수기 마커 패널 (35 유형) 대 현재 slide2 + T (30 유형) 비교 (요약)

- 사용: 1) python (xlrd) scripts/r76_extract_manual_panel_xls.py 로 원본 XLS 형식 추출  2) run-r dcta R/branch_biology/76_gse253355_manual_panel.R --out-dir <새 폴더>.
- 원본: final_unique_gene.xls (sha256 3ae23f8eb1b5…) 두 sheet 389 칸, alias · typo 교정 25 칸, 중복 8 제거 → 381 항목, 양 · 음 충돌 5 쌍은 양쪽 보존 (MPP SLAMF1, CMP KLF1, pDC ITGAX, Pro B cell MS4A1, Cytotoxic CD8 CD28).
- 미측정: current_slide2 3 · manual_original 14. unresolved 이름 (HBA, HBBA, TSAB1, VREB2) 은 원문 유지.
- 점수: R75 와 같은 배경 보정 module score. 대조 유전자는 두 판 합집합 214 개를 빼고 한 번 만들었다 (R75 와 다름: 현재판 유전자 0 / 158 만 같은 대조군; 현재판 1 위가 R75 와 같은 군집 24 / 26).
- current_slide2: 1 위 후보 26 · 동률 0 · 미판정 0 | manual_original: 1 위 후보 26 · 동률 0 · 미판정 0. 표시명 대응 기준 두 판 1 위가 같은 군집 22 / 26.
- 주요 군집 (저자 최다) 현재 / 수기: 6 (Neutrophil) Classical Monocyte / CD14 Mono · 8 (HSC) CDP / CDP · 17 (MPP) Plasma cell / Plasma B · 7 (GMP) GMP / GMP · 12 (Early Myeloid Progenitor) Neutrophil / Neutrophils · 15 (Late Myeloid) Classical Monocyte / CD14 Mono · 16 (Cycling HSPC) CDP / GMP.
- HSC 점수 최고 군집 현재 19 (Fibro-MSC) · 수기 19 (Fibro-MSC); 1 위인 군집 현재 [19] · 수기 [19].
- MPP 점수 최고 군집 현재 5 (CD4+ T-Cell) · 수기 5 (CD4+ T-Cell); 1 위인 군집 현재 [] · 수기 [].
- Neutrophils 점수 최고 군집 현재 6 (Neutrophil) · 수기 6 (Neutrophil); 1 위인 군집 현재 [12] · 수기 [12].
- CD14 Mono 점수 최고 군집 현재 15 (Late Myeloid) · 수기 15 (Late Myeloid); 1 위인 군집 현재 [6,15] · 수기 [6,15].
- 해석: 학습 · 독립 검증이 없어 어느 판이 더 정확하다고 말하지 않는다. 유형 수 (35 vs 30) · 세분화 · 계층 이름 (Mono 와 CD14 / CD16 Mono) 이 달라 단순 비교가 아니다.
- 한계: MSC · 내피 등 패널 없는 계통도 혈액 유형 후보가 붙는다. PM 등 성숙 단계는 다루지 않았다. 저자 라벨은 계산 고정 뒤 사후 대조에만 썼다.
- 표: cluster_side_by_side · cluster_type_scores · core_type_top_clusters · manual_panel_provenance / used / sign_conflicts · markers_unmeasured · control_genes · display_name_map.
- 검사 13 / 13 통과 (checks.csv).
