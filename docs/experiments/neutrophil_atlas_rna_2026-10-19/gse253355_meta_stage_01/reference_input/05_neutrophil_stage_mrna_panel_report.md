# 고전적 호중구 분화 단계의 mRNA 패널 — 문헌 기반 메타분석 보고서

작성 2026-10-20 · 대상 원고: `01_docs/04_neutrophil_review_manuscript.md` · 결과 통합본: `04_results/12_neutrophil_stage_mrna_panel.xlsx`

## 1. 핵심 결과

- **자료**: 체계적 검색(PubMed 50, GEO 163) → 데이터셋 33건 정밀 확인 → **사람 골수 정렬 전사체 12개** 포함(제외 21건, 사유 기록). 데이터셋별 limma 후 요약통계만 결합(원자료 병합·배치보정 없음).
- **고신뢰 변화**: 직접 검정한 8개 이행(GMP→MB는 추론)마다 견고한 HC 유전자 44–1,663개(데이터셋 하나씩 빼기에도 95–100% 유지). 합의 곡선은 개별 데이터셋과 잘 맞음(기능 패널 Spearman ρ 중앙값 0.91).
- **문헌 대조**: 검증된 근거 노트에서 유전자 단위 진술 943건(사람 601·생쥐 342, 43개 출처, 유전자 213개)을 추출해 모든 인용문을 원문 줄과 대조. 판정 가능한 **사람 mRNA 진술 174건 중 148건(85%) 지지, 반대 3건(1.7%)** — 반대 3건은 소수 경로(IFIT3, 인터페론 NCP 경로)와 군집 이름 오프셋(DUSP1·RGS2)으로 설명된다.
- **기능 패널 69개**: 등급 **A 38**(메타분석 + 사람 문헌) · **B 29**(메타분석; 주로 HSC/MPP 유전자 — 원고 참고문헌이 호중구 중심이라 문헌 대조 대상이 적음) · C 2(THY1, FCGR1A) · R 0. 유전자 전체로는 A 76, B 9,236.
- **단계 순서(메타 합의)**: HSC 유전자(HLF·MECOM·CRHBP) 소실(MPP→CMP·CMP→GMP) → 아주로필 과립(MPO·ELANE·AZU1·PRTN3·CTSG) 증가(CMP→GMP), MB에서 정점 → CD34 소실·디펜신·BPI·초기 특이과립(CEACAM8·OLFM4·ARG1) 시작(MB→PM) → 특이과립(LTF·CAMP·LCN2)·CD11b·NADPH 산화효소·CD101 증가, 아주로필·KIT 감소(PM→MY) → 세포주기 이탈·MMP9·FPR1·CXCR1 증가, CD49d 소실(MY→MM) → CD16b·CD10·CR1·CSF3R 증가, MKI67·MYB·C/EBPε 감소(MM→BN) → CD10·CXCR2·CSF3R 최대, 특이과립 mRNA 감소(BN→SN).
- **기존 DCTA 단계 패널 점검**(GSE294330 단독): 60개 중 49개가 ±1단계로 일치. 다만 **PM 패널 10개 중 6개(G0S2·MMP9·FPR1·AQP9·FCGR3B·CCL4L2)는 메타 합의상 BN/SN 정점** — GSE294330 PM 분획의 성숙세포 약 1% 혼입 신호로 보인다(`04_results/10_meta_panel/10_dcta_atlas_panel_check.tsv`).
- 그림: `04_results/11_figures/11a`(참조 그림 형식 쐐기 그림) · `11b`(mRNA 패널 히트맵) · `11c`(이행별 증가/감소 유전자) · `11d`(데이터셋 구성) · `11e`(forest plot).

## 2. 질문과 단계 정의

- 질문: 사람 골수에서 고전적 단계 HSC → MPP → CMP → GMP → Myeloblast(MB) → Promyelocyte(PM) → Myelocyte(MY) → Metamyelocyte(MM) → Band(BN) → Segmented(SN)로 넘어갈 때 **여러 독립 데이터셋에서 일관되게** 증가·감소하는 mRNA는 무엇이고, 문헌 진술과 맞는가.
- 단계 대응: 각 연구의 저자 표기를 따른다. 예외는 두 가지다.
  - MB: 정렬 분획 이름이 "Myeloblast"인 경우(GSE294330)와 NCP1–4(GSE274002; 저자 서술 "likely myeloblasts")를 MB로 둔다. NCP5/6은 MB 후기로 따로 둔다.
  - BLUEPRINT(Grassi 2018·Hoogendijk 2019): 원고 Box 2의 형태학 대응(P/M = PM, "MM" = MY, "BN" = MM, "SN" = BN+SN)을 기본으로 쓰고, 저자 표기는 민감도 분석으로 확인했다.
- 두 단계를 합친 분획(예: MY/MM, BN/SN, MB/PM, CD34⁺ 전체)은 **블록 대비**(B1 HSC/MPP→CMP/GMP, B2 GMP→PM, B3 PM→MY/MM, B4 MY/MM→BN/SN, B34 PM→BN/SN)에만 쓴다.

## 3. 방법

### 3.1 체계적 검색 (2026-10-20)

| DB | 검색식 ID | 건수 |
|---|---|---|
| PubMed | P1 호중구 단계 × 전사체 × 골수 × 사람 | 43 |
| PubMed | P2 HSC × CMP/GMP × 전사체 × 정렬 × 사람 | 7 |
| GEO DataSets | G1 전골수구·골수구·후골수구·호중구 전구체 | 67 |
| GEO DataSets | G2 골수 호중구 성숙/분화 × 정렬 | 35 |
| GEO DataSets | G3 CMP/GMP | 43 |
| GEO DataSets | G4 단일세포 × 골수 호중구 발달 | 54 |

검색식 전문·결과는 `04_results/06_meta_search/` (`06_search_log.tsv`, `06_pubmed_records.tsv`, `06_geo_records.tsv`, `geo_brief/`).

### 3.2 선별 (PRISMA 요약)

- 식별 213건(PubMed 50, GEO 163) → 제목·요약 선별 후 데이터셋 단위 33건을 표본 구성까지 확인(`03_refs/08_neutrophil_meta_study_registry.tsv`).
- 포함 기준: 사람 **골수** 1차 세포, 정상 공여자(또는 대조군), 연속한 고전 단계 2개 이상을 FACS·형태로 정렬한 분획, 공개 발현 행렬, 분획당 반복 2개 이상.
- 제외 20건(적합성): 제대혈·말초혈 HSPC 5, 체외 분화·자극 4, 생쥐/혈액만·세포주 2, mRNA 아님 3, 연속 단계 없음 3, 설계(단일세포·미니벌크) 3.
- 자료 점검 제외 1건: **GSE19556 (Theilgaard-Mönch 2005)** — GEO 처리값에서 하우스키핑 유전자가 상위 분위에 없고(B2M 49%, EEF1A1 41%; 나머지 12개 데이터셋은 모두 95% 이상) 표지 유전자 양상이 생물학과 반대(예: 혈액 PMN에서 MPO > PM)여서 probe–값 어긋남으로 판단. 원시 CEL 재처리는 하지 않았다.
- **포함 12건** — 표와 그림 `11d`. 문헌 진술 추출·검증: `03_refs/07_neutrophil_literature_stage_claims.tsv`, `02_scripts/10a_neutrophil_literature_claims_validate.py`(943행 모두 원문 줄 대조 통과).

| 데이터셋 | 연구 | 플랫폼 | 단계 (n) |
|---|---|---|---|
| GSE42519 | Rapin 2014 / Bagger 2016 | GPL570 | HSC 4 · MPP 2 · CMP 3 · GMP 5 · PM 6 · MY 2 · MM 3 · BN 4 · SN 3 |
| GSE274002 (+GSE175880) | Signoretto 2025 + Calzetti 2022 | RNA-seq FPKM | MB 24 · MB 후기 6 · PM 4 · MY 5 · MM 5 · BN 5 · SN 5 |
| GSE74246 | Corces 2016 | RNA-seq | HSC · MPP · CMP · GMP 각 4 (공여자 짝지음) |
| GSE294330 | Hesse/Klein (GEO 2026) | RNA-seq | MB 4 · PM 4 · MY 5 · MM 4 · BN 5 · SN 3 |
| GSE12662 | Payton 2009 | GPL570 | CD34⁺ 5 · PM 5 · BN/SN 5 |
| GSE19599 | Andersson 2010 | GPL570 | CMP · GMP · PM · MY 각 2 |
| GSE49883 | Rönnerblad 2014 | Illumina | CMP · GMP · PM · BN/SN 각 5 (짝지음) |
| GSE63270 | Jung 2015 | U133Plus2 (Entrez CDF) | HSC · MPP · CMP · GMP 각 7 (짝지음) |
| GSE77439 | Miyawaki 2017 | Agilent | HSC · CMP · GMP 각 4 |
| GSE98946 | Schwarzer 2017 | RNA-seq | MB/PM 4 · MM 3 · SN 2 |
| GSE153263 | Dinh 2020 | RNA-seq | eNeP 3 · PM 3 · MY/MM 3 · BN/SN 3 |
| BLUEPRINT | Grassi 2018 / Hoogendijk 2019 | RNA-seq FPKM (3,455 유전자) | PM · MY · MM · BN/SN 각 3 |

### 3.3 데이터셋별 분석과 결합

1. 데이터셋마다 따로 정규화(count → TMM·voom, FPKM → log2(x+1), 마이크로어레이 → 저자 처리값 log2)하고 유전자 기호를 HGNC 공식 기호로 맞췄다. 원자료 병합·배치보정은 하지 않았다.
2. limma(voom 또는 limma-trend, robust eBayes). 공여자 짝지음이 기록된 데이터셋만 공여자를 고정효과로 넣었다(GSE274002, GSE74246, GSE49883, GSE63270, GSE98946). GSE294330의 공여자 코드는 단계 간 동일인 여부가 확인되지 않아 짝짓지 않았다.
3. 대비별 발현 하한: RNA-seq은 log2 ≥ 1, 마이크로어레이는 유전자별 최대 단계 평균의 중앙값(배경 수준). 비교하는 두 단계 모두 하한 미만이면 그 데이터셋은 그 유전자 대비에서 뺀다.
4. 결합: 가중 Stouffer Z(가중치 √n) → BH FDR, DerSimonian–Laird 무작위효과 Hedges' g·log2FC, I².
5. **고신뢰(HC) 호출**: 데이터셋 k ≥ 2, 메타 FDR < 0.01, 개별 유의(FDR<0.05) ≥ 2, 반대 방향 유의 0, 방향 일치 100%(k ≥ 5면 80%), 중앙 |log2FC| ≥ 1. 데이터셋을 하나씩 빼고 다시 계산해 MOD 이상이 유지되면 "견고(robust)".
6. **GMP→MB (T4)**: 두 단계를 함께 가진 데이터셋이 없다. B2(GMP→PM)가 견고한 HC이고, MB→PM(T5)에 같은 방향 변화가 없으며, 증가는 MB 데이터 모두에서 이미 발현될 때만 "MB 이전 변화"로 추론(△/▽). PM 분획의 후기 세포 혼입 때문에 생기는 가짜 증가(예: MMP9·FCGR3B)는 이 조건에서 걸러졌다.
7. **합의 곡선**: 이행 T1–T9의 결합 log2FC를 누적(HC·MOD, 또는 k = 1에서 그 데이터셋이 유의한 경우만; 나머지 0). 각 단계를 가진 데이터셋 과반에서 발현 하한 미만이면 그 단계를 0으로 두고, 선형 상대 수준 h(단계 최대 = 1)로 바꿨다. 쐐기 그림·히트맵의 모양은 이 값이다.

### 3.4 문헌 근거 대조와 등급

- 출처: 원고 프로젝트의 검증된 근거 노트(`03_refs/02–05_neutrophil_evidence_*.md`; 각 진술에 [stated]/[inferred] 표시와 원문 인용)와 원고 본문의 [@Key] 진술. 참고문헌 key는 모두 `03_refs/01_neutrophil_refs_verified.csv`(PubMed 검증)에 있는 것만 썼다.
- 추출 단위: (유전자, 진술, 출처) 한 행. 단백질·표면 항원 이름은 유전자 기호로 바꿨다(CD16 → FCGR3B, CD10 → MME 등). 단계는 저자 표기를 고전 코드로 옮기고 근거를 함께 적었다. 표: `03_refs/07_neutrophil_literature_stage_claims.tsv`.
- 추출 검증: 행마다 인용문이 지정한 파일·줄(±3줄)에 실제로 있는지, key가 검증 목록에 있는지, 어휘가 규칙대로인지 기계적으로 확인했다.
- 판정(메타 합의 곡선 기준, 문 이름 오프셋을 고려해 ±1단계 허용): peak S → 곡선 최고 단계가 S±1이면 지지, 3단계 이상 떨어지고 h(S) < 0.25면 반대 · onset S → S±1 안에서 h가 처음 0.25를 넘으면 지지 · rises/falls A→B → log2 곡선 차이 ≥ 1 같은 방향이면 지지, 반대 방향 ≥ 1이면 반대 · present/high → 이웃 포함 최대 h ≥ 0.25면 지지 · absent/low → h < 0.1이면 지지. 혈액 PMN 진술은 골수 곡선으로 판정하지 않는다.
- 등급(유전자): **A** = 견고한 HC 호출 + 사람 mRNA 진술 지지 ≥ 1(또는 사람 진술 지지 출처 ≥ 2) + 사람 mRNA 진술 반대 0 · **B** = 견고한 HC 호출만(문헌 없음 또는 교차 모달리티·생쥐만), 또는 문헌 지지 + MOD 호출 · **C** = 근거 약함 · **R** = 견고한 HC와 사람 mRNA 진술이 반대(검토 필요).
- 단계 패널: 이행마다 견고한 HC 유전자(T4는 추론) 중 단백질 암호화 유전자를 등급 A → B, |Stouffer Z| 순으로 증가·감소 각 10개. 호중구 구간(GMP→SN)에서는 다른 계통 표지(적혈구·거핵구·호산구/호염기구/비만세포·림프구) 유전자를 정렬 분획 혼입 신호로 보고 제외했다.

## 4. 결과

### 4.1 이행별 고신뢰 변화 (유전자 수)

| 이행 | k(최대) | HC 증가 | HC 감소 | 견고 |
|---|---|---|---|---|
| T1 HSC→MPP | 3 | 37 | 8 | 44 |
| T2 MPP→CMP | 3 | 159 | 80 | 238 |
| T3 CMP→GMP | 6 | 441 | 367 | 802 |
| T4 GMP→MB (추론) | – | 35 | 194 | – |
| T5 MB→PM | 2 | 448 | 460 | 908 |
| T6 PM→MY | 4 | 546 | 688 | 1,208 |
| T7 MY→MM | 3 | 579 | 1,092 | 1,663 |
| T8 MM→BN | 3 | 500 | 1,050 | 1,546 |
| T9 BN→SN | 3 | 555 | 895 | 1,442 |
| B3 PM→MY/MM (블록) | 7 | 1,764 | 2,442 | 4,189 |
| B4 MY/MM→BN/SN (블록) | 6 | 1,583 | 2,787 | 4,346 |

- 데이터셋 하나씩 빼기: HC 호출의 95–100%가 유지(B2 95%, 나머지 98–100%).
- BLUEPRINT 단계 표기를 저자 기준으로 바꿔도 HC 집합이 거의 같다(B3 4,206 vs 4,147, 공통 4,143).
- 합의 곡선 검증: 개별 데이터셋 단계 평균과의 Spearman ρ 중앙값 — 기능 패널 0.91(IQR 0.87–0.95, 유전자당 데이터셋 중앙 10개), HC 유전자 전체 9,315개 0.87(88%가 ρ > 0.5).

### 4.2 기능 패널과 근거 등급

- 등급별 유전자(기능 패널 69개):
  - **A (38)**: CD34 KIT MKI67 TOP2A MYB CEBPA GFI1 CEBPE CEBPB MPO ELANE AZU1 PRTN3 CTSG DEFA4 DEFA1/3 BPI LTF CAMP LCN2 CEACAM8 OLFM4 MMP9 ARG1 CR1 CYBB CD33 ITGA4 ITGAM FCGR3B MME CD101 CXCR4 CXCR2 CXCR1 FPR1 SELL CSF3R
  - **B (29)**: PROM1 HLF MECOM CRHBP MLLT3 FLT3 CCNB1 CDK1 PCNA MCM2 AURKB BUB1 GATA2 MEIS1 HOXA9 MYC SPI1 IRF8 KLF4 CEBPD MMP8 CHI3L1 FCN1 ALPL CYBA NCF1 NCF2 NCF4 C5AR1
  - **C (2)**: THY1(CD90; HSC에만, 데이터셋 1개), FCGR1A(CD64; 단계 변화 없음)
- 그림 11a의 모듈은 같은 패턴의 유전자만 묶었다(구성 유전자 간 상관 r 0.65–1.00). 상관이 낮던 묶음(HOXA9·MEIS1·GATA2 0.46, 디펜신·BPI 0.43, 피콜린·분비소포 0.41)은 나눴다.
- 진술 판정 요약(사람): mRNA 지지 148 · 불명확 23 · 반대 3 · 판정 불가 52(혈액 PMN) · 메타 자료 없음 2. 표면 단백질 지지 130 · 불명확 30 · 반대 11(대부분 CXCR4처럼 mRNA와 표면 발현이 원래 다른 경우). 분획 정의 자체를 서술한 169건과 질환·자극 맥락 연구 56건은 근거에서 뺐다.
- 단계 정점 마커(보조; 이웃 단계 최대 h < 0.7 + 들어오는 증가 또는 나가는 감소가 HC): PM CTSG* · MY 세포주기 유전자(PTTG1·AURKB·KIF23 — 마지막 증식 단계) · MM CD101*·LTF*·CAMP* · BN MMP9*·HMGB2* · SN MME*·FCGR3A*·CXCR2*·NAMPT* (*: 등급 A). 엑셀 `Stage_peak_markers` 시트.

### 4.3 단계별 증가·감소 유전자 (상위, 단백질 암호화)

상위 10개(증가 ▲ / 감소 ▼, *는 등급 A, 전체는 엑셀 `Stage_transition_panel`):

| 이행 | 증가 | 감소 |
|---|---|---|
| HSC→MPP | TFRC* · KLF1 · GATA1 · CSF2RB (적혈구 준비) | EMCN · JAM2 · GIMAP7 |
| MPP→CMP | MKI67* · TOP2A* · CXCR4* · PLAUR* · HDC · MS4A2 | HLF · CRHBP · HOPX · MN1 |
| CMP→GMP | MPO 계열(PRTN3* · AZU1* · CTSG*) · LYZ* · S100A8/9* · GFI1* · CD33* · CYBB* | TFRC* · SOD2* · HDC · MS4A2 |
| GMP→MB (추론) | ELANE* · PRTN3* · PADI4* · CEBPB* · SERPINB1 | SOCS2 · MN1 · JAM2 · CD52 |
| MB→PM | BPI* · DEFA4* · CEACAM8* · OLFM4* · ARG1* · CD24* · S100A8/9* | CD34* · SELL* · SPINK2 · EGFL7 |
| PM→MY | ITGAM* · CYBB* · CD101* · CAMP* · LCN2* · CD177* · PADI4* | MPO* · ELANE* · PRTN3* · CTSG* · KIT* · CDK4/6* |
| MY→MM | MMP9* · FPR1* · CXCR1* · SELL* · G0S2* · S100A12* | AZU1* · MPO* · ELANE* · ITGA4* · DEFA4* · CDK4/6* |
| MM→BN | MME* · CR1* · CSF3R* · FCGR3B* · NAMPT* · PLAUR* | MKI67* · TOP2A* · MYB* · GFI1* · CEBPE* · CEACAM8* |
| BN→SN | MME* · CXCR2* · CXCR4* · CSF3R* · NAMPT* · DDIT3* | ANXA1* · CEBPE* · MYB* · BPI* · CD24* · ARG1* |

- HSC→CMP 구간의 적혈구·비만세포 계열 유전자는 계통 준비 신호로 그대로 두었고, GMP→SN 구간에서는 다른 계통 표지를 혼입 신호로 보고 목록에서 뺐다.

## 5. 해석상 주의

1. **문(gate) 이름의 오프셋**: 같은 CD11b/CD16 분획을 연구마다 약 한 단계 다르게 부른다(원고 Box 2). 단계 경계는 ±1 단계 불확실성을 가진다.
2. **정렬 분획 혼입**: GSE294330 PM 분획에는 성숙 호중구 유전자(FCGR3B·FPR1·MME·CXCR2)가 MY보다 높게 나타나 약 1% 수준의 성숙세포 혼입이 있다. GSE42519의 HSC·MPP·CMP 분획에는 적혈구·과립구 RNA가 섞여 있다. 데이터셋 간 일치 조건이 이런 단일 데이터셋 신호를 거르지만, 일치 조건을 통과하지 못한 실제 변화도 있을 수 있다.
3. **MB 근거가 얇다**: MB→PM(T5)은 2개 데이터셋뿐이고, GMP→MB는 추론이다.
4. **mRNA ≠ 단백질·표면 발현**: 과립 단백질은 전사가 끝난 뒤에도 저장되고, CXCR4는 mRNA가 성숙할수록 오르지만 표면 발현은 떨어진다. 패널은 mRNA 기준이다.
5. **HSC–GMP 구간의 문헌 근거**: 원고의 검증 참고문헌이 호중구 중심이라, HSC/MPP 유전자는 데이터 메타분석(3–6개 데이터셋)만으로 등급 B가 많다.
6. 마이크로어레이 포화·압축, BLUEPRINT의 단백질 검출 유전자 제한(3,455개), 단일세포 자료 미포함(이번 정량 합성은 정렬 벌크 분획만).

## 6. 파일

| 구분 | 경로 |
|---|---|
| 검색·선별 | `02_scripts/06_neutrophil_meta_search.py` · `04_results/06_meta_search/` · `03_refs/08_neutrophil_meta_study_registry.tsv` |
| 자료 수집 | `02_scripts/07_neutrophil_meta_datasets_fetch.py` · `00_data/geo/` · `04_results/07_meta_datasets/07_download_log.tsv` |
| 표준화 | `02_scripts/08_neutrophil_meta_datasets_build.R` · `04_results/08_meta_datasets/` |
| 메타분석 | `02_scripts/09_neutrophil_meta_analysis.R` · `04_results/09_meta_analysis/` |
| 문헌 대조·등급·패널 | `02_scripts/10_neutrophil_meta_panel.R` · `02_scripts/10a_neutrophil_literature_claims_validate.py` · `03_refs/07_neutrophil_literature_stage_claims.tsv` · `04_results/10_meta_panel/` |
| 그림 | `02_scripts/11_neutrophil_stage_panel_figures.py` · `04_results/11_figures/` |
| 통합 엑셀 | `02_scripts/12_neutrophil_stage_panel_workbook.py` · `04_results/12_neutrophil_stage_mrna_panel.xlsx` |
| 초기 3-데이터셋 분석(대체됨) | `02_scripts/04_*`, `05_*` · `04_results/04_stage_expression/`, `05_stage_transition_de/` |

재현: 프로젝트 루트에서 R 스크립트는 `run-r bi-bulk-enrichment <script>`, Python은 `/Users/edemhw/miniforge3/bin/python <script>`로 08 → 09 → 10 → 11 → 12 순서로 실행(06·07은 검색·다운로드).
