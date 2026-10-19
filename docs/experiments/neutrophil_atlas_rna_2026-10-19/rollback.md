# 실험 분리와 되돌리기

- 실험 branch: `codex/neutrophil-atlas-rna-2026-10-19`. 기본 branch에 merge하지 않는다.
- 기준선 commit: `ee219c5` (기존 NSO source/test21개를 정확히 복제).
- 설계/다운로드 목록 commit: `1fdff9e`. 이 시점은 새 RNA 패널 구현 전이다.
- 기존 작업은 `/Users/edemhw/Honglab/CS/outputs/dcta_neighbor_state_route_2026-10-16`에 그대로 있다. 이 사본을 계속 쓰면 기존 패널/결과로 돌아간다.
- 실험은 별도 `/Users/edemhw/Honglab/CS/outputs/dcta_neutrophil_atlas_rna_2026-10-19`에서만 진행한다. 공유 폴더의 branch를 바꾸지 않는다.
- 구현 후에는 별도 feature commit을 남긴다. 취소할 때는 해당 feature commit만 revert하며, 현재 미커밋 파일에 reset/clean/stash를 사용하지 않는다.
- 원자료는 `data/external/`에서 Git 제외. 다시 받는 URL/SHA256은 `source_downloads.json`에 있다.
- 기존 source21개와 observation/connectivity manifest2개가 바뀌지 않았음을 Codex가 확인했다.
