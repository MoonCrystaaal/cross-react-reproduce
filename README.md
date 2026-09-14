# Cross-React 재현 및 순위 평가

이 저장소는 알레르겐의 표면 패치가 얼마나 비슷한지 살펴보고, 이를 바탕으로 교차 반응성을 예측하는 Cross-React 도구의 계산 흐름을 재현한 코드입니다. 표면 패치를 20×20 connectivity matrix로 나타낸 뒤 PCC 점수를 계산해 서로 비교합니다.


## 환경과 실행

Linux 또는 WSL 환경에서 실행하는 것을 기준으로 정리했습니다. Python 3.11, NumPy, FreeSASA Python 모듈과 CLI가 필요합니다.

```bash
conda env create -f pipeline/environment.yml
conda activate cross-react-pipeline
```

전체 명령과 입력 CSV 형식은 [실행 안내](pipeline/README.md)에 정리했습니다. 저장소 루트 아래 `library/`를 작업 폴더로 사용하며, PDB 구조 파일은 별도로 준비합니다. 기존 로컬 구조 파일은 `pipeline/pdb_all/` 등에 있으므로 실행할 때 `--pdb-dir`와 `--input-dir`를 실제 위치에 맞춥니다.

## 코드와 방법

| 코드 | 역할 |
|---|---|
| `pipeline/1_extract_surface_residues_freesasa.py` | SASA > 10 Å²인 표면 잔기 선정 |
| `pipeline/2_build_patches_10A.py` | 대표 원자(Cβ, Gly는 Cα)를 중심으로 10 Å patch 구성 |
| `pipeline/3_build_connectivity_matrices.py` | 거리 ≤ 8 Å 접촉으로 대칭 20×20 행렬 생성 |
| `pipeline/4_build_query_matrix_and_score_pcc.py` | 지정 epitope와 후보 patch 비교 |
| `pipeline/5_score_allergen_patch_pairs_legacy.py` | query와 target의 모든 patch 쌍 비교, 최고 점수로 후보 순위 산출 |
| `pipeline/legacy_pcc.py` | off-diagonal 합산과 분자 절댓값을 포함하는 legacy 식 |
| `pipeline/rerank_allergen_rankings.py` | 기존 순위 파일 재정리 |
| `pipeline/benchmark_ranking_metrics.py` | 같은 점수 그룹의 순열에 대한 Hits/Recall/NDCG 기댓값 평가 |

여기서 사용하는 PCC는 일반적인 400개 원소 Pearson 상관계수와 계산 방식이 다릅니다. Stage 5의 hit 저장 임계값은 전체 후보 최고 점수 계산을 자르는 기준이 아니라, 중간 hit 저장 범위를 정하는 값입니다. 원 웹 도구 재현용 epitope–patch 검증과 patch–patch 후보 평가는 별도로 다룹니다.

```bash
python -m unittest discover -s pipeline -p 'test_*.py' -v
python pipeline/benchmark_ranking_metrics.py --help
```

## 데이터와 제출

`library/pair_label_dataset.xlsx`와 query 목록은 재현 입력으로 유지합니다. 생성 행렬, patch hit, 순위, 캐시, 백업, PDB 묶음은 Git에서 제외하고 로컬 파일로 보관합니다. 보고서 증빙에 필요한 소형 요약은 `submission/` 아래에 따로 정리했습니다.

평가 결과를 다시 확인할 때는 라벨 SHA-256, query 목록, 후보 집합, 점수 파일 및 동점 처리 기준을 함께 확인합니다.

원격 저장소는 아직 설정되지 않았습니다. 제출 링크에는 본인 저장소를 만든 뒤 실제 URL과 제출 커밋 해시를 기입합니다.

## 최종 보고서 평가 자료

[최종 평가 증빙](submission/benchmark/README.md)에 갱신 라벨(181개·양성 133쌍), 두 방법의 동일한 30 query 순위, 기존 집계 및 재집계 검증 코드를 포함했습니다. 본 저장소의 평가 코드는 최종 집계에 사용한 코드와 SHA-256이 같습니다. [재현 사례 기록](submission/cit_s_7/README.md)에는 PCC 0.85 사례의 저장 결과를 포함했습니다.
