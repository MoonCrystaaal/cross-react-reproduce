# Cross-React Reproduction Pipeline

## Stage 4와 Stage 5

- Stage 4는 기존 방식대로 `query PDB + epitope residue CSV`로 query matrix 하나를
  만든 뒤 target patch들과 비교합니다.
- Stage 5는 epitope 입력 없이 `query_allergens.txt`에서 선택한 allergen의 모든
  patch를 모든 target allergen의 patch와 비교합니다.
- 두 단계는 모두 `legacy_pcc.py`에 있는 동일한 legacy Cross-React 공식을
  사용합니다.

전체 patch-patch 방식은 작업 디렉터리(예: `library`)에서 다음과 같이 실행합니다.

```powershell
python ..\pipeline\5_score_allergen_patch_pairs_legacy.py `
  --flat-dir ".\3_connectivity\connectivity_matrix_flat" `
  --output-root ".\5_allergen_pcc" `
  --query-list ".\query_allergens.txt" `
  --patch-hit-threshold 0.5 `
  --exclude-self `
  --block-size 256 `
  --workers 4
```

`query_allergens.txt`에는 connectivity matrix의 `allergen` 열과 정확히 같은 이름을
한 줄에 하나씩 기록합니다. 빈 줄과 `#` 주석은 무시하며, 중복 이름은 경고 후 한
번만 사용합니다. Excel label 파일을 자동으로 읽지는 않습니다.

`--patch-hit-threshold 0.5`는 legacy PCC가 0.5 이상인 모든 patch-patch 쌍을
`5_allergen_pcc/patch_hits/*.csv.gz`에 저장합니다. 이 threshold는 hit 저장에만
적용됩니다. Target allergen별 최고 PCC와 순위는 threshold 미만을 포함한 전체
유효 점수로 계산됩니다.

주요 출력은 다음과 같습니다.

- `5_allergen_pcc/patch_hits/`: query별 threshold 이상 patch hit
- `5_allergen_pcc/pair_scores/allergen_pair_best_scores.csv`: query-target별 최고 PCC
- `5_allergen_pcc/rankings/allergen_rankings.csv`: query별 allergen 순위
- `5_allergen_pcc/summary/allergen_scoring_summary.csv`: 계산 및 hit 개수 요약
- `5_allergen_pcc/checkpoints/`: `--resume` 재실행용 완료 기록

`--workers`는 동시에 처리할 query allergen 수입니다. 검증할 때는 작은 query
목록과 `--workers 1`을 사용하고, 전체 목록은 공식 및 출력 검증 후 실행하는 것을
권장합니다.

이 디렉터리는 Cross-React 재현용 전처리 3단계와 두 가지 scoring 모드를
제공합니다.

1. PDB별 surface residue 추출
2. surface residue 중심 10 A patch 생성
3. patch별 20x20 connectivity matrix 생성
4. epitope 기반 query matrix 생성 및 legacy Cross-React correlation 계산
5. 선택한 allergen의 전체 patch-patch legacy PCC 계산

## 설치

필요 환경:

- Linux
- Python 3.10 이상
- `freesasa` Python 모듈
- `freesasa` CLI 실행 파일

확인:

```bash
python -c "import freesasa; print(freesasa.__version__)"
freesasa --version
```

원하면 `environment.yml`을 사용할 수 있습니다.

## 실행

스크립트는 `pipeline/` 안에 두고, 작업용 폴더에서 실행합니다.
결과는 현재 작업 폴더 아래에 저장됩니다.

- `1_surface_residues/`
- `2_patches_10A/`
- `3_connectivity/`
- `4_pcc/`

필수 입력 폴더:

- `pdb_all/`
- `query_inputs/`

예시:

```bash
cd run_example
```

1단계:

```bash
python ../pipeline/1_extract_surface_residues_freesasa.py \
  --input-dir ./pdb_all \
  --output-dir ./1_surface_residues/surface_residue_csv \
  --rsa-dir ./1_surface_residues/surface_residue_rsa \
  --summary-dir ./1_surface_residues/surface_residue_summary \
  --log-dir ./1_surface_residues/surface_residue_logs \
  --algorithm LeeRichards \
  --probe-radius 1.4 \
  --surface-cutoff 10 \
  --n-slices 20 \
  --n-threads 1
```

2단계:

```bash
python ../pipeline/2_build_patches_10A.py \
  --pdb-dir ./pdb_all \
  --surface-csv-dir ./1_surface_residues/surface_residue_csv \
  --center-dir ./2_patches_10A/surface_center_coords \
  --members-dir ./2_patches_10A/patch_members \
  --summary-dir ./2_patches_10A/patch_summary \
  --log-dir ./2_patches_10A/patch_logs \
  --radius 10
```

3단계:

```bash
python ../pipeline/3_build_connectivity_matrices.py \
  --patch-members-dir ./2_patches_10A/patch_members \
  --residue-set-dir ./3_connectivity/patch_residue_sets \
  --matrix-dir ./3_connectivity/connectivity_matrices \
  --flat-dir ./3_connectivity/connectivity_matrix_flat \
  --summary-dir ./3_connectivity/connectivity_summary \
  --log-dir ./3_connectivity/connectivity_logs \
  --cutoff 8
```

4단계:

```bash
python ../pipeline/4_build_query_matrix_and_score_pcc.py \
  --query-pdb ./pdb_all/QUERY.pdb \
  --query-epitope-csv ./query_inputs/query_epitope.csv \
  --target-flat-dir ./3_connectivity/connectivity_matrix_flat \
  --query-matrix-dir ./4_pcc/query_matrix \
  --pcc-scores-dir ./4_pcc/pcc_scores \
  --pcc-ranked-dir ./4_pcc/pcc_ranked \
  --log-dir ./4_pcc/pcc_logs \
  --cutoff 8
```

## 입력 형식

query epitope CSV 최소 형식:

```csv
residue_number
42
43
44
```

같은 residue number가 여러 chain에 있으면:

```csv
chain_id,residue_number,residue_name
A,42,GLY
A,43,ARG
A,44,GLU
```

## 현재 구현 규칙

- 대표 원자: Gly 이외 `CB`, Gly 는 `CA`
- surface residue 기준: total SASA `> 10 A^2`
- patch 반경: `10 A`
- contact 기준: 대표 원자 거리 `<= 8 A`
- connected pair가 생기면 항상 `[i,j]`와 `[j,i]`를 모두 `+1`
- 같은 residue type contact는 결과적으로 diagonal이 `+2`
- stage 3 출력은 raw symmetric matrix
- stage 4 점수는 legacy Cross-React correlation 식 사용
- legacy 점수 계산 전, raw matrix의 off-diagonal은 내부 계산용으로 `m[i,j] + m[j,i]`로 변환
