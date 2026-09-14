# 보고서 최종 평가 증빙

갱신 라벨로 재평가한 30 query 비교 자료입니다. 새로 검색을 돌린 자료가 아니라, 저장해 둔 순위 파일을 같은 기준으로 다시 집계한 결과입니다.

`labels_181.xlsx`는 최종 181개 알레르겐 라벨이며 `library/pair_label_dataset.xlsx`의 기존 182개 입력과 구분합니다. `crossreact_rankings/`와 `surfaceid_rankings/`는 같은 30 query와 query별 180 후보로 맞춘 소형 순위 자료입니다. 원 점수 자료에서 Jug r 4를 제외한 집합입니다.

`comparison_summary.csv`, `expected_*/`는 저장 결과이며, `evaluation_configuration.json`과 `input_audit.json`은 입력 파일과 해시를 확인하기 위한 기록입니다. 아래 명령으로 저장소에 포함된 입력을 다시 검증할 수 있습니다.

```bash
python submission/benchmark/verify_and_recompute.py --output-dir outputs/submission-verification
```

이 명령은 파일 해시, query·후보 동일성, 양성 후보 누락 여부, 평균 및 query별 지표 일치를 확인합니다. 모델 검색, 구조 전처리, 임상적 타당성, 두 방법의 입력 구조 동일성까지 검증하는 절차는 아닙니다.

보고서의 K=5 값은 Cross-React 0.88/0.74/0.73, SurfaceID 0.73/0.58/0.62입니다. 각 값은 Hits/Recall/NDCG 순서이며 반올림 전 수치는 CSV에 보존했습니다.
