import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("benchmark", ROOT / "benchmark_ranking_metrics.py")
benchmark = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = benchmark
SPEC.loader.exec_module(benchmark)


class TieExpectedMetricsTest(unittest.TestCase):
    def test_top1_is_independent_of_name_order_with_two_way_tie(self):
        first = benchmark.compute_query_rows(
            "cross-react", "pcc", "q", (1, 1, 3),
            ("positive", "negative", "other"), (1.0, 1.0, 0.5),
            frozenset({"positive"}), (1,),
        )[0]
        second = benchmark.compute_query_rows(
            "cross-react", "pcc", "q", (1, 1, 3),
            ("negative", "positive", "other"), (1.0, 1.0, 0.5),
            frozenset({"positive"}), (1,),
        )[0]
        self.assertEqual(first, second)
        self.assertAlmostEqual(first["hit_at_k"], 0.5)
        self.assertAlmostEqual(first["recall_at_k"], 0.5)
        self.assertAlmostEqual(first["ndcg_at_k"], 0.5)
        self.assertAlmostEqual(first["weighted_hits_at_k"], 0.5)

    def test_complete_tie_group_has_expected_dcg(self):
        row = benchmark.compute_query_rows(
            "cross-react", "pcc", "q", (1, 1, 3),
            ("positive", "negative", "other"), (1.0, 1.0, 0.5),
            frozenset({"positive"}), (2,),
        )[0]
        expected_dcg = 0.5 * (1.0 + 1.0 / benchmark.math.log2(3))
        self.assertAlmostEqual(row["dcg_at_k"], expected_dcg)
        self.assertEqual(row["hit_at_k"], 1.0)
        self.assertEqual(row["recall_at_k"], 1.0)


if __name__ == "__main__":
    unittest.main()
