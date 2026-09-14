#!/usr/bin/env python3
"""Shared implementation of the legacy Cross-React correlation formula."""

from __future__ import annotations

import math
from typing import Sequence, Tuple

import numpy as np


MATRIX_SIZE = 20
MATRIX_ELEMENT_COUNT = MATRIX_SIZE * MATRIX_SIZE


def _as_square_matrix(raw_matrix: Sequence[Sequence[float]]) -> np.ndarray:
    matrix = np.asarray(raw_matrix, dtype=np.float64)
    if matrix.shape != (MATRIX_SIZE, MATRIX_SIZE):
        raise ValueError(
            f"Expected a {MATRIX_SIZE}x{MATRIX_SIZE} matrix, got {matrix.shape}"
        )
    return matrix


def build_legacy_row_sums(
    raw_matrix: Sequence[Sequence[float]],
) -> Tuple[float, ...]:
    """Convert one raw matrix to the 20 legacy symmetric-matrix row sums."""
    matrix = _as_square_matrix(raw_matrix)
    symmetric = matrix + matrix.T
    row_sums = symmetric.sum(axis=1) - np.diag(matrix)
    return tuple(float(value) for value in row_sums)


def build_legacy_stats(
    raw_matrix: Sequence[Sequence[float]],
) -> Tuple[Tuple[float, ...], float, float]:
    row_sums = build_legacy_row_sums(raw_matrix)
    matrix_sum = float(sum(row_sums))
    matrix_sumsq = float(sum(value * value for value in row_sums))
    return row_sums, matrix_sum, matrix_sumsq


def legacy_crossreact_score_from_stats(
    query_row_sums: Sequence[float],
    query_sum: float,
    query_sumsq: float,
    target_row_sums: Sequence[float],
    target_sum: float,
    target_sumsq: float,
) -> float:
    """Calculate one score exactly as the legacy 20x20 implementation does."""
    if len(query_row_sums) != MATRIX_SIZE or len(target_row_sums) != MATRIX_SIZE:
        raise ValueError(f"Legacy row-sum vectors must contain {MATRIX_SIZE} values")

    sum12 = float(
        sum(left * right for left, right in zip(query_row_sums, target_row_sums))
    )
    meanxy = (float(query_sum) * float(target_sum)) / MATRIX_ELEMENT_COUNT
    sumsquare1 = (float(query_sum) * float(query_sum)) / MATRIX_ELEMENT_COUNT
    sumsquare2 = (float(target_sum) * float(target_sum)) / MATRIX_ELEMENT_COUNT

    upper1 = sum12 - meanxy
    denom_left = float(query_sumsq) - sumsquare1
    denom_right = float(target_sumsq) - sumsquare2
    if denom_left <= 0.0 or denom_right <= 0.0:
        return float("nan")

    lower1 = math.sqrt(denom_left * denom_right)
    if lower1 == 0.0:
        return float("nan")
    return abs(upper1) / lower1


def legacy_crossreact_score(
    query_matrix: Sequence[Sequence[float]],
    target_matrix: Sequence[Sequence[float]],
) -> float:
    query_row_sums, query_sum, query_sumsq = build_legacy_stats(query_matrix)
    target_row_sums, target_sum, target_sumsq = build_legacy_stats(target_matrix)
    return legacy_crossreact_score_from_stats(
        query_row_sums=query_row_sums,
        query_sum=query_sum,
        query_sumsq=query_sumsq,
        target_row_sums=target_row_sums,
        target_sum=target_sum,
        target_sumsq=target_sumsq,
    )


def legacy_crossreact_scores_from_stats(
    query_row_sums: np.ndarray,
    query_sums: np.ndarray,
    query_sumsq: np.ndarray,
    target_row_sums: np.ndarray,
    target_sums: np.ndarray,
    target_sumsq: np.ndarray,
) -> np.ndarray:
    """Calculate a query-patch by target-patch score matrix in one NumPy call."""
    query_vectors = np.asarray(query_row_sums, dtype=np.float64)
    target_vectors = np.asarray(target_row_sums, dtype=np.float64)
    query_sum_values = np.asarray(query_sums, dtype=np.float64)
    target_sum_values = np.asarray(target_sums, dtype=np.float64)
    query_sumsq_values = np.asarray(query_sumsq, dtype=np.float64)
    target_sumsq_values = np.asarray(target_sumsq, dtype=np.float64)

    if query_vectors.ndim != 2 or query_vectors.shape[1] != MATRIX_SIZE:
        raise ValueError(
            f"query_row_sums must have shape (n, {MATRIX_SIZE}), got {query_vectors.shape}"
        )
    if target_vectors.ndim != 2 or target_vectors.shape[1] != MATRIX_SIZE:
        raise ValueError(
            f"target_row_sums must have shape (n, {MATRIX_SIZE}), got {target_vectors.shape}"
        )
    if query_sum_values.shape != (query_vectors.shape[0],):
        raise ValueError("query_sums length does not match query_row_sums")
    if query_sumsq_values.shape != (query_vectors.shape[0],):
        raise ValueError("query_sumsq length does not match query_row_sums")
    if target_sum_values.shape != (target_vectors.shape[0],):
        raise ValueError("target_sums length does not match target_row_sums")
    if target_sumsq_values.shape != (target_vectors.shape[0],):
        raise ValueError("target_sumsq length does not match target_row_sums")

    sum12 = query_vectors @ target_vectors.T
    numerator = np.abs(
        sum12
        - np.multiply.outer(query_sum_values, target_sum_values)
        / MATRIX_ELEMENT_COUNT
    )
    query_denominator = (
        query_sumsq_values
        - np.square(query_sum_values) / MATRIX_ELEMENT_COUNT
    )
    target_denominator = (
        target_sumsq_values
        - np.square(target_sum_values) / MATRIX_ELEMENT_COUNT
    )
    valid = np.logical_and.outer(query_denominator > 0.0, target_denominator > 0.0)
    denominator = np.sqrt(
        np.maximum(query_denominator[:, None], 0.0)
        * np.maximum(target_denominator[None, :], 0.0)
    )
    scores = np.full(numerator.shape, np.nan, dtype=np.float64)
    np.divide(numerator, denominator, out=scores, where=valid & (denominator > 0.0))
    return scores
