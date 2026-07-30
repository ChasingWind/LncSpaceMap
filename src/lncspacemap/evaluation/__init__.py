"""Evaluation utilities."""

from .minimal import evaluate_predictions
from .multifold import aggregate_multifold

__all__ = ["aggregate_multifold", "evaluate_predictions"]
