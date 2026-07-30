"""Mapping backend adapters."""

from .spage import run_spage
from .tangram import run_tangram

__all__ = ["run_spage", "run_tangram"]
