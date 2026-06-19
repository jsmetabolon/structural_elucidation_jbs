"""HTML compound exploration reports for X-compound structural elucidation.

The ``compound`` module builds self-contained public/internal HTML dashboards
from rerank outputs, metabolite DB lookups, and similarity networks.
"""

from __future__ import annotations

from structural_elucidation.reports.compound import build_report, main

__all__ = ["build_report", "main"]
