#!/usr/bin/env python3
"""Lightweight wall-clock stage timing for pipeline CLI modules.

Emits one ``[timing]`` JSON line per stage on stderr so orchestration scripts
(and humans reading logs) can attribute wall-clock to pipeline stages without
a profiler. Timing must never alter pipeline behavior: the context manager
re-raises everything and adds no work on the hot path.
"""

from __future__ import annotations

import json
import sys
import time
from contextlib import contextmanager
from typing import Any, Iterator


@contextmanager
def stage_timer(stage: str, **extra: Any) -> Iterator[None]:
    """Print a ``[timing]`` JSON line with elapsed seconds for *stage*."""
    start = time.perf_counter()
    try:
        yield
    finally:
        payload = {
            "event": "stage_timing",
            "stage": str(stage),
            "seconds": round(time.perf_counter() - start, 3),
        }
        for key, value in extra.items():
            payload[str(key)] = value if isinstance(value, (str, int, float, bool)) else str(value)
        print(f"[timing] {json.dumps(payload, sort_keys=True)}", file=sys.stderr, flush=True)
