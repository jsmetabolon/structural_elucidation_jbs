#!/usr/bin/env python3
"""PyTorch device resolution for inference CLIs.

Normalizes ``--device`` flags (``auto``, ``cpu``, ``cuda``, ``mps``) so
encoder and rerank scripts pick a consistent runtime without duplicating
availability checks in every entry point.
"""

from __future__ import annotations

import torch

# ── CLI help text ────────────────────────────────────────────────────
DEVICE_HELP_TEXT = "Device to run on: auto, cpu, cuda, or mps."


def _mps_available() -> bool:
    """Return whether Apple Metal (MPS) is available in the current PyTorch build."""
    mps_backend = getattr(torch.backends, "mps", None)
    return bool(mps_backend is not None and mps_backend.is_available())


def resolve_device(raw_device: str) -> torch.device:
    """Resolve a user device string to a validated ``torch.device``.

    ``auto`` prefers CUDA, then MPS, then CPU. Explicit ``cuda`` or ``mps``
    requests raise when the backend is unavailable so callers fail fast
    instead of silently falling back.

    Args:
        raw_device: One of ``auto``, ``cpu``, ``cuda``, ``mps``, or a
            torch-accepted device string.

    Returns:
        Resolved device ready for tensor allocation.

    Raises:
        ValueError: When ``raw_device`` is not a supported device name.
        RuntimeError: When CUDA or MPS was requested but is unavailable.
    """
    requested = str(raw_device).strip().lower()
    if requested == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if _mps_available():
            return torch.device("mps")
        return torch.device("cpu")

    try:
        device = torch.device(requested)
    except (TypeError, RuntimeError) as exc:
        raise ValueError(f"Unsupported device '{raw_device}'. " "Use one of: auto, cpu, cuda, mps.") from exc

    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("Requested CUDA device, but CUDA is not available.")
    if device.type == "mps" and not _mps_available():
        raise RuntimeError("Requested MPS device, but MPS is not available.")
    return device
