#!/usr/bin/env python3
from __future__ import annotations

import torch

DEVICE_HELP_TEXT = "Device to run on: auto, cpu, cuda, or mps."


def _mps_available() -> bool:
    mps_backend = getattr(torch.backends, "mps", None)
    return bool(mps_backend is not None and mps_backend.is_available())


def resolve_device(raw_device: str) -> torch.device:
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
