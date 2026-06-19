#!/usr/bin/env python3
"""Affine latent-space calibration applied after encoder inference.

Optional post-encoder transform that aligns query latents with the
embedding bank used for Track A retrieval. Checkpoints are identity-
initialized so missing calibrators are a no-op.
"""

from __future__ import annotations

from pathlib import Path

import torch
import torch.nn as nn


class AffineLatentCalibrator(nn.Module):
    """Learned affine map ``z -> W z + b`` over 512-d encoder latents."""

    def __init__(self, dim: int = 512) -> None:
        super().__init__()
        self.linear = nn.Linear(dim, dim)
        self.reset_parameters()

    def reset_parameters(self) -> None:
        """Initialize to identity so unloaded calibrators are neutral."""
        with torch.no_grad():
            self.linear.weight.copy_(torch.eye(self.linear.in_features, dtype=torch.float32))
            self.linear.bias.zero_()

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.linear(z)


def _extract_state_dict(payload: object) -> dict[str, torch.Tensor]:
    """Unwrap ``state_dict`` from a full training checkpoint or raw dict."""
    if isinstance(payload, dict) and "state_dict" in payload and isinstance(payload["state_dict"], dict):
        return payload["state_dict"]
    if isinstance(payload, dict):
        return payload
    raise TypeError(f"Unsupported calibrator payload type: {type(payload)}")


def load_latent_calibrator(path: Path, device: torch.device) -> AffineLatentCalibrator:
    """Load an affine calibrator checkpoint onto ``device`` in eval mode.

    Args:
        path: ``.pt`` checkpoint with linear weights (or nested ``state_dict``).
        device: Target inference device.

    Returns:
        Eval-mode ``AffineLatentCalibrator``.
    """
    payload = torch.load(path, map_location="cpu")
    state_dict = _extract_state_dict(payload)
    calibrator = AffineLatentCalibrator()
    calibrator.load_state_dict(state_dict, strict=False)
    calibrator.to(device)
    calibrator.eval()
    return calibrator


@torch.no_grad()
def apply_latent_calibrator(
    latent_tensor: torch.Tensor,
    calibrator: AffineLatentCalibrator | None,
) -> torch.Tensor:
    """Apply calibrator when present; otherwise return latents unchanged."""
    if calibrator is None:
        return latent_tensor
    return calibrator(latent_tensor)
