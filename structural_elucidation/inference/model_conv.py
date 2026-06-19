#!/usr/bin/env python3
"""1-D CNN spectrum encoder (``Net1D``) and checkpoint compatibility helpers.

Supports both legacy config-driven constructors and the mass-conditioned
encoder used by X209 / merged-library inference. Shape inference utilities
recover compatible ``length_in`` from stored ``fc1`` / ``conv3`` weights.
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


def conv_out_dim(
    length_in: int,
    kernel: int,
    stride: int,
    padding: int,
    dilation: int,
) -> int:
    """Output length after a 1-D convolution or pooling layer."""
    return (length_in + 2 * padding - dilation * (kernel - 1) - 1) // stride + 1


def _pool_halve_dim(length_in: int) -> int:
    """Length after ``MaxPool1d(kernel=2, stride=2)``."""
    return conv_out_dim(
        length_in=length_in,
        kernel=2,
        stride=2,
        padding=0,
        dilation=1,
    )


def pooled_length_after_three_halves(length_in: int) -> int:
    """Spectral length after three successive halving pools in ``Net1D``."""
    length = int(length_in)
    for _ in range(3):
        length = _pool_halve_dim(length)
    return int(length)


def infer_mass_conditioned_compatible_length(
    state_dict: dict[str, torch.Tensor],
    preferred_length_in: int | None = None,
) -> int | None:
    """Infer ``length_in`` compatible with a mass-conditioned checkpoint.

    Derives pooled length from ``fc1.weight`` in-features and ``conv3`` channel
    count, then maps back to an 8×-upsampled input width band.

    Args:
        state_dict: Encoder checkpoint tensors.
        preferred_length_in: Caller-requested bin count; clamped into the
            compatible band when present.

    Returns:
        Compatible ``length_in``, or ``preferred_length_in`` when tensors are
        missing or shapes are inconsistent.
    """
    fc1_weight = state_dict.get("fc1.weight")
    conv3_weight = state_dict.get("conv3.weight")
    if not isinstance(fc1_weight, torch.Tensor):
        return preferred_length_in
    if not isinstance(conv3_weight, torch.Tensor):
        return preferred_length_in
    if fc1_weight.ndim != 2 or conv3_weight.ndim != 3:
        return preferred_length_in

    channels_out = int(conv3_weight.shape[0])
    fc1_in_features = int(fc1_weight.shape[1])
    if channels_out <= 0 or fc1_in_features <= 0:
        return preferred_length_in
    if fc1_in_features % channels_out != 0:
        return preferred_length_in

    pooled_length = fc1_in_features // channels_out
    min_length = int(pooled_length * 8)
    max_length = int(min_length + 7)
    if preferred_length_in is None:
        return max_length
    preferred = int(preferred_length_in)
    if preferred < min_length:
        return min_length
    if preferred > max_length:
        return max_length
    return preferred


def load_compatible_state_dict(
    model: nn.Module,
    state_dict: dict[str, torch.Tensor],
) -> dict[str, object]:
    """Load only shape-matching tensors; report skips and missing keys.

    Args:
        model: Target module (typically ``Net1D``).
        state_dict: Checkpoint weights.

    Returns:
        Load report with counts, key lists, and per-key shape mismatches.
    """
    model_state = model.state_dict()
    filtered: dict[str, torch.Tensor] = {}
    skipped_mismatched: dict[str, dict[str, list[int]]] = {}
    unexpected_keys: list[str] = []

    for key, value in state_dict.items():
        if key not in model_state:
            unexpected_keys.append(key)
            continue
        target = model_state[key]
        if not isinstance(value, torch.Tensor):
            unexpected_keys.append(key)
            continue
        if tuple(value.shape) != tuple(target.shape):
            skipped_mismatched[key] = {
                "checkpoint": list(value.shape),
                "model": list(target.shape),
            }
            continue
        filtered[key] = value

    incompatible = model.load_state_dict(filtered, strict=False)
    return {
        "loaded_key_count": int(len(filtered)),
        "loaded_keys": sorted(filtered.keys()),
        "missing_keys": list(incompatible.missing_keys),
        "unexpected_keys": sorted(set(unexpected_keys) | set(incompatible.unexpected_keys)),
        "skipped_mismatched": skipped_mismatched,
    }


class Net1D(nn.Module):
    """Encoder used by legacy and mass-conditioned spectral workflows.

    Supported constructor patterns:

    - **Mass-conditioned (current):**
      ``Net1D(length_in=..., minmass=..., maxmass=..., ...)``
    - **Legacy:** ``Net1D(args, length_in)`` with a config object as first arg.
    """

    def __init__(self, *args, **kwargs):
        super().__init__()
        self.legacy_mode = False

        if args and not isinstance(args[0], (int, float)):
            self._init_legacy(*args, **kwargs)
            self.legacy_mode = True
            return

        self._init_mass_conditioned(*args, **kwargs)

    def _init_legacy(self, *args, **kwargs) -> None:
        """Build encoder from historical config objects (pre-mass-conditioning)."""
        cfg = args[0] if args else kwargs.pop("args")
        if len(args) >= 2:
            length_in = int(args[1])
        else:
            length_in = int(kwargs.pop("length_in"))
        self.length_in = int(length_in)

        out_1 = conv_out_dim(
            length_in,
            int(cfg.conv_kernel_dim_1),
            int(cfg.conv_stride_1),
            int(cfg.conv_padding_1),
            int(cfg.conv_dilation),
        )
        out_2 = conv_out_dim(
            out_1,
            int(cfg.pool_kernel_dim_1),
            int(cfg.pool_stride_1),
            int(cfg.pool_padding_1),
            int(cfg.pool_dilation),
        )
        out_3 = conv_out_dim(
            out_2,
            int(cfg.conv_kernel_dim_2),
            int(cfg.conv_stride_2),
            int(cfg.conv_padding_2),
            int(cfg.conv_dilation),
        )
        out_4 = conv_out_dim(
            out_3,
            int(cfg.pool_kernel_dim_2),
            int(cfg.pool_stride_2),
            int(cfg.pool_padding_2),
            int(cfg.pool_dilation),
        )
        self.cnn_out = int(cfg.channels_out) * out_4
        self.conv1 = nn.Conv1d(
            int(cfg.channels_in),
            int(cfg.channels_med_1),
            int(cfg.conv_kernel_dim_1),
            stride=int(cfg.conv_stride_1),
            padding=int(cfg.conv_padding_1),
            dilation=int(cfg.conv_dilation),
        )
        self.norm1 = nn.BatchNorm1d(int(cfg.channels_med_1))
        self.pool1 = nn.MaxPool1d(
            int(cfg.pool_kernel_dim_1),
            stride=int(cfg.pool_stride_1),
            padding=int(cfg.pool_padding_1),
            dilation=int(cfg.pool_dilation),
        )
        self.conv2 = nn.Conv1d(
            int(cfg.channels_med_1),
            int(cfg.channels_out),
            int(cfg.conv_kernel_dim_2),
            stride=int(cfg.conv_stride_2),
            padding=int(cfg.conv_padding_2),
            dilation=int(cfg.conv_dilation),
        )
        self.norm2 = nn.BatchNorm1d(int(cfg.channels_out))
        self.pool2 = nn.MaxPool1d(
            int(cfg.pool_kernel_dim_2),
            stride=int(cfg.pool_stride_2),
            padding=int(cfg.pool_padding_2),
            dilation=int(cfg.pool_dilation),
        )
        self.fc1 = nn.Linear(int(cfg.channels_out) * out_4, int(cfg.fc_dim_1))
        self.fc2 = nn.Linear(int(cfg.fc_dim_1), int(cfg.emb_dim))
        self.norm3 = nn.BatchNorm1d(int(cfg.fc_dim_1))

        # Legacy scripts expect these attributes even without mass conditioning.
        self.apply_tanh = True
        self.minmass = 0.0
        self.maxmass = 1.0
        self.parent_mz_min = 0.0
        self.parent_mz_max = 1.0

    def _init_mass_conditioned(self, *args, **kwargs) -> None:
        """Build three-conv mass-conditioned encoder matching augmented checkpoints."""
        if "length_in" in kwargs:
            length_in = int(kwargs.pop("length_in"))
        elif args:
            length_in = int(args[0])
            args = args[1:]
        else:
            raise TypeError("Net1D requires `length_in` " "for mass-conditioned mode.")
        if args:
            raise TypeError(f"Unexpected positional arguments: {args}")
        self.length_in = int(length_in)

        self.minmass = float(kwargs.pop("minmass", 50.0))
        self.maxmass = float(kwargs.pop("maxmass", 499.97))
        self.parent_mz_min = float(kwargs.pop("parent_mz_min", self.minmass))
        self.parent_mz_max = float(kwargs.pop("parent_mz_max", self.maxmass))
        self.apply_tanh = bool(kwargs.pop("apply_tanh", True))
        if kwargs:
            raise TypeError("Unexpected keyword arguments: " f"{sorted(kwargs.keys())}")

        # Shape-compatible with models/spectra_encoder_augmented_209_400/best_encoder.pt
        self.conv1 = nn.Conv1d(1, 32, kernel_size=5, stride=1, padding=2)
        self.bn1 = nn.BatchNorm1d(32)
        self.conv2 = nn.Conv1d(32, 64, kernel_size=5, stride=1, padding=2)
        self.bn2 = nn.BatchNorm1d(64)
        self.conv3 = nn.Conv1d(64, 128, kernel_size=5, stride=1, padding=2)
        self.bn3 = nn.BatchNorm1d(128)
        self.pool = nn.MaxPool1d(kernel_size=2, stride=2)

        out_1 = _pool_halve_dim(length_in)
        out_2 = _pool_halve_dim(out_1)
        out_3 = _pool_halve_dim(out_2)
        if out_3 <= 0:
            raise ValueError("Net1D pooled length is non-positive. " f"Got length_in={length_in}, pooled={out_3}.")
        self.pooled_length = int(out_3)
        self.cnn_out = 128 * out_3

        self.fc1 = nn.Linear(self.cnn_out, 512)
        self.bn4 = nn.BatchNorm1d(512)
        self.fc2 = nn.Linear(512, 512)

        self.mz_mlp = nn.Sequential(
            nn.Linear(1, 64),
            nn.ReLU(),
            nn.Linear(64, 512),
        )
        self.mz_alpha = nn.Parameter(torch.tensor(1.0, dtype=torch.float32))

    def _normalize_parent_mz(self, parent_mz: torch.Tensor) -> torch.Tensor:
        """Scale precursor m/z to [0, 1] using training normalization bounds."""
        denom = max(self.parent_mz_max - self.parent_mz_min, 1e-6)
        normalized = (parent_mz - self.parent_mz_min) / denom
        return normalized.clamp(0.0, 1.0).unsqueeze(-1)

    def forward(
        self,
        x: torch.Tensor,
        parent_mz: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Encode 1-D spectra into a latent embedding with optional parent m/z conditioning."""
        if self.legacy_mode:
            x = self.pool1(F.relu(self.norm1(self.conv1(x))))
            x = self.pool2(F.relu(self.norm2(self.conv2(x))))
            x = x.view(-1, self.cnn_out)
            x = F.relu(self.norm3(self.fc1(x)))
            x = self.fc2(x)
            return torch.tanh(x) if self.apply_tanh else x

        x = self.pool(F.relu(self.bn1(self.conv1(x))))
        x = self.pool(F.relu(self.bn2(self.conv2(x))))
        x = self.pool(F.relu(self.bn3(self.conv3(x))))
        x = x.view(x.size(0), -1)
        x = F.relu(self.bn4(self.fc1(x)))
        x = self.fc2(x)

        if parent_mz is not None:
            parent_mz = parent_mz.to(dtype=x.dtype, device=x.device)
            if parent_mz.ndim == 0:
                parent_mz = parent_mz.unsqueeze(0)
            mz_feat = self.mz_mlp(self._normalize_parent_mz(parent_mz))
            x = x + self.mz_alpha * mz_feat

        return torch.tanh(x) if self.apply_tanh else x
