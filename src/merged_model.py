from __future__ import annotations

from typing import Any

import torch
import torch.nn as nn
import torch.nn.functional as F

from model1D2Conv import Net1D, load_compatible_state_dict


ION_MODE_TO_INDEX = {
    "negative": 0,
    "neg": 0,
    "n": 0,
    "0": 0,
    "0.0": 0,
    "positive": 1,
    "pos": 1,
    "p": 1,
    "1": 1,
    "1.0": 1,
}

LIB_ID_TO_INDEX = {
    209: 0,
    305: 1,
    400: 2,
    402: 3,
}
LIB_INDEX_TO_ID = {v: k for k, v in LIB_ID_TO_INDEX.items()}


def ion_mode_label_to_scalar(value: str | float | int | None) -> float | None:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)

    text = str(value).strip().lower()
    if text in ION_MODE_TO_INDEX:
        return float(ION_MODE_TO_INDEX[text])

    try:
        return float(text)
    except ValueError as exc:
        raise ValueError(
            "Could not parse ion mode value. Expected negative/positive or numeric scalar, "
            f"received: {value!r}"
        ) from exc


def lib_id_to_index(value: int | str) -> int:
    try:
        raw = int(value)
    except Exception as exc:
        raise ValueError(f"Invalid lib_id value: {value!r}") from exc

    if raw in LIB_ID_TO_INDEX:
        return int(LIB_ID_TO_INDEX[raw])
    if raw in LIB_INDEX_TO_ID:
        return int(raw)

    raise ValueError(
        f"Unsupported lib_id {raw}. Expected one of raw ids {sorted(LIB_ID_TO_INDEX.keys())} "
        f"or indices {sorted(LIB_INDEX_TO_ID.keys())}."
    )


class MergedNet1D(Net1D):
    """Mass-conditioned Net1D extended with ion_mode + lib_id embeddings."""

    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__(*args, **kwargs)
        if self.legacy_mode:
            raise ValueError(
                "MergedNet1D only supports the mass-conditioned Net1D constructor. "
                "Legacy mode is not supported."
            )

        self.mode_embedding = nn.Embedding(2, 512)
        self.lib_embedding = nn.Embedding(4, 512)
        self.mode_alpha = nn.Parameter(torch.tensor(1.0, dtype=torch.float32))
        self.lib_alpha = nn.Parameter(torch.tensor(1.0, dtype=torch.float32))

        self.default_ion_mode_idx: int | None = None
        self.default_lib_idx: int | None = None

    def _resolve_ion_mode_indices(
        self,
        x: torch.Tensor,
        ion_mode: torch.Tensor | float | int | None,
    ) -> torch.Tensor | None:
        batch_size = int(x.shape[0])

        if ion_mode is None:
            if self.default_ion_mode_idx is None:
                return None
            return torch.full(
                (batch_size,),
                int(self.default_ion_mode_idx),
                dtype=torch.long,
                device=x.device,
            )

        if isinstance(ion_mode, (float, int)):
            value = int(round(float(ion_mode)))
            value = max(0, min(1, value))
            return torch.full((batch_size,), value, dtype=torch.long, device=x.device)

        mode_t = ion_mode.to(device=x.device)
        if mode_t.ndim == 0:
            mode_t = mode_t.unsqueeze(0)
        if mode_t.ndim == 2 and mode_t.shape[1] == 1:
            mode_t = mode_t.squeeze(1)
        if mode_t.ndim != 1:
            raise ValueError(
                "Expected ion_mode tensor to be scalar, shape (B,), or shape (B,1). "
                f"Received shape {tuple(mode_t.shape)}."
            )

        if mode_t.shape[0] == 1 and batch_size > 1:
            mode_t = mode_t.expand(batch_size)
        if mode_t.shape[0] != batch_size:
            raise ValueError(
                f"Batch mismatch: spectra batch={batch_size}, ion_mode batch={mode_t.shape[0]}."
            )

        if mode_t.dtype.is_floating_point:
            mode_t = torch.round(mode_t).to(dtype=torch.long)
        else:
            mode_t = mode_t.to(dtype=torch.long)

        mode_t = mode_t.clamp(0, 1)
        return mode_t

    def _resolve_lib_indices(
        self,
        x: torch.Tensor,
        lib_id: torch.Tensor | int | None,
    ) -> torch.Tensor | None:
        batch_size = int(x.shape[0])

        if lib_id is None:
            if self.default_lib_idx is None:
                return None
            return torch.full(
                (batch_size,),
                int(self.default_lib_idx),
                dtype=torch.long,
                device=x.device,
            )

        if isinstance(lib_id, int):
            idx = lib_id_to_index(lib_id)
            return torch.full((batch_size,), int(idx), dtype=torch.long, device=x.device)

        lib_t = lib_id.to(device=x.device)
        if lib_t.ndim == 0:
            lib_t = lib_t.unsqueeze(0)
        if lib_t.ndim == 2 and lib_t.shape[1] == 1:
            lib_t = lib_t.squeeze(1)
        if lib_t.ndim != 1:
            raise ValueError(
                "Expected lib_id tensor to be scalar, shape (B,), or shape (B,1). "
                f"Received shape {tuple(lib_t.shape)}."
            )

        if lib_t.shape[0] == 1 and batch_size > 1:
            lib_t = lib_t.expand(batch_size)
        if lib_t.shape[0] != batch_size:
            raise ValueError(
                f"Batch mismatch: spectra batch={batch_size}, lib_id batch={lib_t.shape[0]}."
            )

        if lib_t.dtype.is_floating_point:
            lib_t = torch.round(lib_t).to(dtype=torch.long)
        else:
            lib_t = lib_t.to(dtype=torch.long)

        out = torch.empty_like(lib_t)
        # MPS does not implement aten::_unique2 on some torch builds.
        unique_values = sorted({int(v) for v in lib_t.detach().cpu().tolist()})
        for value in unique_values:
            idx = lib_id_to_index(int(value))
            out[lib_t == int(value)] = int(idx)
        return out

    def forward(
        self,
        x: torch.Tensor,
        parent_mz: torch.Tensor | None = None,
        ion_mode: torch.Tensor | float | int | None = None,
        lib_id: torch.Tensor | int | None = None,
    ) -> torch.Tensor:
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
            if parent_mz.ndim == 2 and parent_mz.shape[1] == 1:
                parent_mz = parent_mz.squeeze(1)
            if parent_mz.ndim != 1:
                raise ValueError(
                    "Expected parent_mz tensor to be scalar, shape (B,), or shape (B,1). "
                    f"Received shape {tuple(parent_mz.shape)}."
                )
            if parent_mz.shape[0] == 1 and x.shape[0] > 1:
                parent_mz = parent_mz.expand(x.shape[0])
            if parent_mz.shape[0] != x.shape[0]:
                raise ValueError(
                    f"Batch mismatch: spectra batch={x.shape[0]}, parent_mz batch={parent_mz.shape[0]}."
                )
            mz_feat = self.mz_mlp(self._normalize_parent_mz(parent_mz))
            x = x + self.mz_alpha * mz_feat

        mode_idx = self._resolve_ion_mode_indices(x, ion_mode)
        if mode_idx is not None:
            mode_feat = self.mode_embedding(mode_idx)
            x = x + self.mode_alpha * mode_feat

        lib_idx = self._resolve_lib_indices(x, lib_id)
        if lib_idx is not None:
            lib_feat = self.lib_embedding(lib_idx)
            x = x + self.lib_alpha * lib_feat

        return torch.tanh(x) if self.apply_tanh else x


def load_compatible_state_dict_for_merged(
    model: MergedNet1D,
    state_dict: dict[str, torch.Tensor],
    *,
    ignore_new_conditioning_params: bool = True,
) -> dict[str, object]:
    """Load checkpoint weights while optionally preserving fresh conditioning params."""
    skipped_conditioning_keys: list[str] = []

    if ignore_new_conditioning_params:
        filtered: dict[str, torch.Tensor] = {}
        for key, value in state_dict.items():
            if key in {"mode_alpha", "lib_alpha"}:
                skipped_conditioning_keys.append(key)
                continue
            if key.startswith("mode_embedding.") or key.startswith("lib_embedding."):
                skipped_conditioning_keys.append(key)
                continue
            filtered[key] = value
    else:
        filtered = state_dict

    info = load_compatible_state_dict(model, filtered)
    if skipped_conditioning_keys:
        info["skipped_conditioning_keys"] = sorted(skipped_conditioning_keys)
    info["ignore_new_conditioning_params"] = bool(ignore_new_conditioning_params)
    return info
