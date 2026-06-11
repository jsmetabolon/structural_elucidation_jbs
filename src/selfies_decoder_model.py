#!/usr/bin/env python3
"""Standalone SELFIES decoder with direct mass conditioning."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    import selfies as sf
except ImportError as exc:  # pragma: no cover - runtime dependency
    raise RuntimeError(
        "Package `selfies` is required for SELFIES decoder."
    ) from exc

try:
    from rdkit import Chem

    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


def canonicalize_smiles(smiles: str) -> str | None:
    if not HAS_RDKIT or not smiles or not str(smiles).strip():
        return None
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return None


def tokenize_selfies(selfies_text: str) -> list[str]:
    if not selfies_text:
        return []
    return list(sf.split_selfies(selfies_text))


@dataclass(frozen=True)
class SelfiesVocab:
    token_to_id: dict[str, int]
    id_to_token: list[str]
    pad_token: str = "[PAD]"
    bos_token: str = "[BOS]"
    eos_token: str = "[EOS]"
    unk_token: str = "[UNK]"

    @property
    def pad_id(self) -> int:
        return int(self.token_to_id[self.pad_token])

    @property
    def bos_id(self) -> int:
        return int(self.token_to_id[self.bos_token])

    @property
    def eos_id(self) -> int:
        return int(self.token_to_id[self.eos_token])

    @property
    def unk_id(self) -> int:
        return int(self.token_to_id[self.unk_token])

    def __len__(self) -> int:
        return int(len(self.id_to_token))

    def encode(
        self, selfies_text: str, add_bos: bool = True, add_eos: bool = True
    ) -> list[int]:
        tokens = tokenize_selfies(selfies_text)
        ids = [self.token_to_id.get(token, self.unk_id) for token in tokens]
        if add_bos:
            ids = [self.bos_id] + ids
        if add_eos:
            ids = ids + [self.eos_id]
        return ids

    def decode(self, ids: Sequence[int], skip_special: bool = True) -> str:
        specials = {self.pad_token, self.bos_token, self.eos_token}
        tokens: list[str] = []
        for idx in ids:
            token = self.id_to_token[int(idx)]
            if token == self.eos_token:
                break
            if skip_special and token in specials:
                continue
            tokens.append(token)
        return "".join(tokens)

    def to_dict(self) -> dict[str, object]:
        return {
            "id_to_token": list(self.id_to_token),
            "pad_token": self.pad_token,
            "bos_token": self.bos_token,
            "eos_token": self.eos_token,
            "unk_token": self.unk_token,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, object]) -> "SelfiesVocab":
        id_to_token = [str(token) for token in payload["id_to_token"]]
        token_to_id = {token: idx for idx, token in enumerate(id_to_token)}
        pad_token = str(payload.get("pad_token", "[PAD]"))
        bos_token = str(payload.get("bos_token", "[BOS]"))
        eos_token = str(payload.get("eos_token", "[EOS]"))
        unk_token = str(payload.get("unk_token", "[UNK]"))
        for token in (pad_token, bos_token, eos_token, unk_token):
            if token not in token_to_id:
                raise ValueError(f"Vocabulary missing special token: {token}")
        return cls(
            token_to_id=token_to_id,
            id_to_token=id_to_token,
            pad_token=pad_token,
            bos_token=bos_token,
            eos_token=eos_token,
            unk_token=unk_token,
        )

    @classmethod
    def build(cls, selfies_strings: Iterable[str]) -> "SelfiesVocab":
        unique_tokens: set[str] = set()
        for text in selfies_strings:
            unique_tokens.update(tokenize_selfies(text))
        ordered_tokens = sorted(unique_tokens)
        id_to_token = ["[PAD]", "[BOS]", "[EOS]", "[UNK]"] + ordered_tokens
        token_to_id = {token: idx for idx, token in enumerate(id_to_token)}
        return cls(token_to_id=token_to_id, id_to_token=id_to_token)

    def save_json(self, path: Path) -> None:
        path = Path(path).expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2))

    @classmethod
    def load_json(cls, path: Path) -> "SelfiesVocab":
        payload = json.loads(Path(path).expanduser().resolve().read_text())
        if not isinstance(payload, dict):
            raise ValueError(f"Vocabulary JSON must contain an object: {path}")
        return cls.from_dict(payload)


class SelfiesDecoder(nn.Module):
    def __init__(
        self,
        vocab: SelfiesVocab,
        d_z: int = 512,
        d_z_proj: int = 384,
        d_mass_emb: int = 128,
        d_ri_emb: int = 0,
        use_ri_conditioning: bool = False,
        formula_input_dim: int = 0,
        adduct_input_dim: int = 0,
        substructure_input_dim: int = 0,
        use_formula_conditioning: bool = False,
        use_adduct_conditioning: bool = False,
        use_substructure_conditioning: bool = False,
        d_substructure_emb: int = 0,
        use_class_conditioning: bool = False,
        num_classes: int = 9,
        d_class_emb: int = 32,
        d_token_emb: int = 64,
        d_gru_hidden: int = 512,
        n_gru_layers: int = 2,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.vocab = vocab
        self.pad_id = int(vocab.pad_id)
        self.bos_id = int(vocab.bos_id)
        self.eos_id = int(vocab.eos_id)

        self.token_embedding = nn.Embedding(
            len(vocab), d_token_emb, padding_idx=self.pad_id
        )
        self.z_projection = nn.Linear(int(d_z), int(d_z_proj))
        self.use_formula_conditioning = bool(use_formula_conditioning)
        self.formula_input_dim = (
            int(formula_input_dim) if self.use_formula_conditioning else 0
        )
        self.use_adduct_conditioning = bool(use_adduct_conditioning)
        self.adduct_input_dim = (
            int(adduct_input_dim) if self.use_adduct_conditioning else 0
        )
        self.use_substructure_conditioning = bool(use_substructure_conditioning)
        self.substructure_input_dim = (
            int(substructure_input_dim) if self.use_substructure_conditioning else 0
        )
        self.d_substructure_emb = (
            int(d_substructure_emb) if self.use_substructure_conditioning else 0
        )
        self.use_class_conditioning = bool(use_class_conditioning)
        self.num_classes = int(num_classes) if self.use_class_conditioning else 0
        self.d_class_emb = int(d_class_emb) if self.use_class_conditioning else 0
        if self.use_class_conditioning:
            if self.num_classes <= 0:
                raise ValueError("num_classes must be >0 when class conditioning is enabled.")
            if self.d_class_emb <= 0:
                self.d_class_emb = 32
            self.class_embedding = nn.Embedding(self.num_classes, self.d_class_emb)
            self.class_projection = nn.Sequential(
                nn.Linear(self.d_class_emb, self.d_class_emb),
                nn.ReLU(),
            )
        else:
            self.class_embedding = None
            self.class_projection = None
        if self.use_substructure_conditioning:
            if self.substructure_input_dim <= 0:
                raise ValueError(
                    "substructure_input_dim must be >0 when substructure conditioning is enabled."
                )
            if self.d_substructure_emb <= 0:
                self.d_substructure_emb = min(128, max(32, self.substructure_input_dim))
            self.substructure_mlp = nn.Sequential(
                nn.Linear(self.substructure_input_dim, self.d_substructure_emb),
                nn.ReLU(),
                nn.Linear(self.d_substructure_emb, self.d_substructure_emb),
                nn.ReLU(),
            )
        else:
            self.substructure_mlp = None
        condition_input_dim = 1 + int(self.formula_input_dim) + int(self.adduct_input_dim)
        self.mass_mlp = nn.Sequential(
            nn.Linear(condition_input_dim, 64),
            nn.ReLU(),
            nn.Linear(64, int(d_mass_emb)),
            nn.ReLU(),
        )
        self.formula_head = (
            nn.Linear(int(d_z_proj), int(self.formula_input_dim))
            if self.use_formula_conditioning and self.formula_input_dim > 0
            else None
        )
        self.use_ri_conditioning = bool(use_ri_conditioning)
        self.d_ri_emb = int(d_ri_emb) if self.use_ri_conditioning else 0
        if self.use_ri_conditioning:
            if self.d_ri_emb <= 0:
                self.d_ri_emb = 32
            self.ri_mlp = nn.Sequential(
                nn.Linear(1, 32),
                nn.ReLU(),
                nn.Linear(32, int(self.d_ri_emb)),
                nn.ReLU(),
            )
        else:
            self.ri_mlp = None
        self.gru_input_dim = (
            int(d_token_emb)
            + int(d_z_proj)
            + int(d_mass_emb)
            + int(self.d_ri_emb)
            + int(self.d_substructure_emb)
            + int(self.d_class_emb)
        )
        self.gru = nn.GRU(
            input_size=self.gru_input_dim,
            hidden_size=int(d_gru_hidden),
            num_layers=int(n_gru_layers),
            dropout=float(dropout) if int(n_gru_layers) > 1 else 0.0,
            batch_first=True,
        )
        self.output_head = nn.Linear(int(d_gru_hidden), len(vocab))
        self.mass_head = nn.Linear(int(d_z_proj), 1)

    @property
    def device(self) -> torch.device:
        return next(self.parameters()).device

    def _mass_embedding(
        self,
        mass: torch.Tensor,
        formula: torch.Tensor | None = None,
        adduct: torch.Tensor | None = None,
    ) -> torch.Tensor:
        mass = mass.float()
        if mass.ndim == 1:
            mass = mass.unsqueeze(1)
        if mass.ndim != 2 or mass.shape[1] != 1:
            raise ValueError(
                f"Expected mass shape (B,1) or (B,), got {tuple(mass.shape)}"
            )
        parts = [mass]
        batch_size = int(mass.shape[0])
        if self.use_formula_conditioning:
            if formula is None:
                formula = torch.zeros(
                    (batch_size, int(self.formula_input_dim)),
                    dtype=torch.float32,
                    device=mass.device,
                )
            formula = formula.float().to(mass.device)
            if formula.ndim == 1:
                formula = formula.unsqueeze(0)
            if formula.ndim != 2 or formula.shape != (
                batch_size,
                int(self.formula_input_dim),
            ):
                raise ValueError(
                    "Expected formula shape "
                    f"{(batch_size, int(self.formula_input_dim))}, got {tuple(formula.shape)}"
                )
            parts.append(formula)
        if self.use_adduct_conditioning:
            if adduct is None:
                adduct = torch.zeros(
                    (batch_size, int(self.adduct_input_dim)),
                    dtype=torch.float32,
                    device=mass.device,
                )
            adduct = adduct.float().to(mass.device)
            if adduct.ndim == 1:
                adduct = adduct.unsqueeze(0)
            if adduct.ndim != 2 or adduct.shape != (
                batch_size,
                int(self.adduct_input_dim),
            ):
                raise ValueError(
                    "Expected adduct shape "
                    f"{(batch_size, int(self.adduct_input_dim))}, got {tuple(adduct.shape)}"
                )
            parts.append(adduct)
        return self.mass_mlp(torch.cat(parts, dim=1))

    def _ri_embedding(
        self,
        ri: torch.Tensor | None,
        batch_size: int,
    ) -> torch.Tensor | None:
        if not bool(self.use_ri_conditioning):
            return None
        if ri is None:
            ri = torch.zeros((int(batch_size), 1), dtype=torch.float32, device=self.device)
        ri = ri.float().to(self.device)
        if ri.ndim == 1:
            ri = ri.unsqueeze(1)
        if ri.ndim != 2 or ri.shape != (int(batch_size), 1):
            raise ValueError(
                f"Expected ri shape {(int(batch_size), 1)} or ({int(batch_size)},), got {tuple(ri.shape)}"
            )
        assert self.ri_mlp is not None
        return self.ri_mlp(ri)

    def _substructure_embedding(
        self,
        substructure: torch.Tensor | None,
        batch_size: int,
    ) -> torch.Tensor | None:
        if not bool(self.use_substructure_conditioning):
            return None
        if substructure is None:
            substructure = torch.zeros(
                (int(batch_size), int(self.substructure_input_dim)),
                dtype=torch.float32,
                device=self.device,
            )
        substructure = substructure.float().to(self.device)
        if substructure.ndim == 1:
            substructure = substructure.unsqueeze(0)
        if substructure.ndim != 2 or substructure.shape != (
            int(batch_size),
            int(self.substructure_input_dim),
        ):
            raise ValueError(
                "Expected substructure shape "
                f"{(int(batch_size), int(self.substructure_input_dim))}, got {tuple(substructure.shape)}"
            )
        assert self.substructure_mlp is not None
        return self.substructure_mlp(substructure)

    def _class_condition_embedding(
        self,
        class_token: torch.Tensor | None,
        batch_size: int,
    ) -> torch.Tensor | None:
        if not bool(self.use_class_conditioning):
            return None
        if class_token is None:
            class_token = torch.full(
                (int(batch_size),),
                fill_value=int(self.num_classes) - 1,
                dtype=torch.long,
                device=self.device,
            )
        class_token = class_token.long().to(self.device)
        if class_token.ndim == 2 and class_token.shape[1] == 1:
            class_token = class_token.reshape(-1)
        if class_token.ndim != 1 or class_token.shape[0] != int(batch_size):
            raise ValueError(
                f"Expected class_token shape ({int(batch_size)},), got {tuple(class_token.shape)}"
            )
        class_token = class_token.clamp(0, int(self.num_classes) - 1)
        assert self.class_embedding is not None
        assert self.class_projection is not None
        return self.class_projection(self.class_embedding(class_token))

    def forward(
        self,
        z: torch.Tensor,
        mass: torch.Tensor,
        target_tokens: torch.Tensor,
        mass_target: torch.Tensor | None = None,
        ri: torch.Tensor | None = None,
        formula: torch.Tensor | None = None,
        formula_target: torch.Tensor | None = None,
        adduct: torch.Tensor | None = None,
        substructure: torch.Tensor | None = None,
        class_token: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        if z.ndim != 2:
            raise ValueError(f"Expected z shape (B, D), got {tuple(z.shape)}")
        if target_tokens.ndim != 2:
            raise ValueError(
                f"Expected target_tokens shape (B, L), got {tuple(target_tokens.shape)}"
            )
        if target_tokens.shape[1] < 2:
            raise ValueError(
                "target_tokens must include at least BOS and EOS."
            )

        z_proj = self.z_projection(z.float())
        mass_emb = self._mass_embedding(mass, formula=formula, adduct=adduct)
        ri_emb = self._ri_embedding(ri, int(z.shape[0]))
        substructure_emb = self._substructure_embedding(substructure, int(z.shape[0]))
        class_emb = self._class_condition_embedding(class_token, int(z.shape[0]))

        decoder_in = target_tokens[:, :-1].long()
        decoder_out = target_tokens[:, 1:].long()
        token_emb = self.token_embedding(decoder_in)

        seq_len = token_emb.shape[1]
        z_rep = z_proj.unsqueeze(1).expand(-1, seq_len, -1)
        mass_rep = mass_emb.unsqueeze(1).expand(-1, seq_len, -1)
        parts = [token_emb, z_rep, mass_rep]
        if ri_emb is not None:
            parts.append(ri_emb.unsqueeze(1).expand(-1, seq_len, -1))
        if substructure_emb is not None:
            parts.append(substructure_emb.unsqueeze(1).expand(-1, seq_len, -1))
        if class_emb is not None:
            parts.append(class_emb.unsqueeze(1).expand(-1, seq_len, -1))
        gru_input = torch.cat(parts, dim=-1)

        gru_out, _ = self.gru(gru_input)
        logits = self.output_head(gru_out)

        recon_loss = F.cross_entropy(
            logits.reshape(-1, logits.shape[-1]),
            decoder_out.reshape(-1),
            ignore_index=self.pad_id,
        )
        if mass_target is None:
            mass_target = mass
        mass_target = mass_target.reshape(-1).float()
        mass_pred = self.mass_head(z_proj).reshape(-1)
        mass_loss = F.mse_loss(mass_pred, mass_target)
        if self.formula_head is not None and formula_target is not None:
            formula_target = formula_target.float().to(z_proj.device)
            if formula_target.ndim == 1:
                formula_target = formula_target.unsqueeze(0)
            formula_pred = self.formula_head(z_proj)
            formula_loss = F.mse_loss(formula_pred, formula_target)
        else:
            formula_loss = z_proj.new_tensor(0.0)
        return recon_loss, mass_loss, formula_loss

    @torch.no_grad()
    def generate(
        self,
        z: torch.Tensor,
        mass: torch.Tensor,
        temp: float = 1.0,
        max_len: int = 100,
        top_p: float = 0.0,
        ri: torch.Tensor | None = None,
        formula: torch.Tensor | None = None,
        adduct: torch.Tensor | None = None,
        substructure: torch.Tensor | None = None,
        class_token: torch.Tensor | None = None,
        return_log_probs: bool = False,
    ) -> tuple[list[str], list[float]] | list[str]:
        self.eval()
        z = z.float().to(self.device)
        if z.ndim != 2:
            raise ValueError(f"Expected z shape (B, D), got {tuple(z.shape)}")

        mass = mass.float().to(self.device)
        if mass.ndim == 1:
            mass = mass.unsqueeze(1)
        if mass.ndim != 2 or mass.shape != (z.shape[0], 1):
            raise ValueError(
                f"Expected mass shape {(z.shape[0], 1)} or ({z.shape[0]},), got {tuple(mass.shape)}"
            )

        batch_size = z.shape[0]
        z_proj = self.z_projection(z)
        mass_emb = self._mass_embedding(mass, formula=formula, adduct=adduct)
        ri_emb = self._ri_embedding(ri, int(batch_size))
        substructure_emb = self._substructure_embedding(substructure, int(batch_size))
        class_emb = self._class_condition_embedding(class_token, int(batch_size))

        current_tokens = torch.full(
            (batch_size,),
            fill_value=self.bos_id,
            dtype=torch.long,
            device=self.device,
        )
        finished = torch.zeros(
            batch_size, dtype=torch.bool, device=self.device
        )
        sequences: list[list[int]] = [[] for _ in range(batch_size)]
        hidden = None

        log_prob_sums = torch.zeros(
            batch_size, dtype=torch.float32, device=self.device
        )
        token_counts = torch.zeros(
            batch_size, dtype=torch.float32, device=self.device
        )

        for _ in range(int(max_len)):
            token_emb = self.token_embedding(current_tokens).unsqueeze(1)
            parts = [token_emb, z_proj.unsqueeze(1), mass_emb.unsqueeze(1)]
            if ri_emb is not None:
                parts.append(ri_emb.unsqueeze(1))
            if substructure_emb is not None:
                parts.append(substructure_emb.unsqueeze(1))
            if class_emb is not None:
                parts.append(class_emb.unsqueeze(1))
            step_input = torch.cat(parts, dim=-1)
            gru_out, hidden = self.gru(step_input, hidden)
            logits = self.output_head(gru_out.squeeze(1))

            if float(temp) <= 0:
                next_tokens = torch.argmax(logits, dim=-1)
                step_log_probs = torch.zeros_like(log_prob_sums)
            else:
                scaled_logits = logits / max(float(temp), 1e-6)
                probs = torch.softmax(scaled_logits, dim=-1)
                nucleus_p = float(top_p)
                if 0.0 < nucleus_p < 1.0:
                    sorted_probs, sorted_idx = torch.sort(probs, dim=-1, descending=True)
                    cumulative = torch.cumsum(sorted_probs, dim=-1)
                    drop_mask = cumulative > nucleus_p
                    drop_mask[:, 0] = False
                    filtered_sorted_probs = sorted_probs.masked_fill(drop_mask, 0.0)
                    filtered_probs = torch.zeros_like(probs).scatter(1, sorted_idx, filtered_sorted_probs)
                    filtered_sums = filtered_probs.sum(dim=-1, keepdim=True)
                    probs = torch.where(
                        filtered_sums > 0,
                        filtered_probs / filtered_sums.clamp_min(1e-12),
                        probs,
                    )
                next_tokens = torch.multinomial(probs, num_samples=1).squeeze(
                    1
                )
                chosen = torch.gather(
                    probs, 1, next_tokens.unsqueeze(1)
                ).squeeze(1)
                step_log_probs = torch.log(chosen.clamp_min(1e-12))

            active = ~finished
            if active.any():
                log_prob_sums[active] += step_log_probs[active]
                token_counts[active] += 1.0

            next_tokens = torch.where(
                finished,
                torch.full_like(next_tokens, self.eos_id),
                next_tokens,
            )
            for row_idx in range(batch_size):
                token_id = int(next_tokens[row_idx].item())
                sequences[row_idx].append(token_id)

            finished = finished | (next_tokens == self.eos_id)
            current_tokens = next_tokens
            if bool(finished.all()):
                break

        selfies_outputs = [
            self.vocab.decode(token_ids, skip_special=True)
            for token_ids in sequences
        ]
        if not return_log_probs:
            return selfies_outputs

        mean_log_probs = (
            (log_prob_sums / token_counts.clamp_min(1.0))
            .detach()
            .cpu()
            .tolist()
        )
        return selfies_outputs, [float(value) for value in mean_log_probs]

    def selfies_to_smiles(self, selfies_text: str) -> str | None:
        if not selfies_text:
            return None
        try:
            smiles = sf.decoder(selfies_text)
        except Exception:
            return None
        return canonicalize_smiles(smiles)
