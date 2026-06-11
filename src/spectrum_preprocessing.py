from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Iterable

import numpy as np

from utils import spec2vec


@dataclass(frozen=True)
class SpectrumPreprocessingConfig:
    min_rel_intensity_pct: float = 2.0
    top_peaks: int = 25
    disable_basepeak_norm: bool = False

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["basepeak_normalization_enabled"] = not bool(self.disable_basepeak_norm)
        return payload


def _clean_peak_tuple(value: tuple[object, object]) -> tuple[float, float] | None:
    try:
        mz = float(value[0])
        intensity = float(value[1])
    except Exception:
        return None
    if not math.isfinite(mz) or not math.isfinite(intensity):
        return None
    if intensity <= 0:
        return None
    return mz, intensity


def sanitize_peaks(peaks: Iterable[tuple[object, object]]) -> list[tuple[float, float]]:
    cleaned: list[tuple[float, float]] = []
    for peak in peaks:
        parsed = _clean_peak_tuple(peak)
        if parsed is not None:
            cleaned.append(parsed)
    return cleaned


def basepeak_normalize(peaks: Iterable[tuple[object, object]]) -> list[tuple[float, float]]:
    cleaned = sanitize_peaks(peaks)
    if not cleaned:
        return []
    basepeak = max(intensity for _, intensity in cleaned)
    if basepeak <= 0:
        return []
    return [(mz, intensity / basepeak) for mz, intensity in cleaned]


def preprocess_peaks(
    peaks: Iterable[tuple[object, object]],
    config: SpectrumPreprocessingConfig,
) -> dict[str, object]:
    raw_peaks = sanitize_peaks(peaks)
    if not raw_peaks:
        return {
            "raw_peaks": [],
            "peaks": [],
            "n_raw_peaks": 0,
            "n_retained_peaks": 0,
            "basepeak_intensity": 0.0,
        }

    basepeak = max(intensity for _, intensity in raw_peaks)
    threshold = max(float(config.min_rel_intensity_pct), 0.0) / 100.0 * basepeak
    retained = [(mz, intensity) for mz, intensity in raw_peaks if intensity >= threshold]

    if not bool(config.disable_basepeak_norm):
        retained = [(mz, intensity / basepeak) for mz, intensity in retained]

    retained.sort(key=lambda item: (-item[1], item[0]))
    if int(config.top_peaks) > 0:
        retained = retained[: int(config.top_peaks)]

    return {
        "raw_peaks": raw_peaks,
        "peaks": retained,
        "n_raw_peaks": len(raw_peaks),
        "n_retained_peaks": len(retained),
        "basepeak_intensity": float(basepeak),
    }


def vectorize_processed_peaks(
    processed_peaks: Iterable[tuple[float, float]],
    minmass: float,
    maxmass: float,
    resolution: int,
) -> dict[str, object]:
    n_bins = int((float(maxmass) - float(minmass)) * int(resolution)) + 1
    peaks = [(float(mz), float(intensity)) for mz, intensity in processed_peaks]
    if not peaks:
        return {
            "vector": np.zeros(n_bins, dtype=np.float32),
            "peaks": [],
            "n_vectorized_peaks": 0,
            "n_out_of_range_retained_peaks": 0,
            "dropped_out_of_range_peak_fraction": 0.0,
        }

    in_range = [
        (mz, intensity)
        for mz, intensity in peaks
        if float(minmass) <= float(mz) <= float(maxmass)
    ]
    n_out_of_range = int(len(peaks) - len(in_range))
    dropped_fraction = float(n_out_of_range / max(len(peaks), 1))
    if not in_range:
        return {
            "vector": np.zeros(n_bins, dtype=np.float32),
            "peaks": [],
            "n_vectorized_peaks": 0,
            "n_out_of_range_retained_peaks": n_out_of_range,
            "dropped_out_of_range_peak_fraction": dropped_fraction,
        }

    peaks_set = set(in_range)
    if not peaks_set:
        return {
            "vector": np.zeros(n_bins, dtype=np.float32),
            "peaks": [],
            "n_vectorized_peaks": 0,
            "n_out_of_range_retained_peaks": n_out_of_range,
            "dropped_out_of_range_peak_fraction": dropped_fraction,
        }

    vec = spec2vec(
        peaks_set,
        minmass=float(minmass),
        maxmass=float(maxmass),
        resolution=int(resolution),
    )
    return {
        "vector": vec.astype(np.float32),
        "peaks": sorted(peaks_set, key=lambda item: (-item[1], item[0])),
        "n_vectorized_peaks": int(len(peaks_set)),
        "n_out_of_range_retained_peaks": n_out_of_range,
        "dropped_out_of_range_peak_fraction": dropped_fraction,
    }


def preprocess_and_vectorize(
    peaks: Iterable[tuple[object, object]],
    config: SpectrumPreprocessingConfig,
    minmass: float,
    maxmass: float,
    resolution: int,
) -> dict[str, object]:
    result = preprocess_peaks(peaks=peaks, config=config)
    retained_before_range = list(result["peaks"])
    vectorized = vectorize_processed_peaks(
        processed_peaks=result["peaks"],
        minmass=minmass,
        maxmass=maxmass,
        resolution=resolution,
    )
    result["peaks_before_mz_range"] = retained_before_range
    result["peaks"] = list(vectorized["peaks"])
    result["vector"] = vectorized["vector"]
    result["n_vectorized_peaks"] = int(vectorized["n_vectorized_peaks"])
    result["n_out_of_range_retained_peaks"] = int(
        vectorized["n_out_of_range_retained_peaks"]
    )
    result["dropped_out_of_range_peak_fraction"] = float(
        vectorized["dropped_out_of_range_peak_fraction"]
    )
    return result


def summarize_peak_counts(values: Iterable[int]) -> dict[str, float | int]:
    counts = np.asarray([int(v) for v in values], dtype=np.int64)
    if counts.size == 0:
        return {
            "count": 0,
            "min": 0,
            "q25": 0,
            "median": 0,
            "q75": 0,
            "q90": 0,
            "max": 0,
            "mean": 0.0,
        }
    return {
        "count": int(counts.size),
        "min": int(np.min(counts)),
        "q25": int(np.percentile(counts, 25)),
        "median": int(np.percentile(counts, 50)),
        "q75": int(np.percentile(counts, 75)),
        "q90": int(np.percentile(counts, 90)),
        "max": int(np.max(counts)),
        "mean": float(np.mean(counts)),
    }
