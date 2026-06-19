#!/usr/bin/env python3
"""Batch inference for merged-library spectral features.

Exposes the encoder–decoder prediction pipeline (``predict``) and the
merged-ionization CLI wrapper (``merged``) that injects ``lib_id`` /
``ion_mode`` conditioning into ``MergedNet1D``.
"""
