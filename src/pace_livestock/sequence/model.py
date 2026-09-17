"""Safe quantitative adapters. The fixed linear adapter is for deterministic software demos."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..core import bulk_mean
from ..errors import PaceError


def one_hot(sequences: list[str]) -> np.ndarray:
    if not sequences or len({len(s) for s in sequences}) != 1:
        raise PaceError("one_hot requires equally sized, nonempty DNA sequences")
    table = np.zeros((256, 4), dtype=np.float32)
    for index, base in enumerate("ACGT"):
        table[ord(base), index] = 1
    encoded = []
    for seq in sequences:
        seq = seq.upper()
        if set(seq) - set("ACGTN"):
            raise PaceError("DNA input must contain only A/C/G/T/N")
        encoded.append(table[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)].T)
    return np.stack(encoded)


def predict_windows(windows, manifest, *, max_n_fraction=0.05):
    adapter = manifest["adapter"]
    root = Path(manifest["asset_directory"])
    model = None
    if adapter == "fixed_linear":
        if not manifest["is_synthetic"]:
            raise PaceError("fixed_linear is a demonstration-only fixture adapter")
        weights = np.load(root / manifest["weights_file"], allow_pickle=False)
        coefficients, bias = weights["coefficients"], weights["bias"]
        if coefficients.shape != (len(manifest["assays"]), 4) or bias.shape != (
            len(manifest["assays"]),
        ):
            raise PaceError("Fixed model weights have incorrect shape")
    elif adapter == "cnn":
        try:
            import torch
            from safetensors.torch import load_file
        except ImportError as exc:
            raise PaceError(
                "CNN inference requires pip install 'pace-livestock[sequence]'"
            ) from exc
        from .network import QuantitativeCNN

        model = QuantitativeCNN(
            manifest["input_length"],
            manifest["output_window"],
            len(manifest["assays"]),
            channels=manifest.get("channels", [64, 128, 128]),
        )
        model.load_state_dict(load_file(root / manifest["weights_file"]))
        model.eval()
    else:
        raise PaceError(
            "external_table adapter requires inputs.predictions; it cannot infer sequence"
        )
    rows = []
    for window in windows:
        reason = window["reason"]
        values = np.full(len(manifest["assays"]), np.nan)
        seqs = window["sequences"]
        if window["status"] == "resolved":
            if not seqs or any(len(s) != manifest["input_length"] for s in seqs):
                raise PaceError("Prepared window length differs from model input_length")
            if any(s.count("N") / len(s) > max_n_fraction for s in seqs):
                reason = "excess_unknown_bases"
            else:
                x = one_hot(seqs)
                if model is None:
                    left = (manifest["input_length"] - manifest["output_window"]) // 2
                    center = x[:, :, left : left + manifest["output_window"]].mean(axis=2)
                    per_copy = np.logaddexp(0, center @ coefficients.T + bias)
                else:
                    with torch.no_grad():
                        per_copy = model(torch.from_numpy(x)).numpy()
                values = bulk_mean(per_copy)
        for assay, value in zip(manifest["assays"], values, strict=True):
            rows.append(
                {
                    "element_id": window["element_id"],
                    "assay": assay,
                    "predicted_value": float(value),
                    "model_id": manifest["model_id"],
                    "unit": manifest["signal_unit"],
                    "normalization_id": manifest["normalization_id"],
                    "window_id": f"grid:{manifest['output_window']}:mean",
                    "status": "resolved" if np.isfinite(value) else "unresolved",
                    "reason": reason,
                    "output_target_id": window["output_target_id"],
                    "structural_status": window["structural_status"],
                    "callable_fraction": window["callable_fraction"],
                    "assumed_reference_fraction": window["assumed_reference_fraction"],
                }
            )
    return rows
