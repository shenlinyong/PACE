"""Deterministic CPU CNN training, frozen scales, group isolation and safe persistence."""

from __future__ import annotations

import copy

import numpy as np

from ..config import operation_config
from ..errors import PaceError
from ..io.tables import integer, number, read_table
from ..provenance import file_hash, output_directory, write_json
from .model import one_hot


def validate_splits(rows, *, input_length, output_window):
    groups, intervals = {}, {}
    for row in rows:
        role = row["split"]
        if role not in ("train", "validation", "test"):
            raise PaceError("Sequence split must be train, validation or test")
        group = row["group_id"]
        if group in groups and groups[group] != role:
            raise PaceError("Sequence group leaks across train/validation/test")
        groups[group] = role
        start, end = integer(row["start"], "training start"), integer(row["end"], "training end")
        if end - start != output_window:
            raise PaceError("Training target interval differs from output_window")
        flank = (input_length - output_window) // 2
        intervals.setdefault(row["chrom"], []).append((start - flank, end + flank, role))
    for values in intervals.values():
        active = []
        for start, end, role in sorted(values):
            active = [(e, r) for e, r in active if e > start]
            if any(r != role for _, r in active):
                raise PaceError("Overlapping sequence context crosses a split boundary")
            active.append((end, role))


def train_sequence(config_path, out):
    try:
        import torch
        from safetensors.torch import save_file
    except ImportError as exc:
        raise PaceError("train-sequence requires pace-livestock[sequence]") from exc
    from .network import QuantitativeCNN, masked_huber

    allowed = {
        "data",
        "model_id",
        "species",
        "assembly",
        "context_id",
        "target_level",
        "assays",
        "signal_unit",
        "normalization_id",
        "is_synthetic",
        "input_length",
        "output_window",
        "seed",
        "epochs",
        "batch_size",
        "learning_rate",
        "patience",
        "channels",
        "threads",
        "max_n_fraction",
    }
    cfg = operation_config(
        config_path,
        allowed=allowed,
        required=allowed
        - {
            "input_length",
            "output_window",
            "seed",
            "epochs",
            "batch_size",
            "learning_rate",
            "patience",
            "channels",
            "threads",
            "max_n_fraction",
        },
        paths=["data"],
    )
    seed = integer(cfg.get("seed", 17), "seed")
    torch.manual_seed(seed)
    torch.set_num_threads(integer(cfg.get("threads", 1), "threads", minimum=1))
    torch.use_deterministic_algorithms(True)
    length, window = cfg.get("input_length", 8192), cfg.get("output_window", 500)
    assays = cfg["assays"]
    if type(cfg["is_synthetic"]) is not bool or not assays or len(assays) != len(set(assays)):
        raise PaceError("Training requires explicit is_synthetic and unique assay heads")
    rows = read_table(
        cfg["data"],
        required=["sequence_id", "sequence", "split", "group_id", "chrom", "start", "end", *assays],
    )
    validate_splits(rows, input_length=length, output_window=window)
    if any(len(r["sequence"]) != length for r in rows):
        raise PaceError("Training sequences must match input_length")
    max_n = number(cfg.get("max_n_fraction", 0.05), "max_n_fraction", minimum=0, maximum=1)
    if any(r["sequence"].upper().count("N") / length > max_n for r in rows):
        raise PaceError("Training sequence exceeds max_n_fraction")
    labels = np.array(
        [[number(r[a], f"label {a}", missing=True, minimum=0) for a in assays] for r in rows],
        dtype=np.float32,
    )
    train = np.array([i for i, r in enumerate(rows) if r["split"] == "train"])
    val = np.array([i for i, r in enumerate(rows) if r["split"] == "validation"])
    if not len(train) or not len(val):
        raise PaceError("Training and held-out validation rows are required")
    scales = []
    for j, assay in enumerate(assays):
        positives = labels[train, j][np.isfinite(labels[train, j]) & (labels[train, j] > 0)]
        if not len(positives):
            raise PaceError(f"Assay head {assay} has no positive training labels")
        scales.append(float(np.median(positives)))
    x = torch.from_numpy(one_hot([r["sequence"] for r in rows]))
    y, s = torch.from_numpy(labels), torch.tensor(scales)
    channels = cfg.get("channels", [64, 128, 128])
    model = QuantitativeCNN(length, window, len(assays), channels=channels)
    lr = number(cfg.get("learning_rate", 0.001), "learning_rate", minimum=0)
    if lr == 0:
        raise PaceError("learning_rate must be positive")
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr)
    batch = integer(cfg.get("batch_size", 64), "batch_size", minimum=1)
    epochs = integer(cfg.get("epochs", 50), "epochs", minimum=1)
    patience = integer(cfg.get("patience", 5), "patience", minimum=1)
    initial = {k: v.detach().clone() for k, v in model.state_dict().items()}
    best, best_loss, stale, history = None, float("inf"), 0, []
    rng = np.random.default_rng(seed)
    for epoch in range(epochs):
        model.train()
        losses = []
        order = rng.permutation(train)
        for offset in range(0, len(order), batch):
            idx = order[offset : offset + batch]
            if not torch.isfinite(y[idx]).any():
                continue
            xb = x[idx].clone()
            reverse = torch.from_numpy(rng.random(len(idx)) < 0.5)
            xb[reverse] = torch.flip(xb[reverse], dims=[1, 2])
            optimizer.zero_grad()
            loss = masked_huber(model(xb), y[idx], s)
            loss.backward()
            if not torch.isfinite(loss) or any(
                not torch.isfinite(p.grad).all() for p in model.parameters() if p.grad is not None
            ):
                raise PaceError("Nonfinite sequence training loss/gradient")
            optimizer.step()
            losses.append(float(loss.detach()))
        model.eval()
        with torch.no_grad():
            val_loss = float(masked_huber(model(x[val]), y[val], s))
        history.append(
            {"epoch": epoch + 1, "train_loss": float(np.mean(losses)), "validation_loss": val_loss}
        )
        if val_loss < best_loss:
            best, best_loss, stale = copy.deepcopy(model.state_dict()), val_loss, 0
        else:
            stale += 1
        if stale >= patience:
            break
    if best is None:
        raise PaceError("Sequence training produced no finite checkpoint")
    changed = any(not torch.equal(initial[k], best[k]) for k in initial)
    if not changed:
        raise PaceError("Training did not update any model parameters")
    m = {
        k: cfg[k]
        for k in (
            "model_id",
            "species",
            "assembly",
            "context_id",
            "target_level",
            "assays",
            "signal_unit",
            "normalization_id",
            "is_synthetic",
        )
    }
    m.update(
        {
            "kind": "sequence",
            "adapter": "cnn",
            "input_length": length,
            "output_window": window,
            "output_type": "normalized_signal",
            "channels": channels,
            "scales": scales,
            "weights_file": "weights.safetensors",
            "training_sources": [file_hash(cfg["data"])],
            "calibration_sources": [],
            "test_sources": [],
            "validation": {},
            "seed": seed,
            "padding": "same_convolution; nonoverlapping_avgpool4; stride64",
            "target_pooling": "feature_bin_overlap_with_central_target",
            "software_test_only": cfg["is_synthetic"],
        }
    )
    with output_directory(out) as destination:
        save_file({k: v.contiguous() for k, v in best.items()}, destination / "weights.safetensors")
        m["weights_sha256"] = file_hash(destination / "weights.safetensors")
        write_json(destination / "manifest.json", m)
        write_json(
            destination / "training_report.json",
            {
                "history": history,
                "parameters_updated": changed,
                "training_config": cfg,
                "n_train": len(train),
                "n_validation": len(val),
                "test_rows_not_used": sum(r["split"] == "test" for r in rows),
                "biological_validation": "not_assessed",
            },
        )
    return m
