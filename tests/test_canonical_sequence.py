"""CPU training smoke and central pooling verification; no biological accuracy threshold."""

import importlib
import json

import numpy as np
import pytest
import yaml

from pace_livestock.errors import PaceError
from pace_livestock.io.tables import write_table
from pace_livestock.sequence.training import train_sequence, validate_splits

torch = pytest.importorskip(
    "torch", reason="install pace-livestock[sequence] for quantitative CNN tests"
)
pytest.importorskip("safetensors", reason="safe sequence weights require safetensors")
network = importlib.import_module("pace_livestock.sequence.network")
QuantitativeCNN, masked_huber = network.QuantitativeCNN, network.masked_huber


def test_masked_loss_and_gradients():
    predicted = torch.tensor([[2.0, 3.0], [4.0, 5.0]], requires_grad=True)
    labels = torch.tensor([[2.0, float("nan")], [1.0, float("nan")]])
    scales = torch.ones(2)
    loss = masked_huber(predicted, labels, scales)
    loss.backward()
    assert torch.isfinite(loss) and torch.isfinite(predicted.grad).all()
    # No label in head 2: its loss gradient must be exactly zero.
    assert torch.equal(predicted.grad[:, 1], torch.zeros(2))
    assert predicted.grad[1, 0] != 0
    with pytest.raises(PaceError, match="no finite labels"):
        masked_huber(predicted, torch.full_like(labels, float("nan")), scales)


def test_central_pooling_responds_to_position():
    torch.set_num_threads(1)
    model = QuantitativeCNN(1024, 128, 1, channels=[2, 2, 2])
    with torch.no_grad():
        for p in model.parameters():
            p.fill_(0.1)
        # Use only central features, remove the context branch in the final linear layer.
        model.head[0].weight[:, 2:] = 0
    center, distal = torch.zeros(1, 4, 1024), torch.zeros(1, 4, 1024)
    center[:, :, 448:576] = 1
    distal[:, :, 0:128] = 1
    assert model(center).item() > model(distal).item()


def test_train_save_load(tmp_path):
    from safetensors.torch import load_file

    torch.set_num_threads(1)
    data = tmp_path / "train.tsv"
    rows = []
    for split, n in [("train", 4), ("validation", 2), ("test", 2)]:
        for i in range(n):
            rows.append(
                dict(
                    sequence_id=f"{split}{i}",
                    sequence=("ACGT" if i % 2 else "AAAA") * 64,
                    split=split,
                    group_id=f"{split}{i}",
                    chrom=f"{split}{i}",
                    start=1000,
                    end=1064,
                    ATAC=float(i + 1),
                    H3K27ac=None if i == 0 else float(i + 2),
                )
            )
    write_table(data, rows)
    cfg = dict(
        data="train.tsv",
        model_id="test_cnn",
        species="synthetic",
        assembly="toy",
        context_id="test",
        target_level="individual",
        assays=["ATAC", "H3K27ac"],
        signal_unit="toy",
        normalization_id="toy",
        is_synthetic=True,
        input_length=256,
        output_window=64,
        epochs=2,
        batch_size=2,
        channels=[4, 4, 4],
    )
    config = tmp_path / "train.yaml"
    config.write_text(yaml.safe_dump(cfg))
    manifest = train_sequence(config, tmp_path / "model")
    report = json.loads((tmp_path / "model/training_report.json").read_text())
    assert report["parameters_updated"]
    assert all(
        np.isfinite(r["train_loss"]) and np.isfinite(r["validation_loss"])
        for r in report["history"]
    )
    # Training-only positive median excludes val/test and masked labels: ATAC=(1,2,3,4), H3=(3,4,5).
    assert manifest["scales"] == [2.5, 4.0]
    a, b = [QuantitativeCNN(256, 64, 2, channels=[4, 4, 4]) for _ in range(2)]
    state = load_file(tmp_path / "model/weights.safetensors")
    a.load_state_dict(state)
    b.load_state_dict(state)
    x = torch.ones(2, 4, 256)
    torch.testing.assert_close(a(x), b(x), rtol=0, atol=0)


def test_split_context_leak_rejected():
    rows = [
        dict(split="train", group_id="A", chrom="chr1", start=1000, end=1500),
        dict(split="test", group_id="B", chrom="chr1", start=2000, end=2500),
    ]
    with pytest.raises(PaceError, match="Overlapping sequence context"):
        validate_splits(rows, input_length=8192, output_window=500)
