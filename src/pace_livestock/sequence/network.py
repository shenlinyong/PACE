"""Reference CNN with central-target pooling and separately pooled context."""

import torch
from torch import nn

from ..errors import PaceError


class QuantitativeCNN(nn.Module):
    def __init__(self, input_length=8192, output_window=500, heads=2, *, channels=(64, 128, 128)):
        super().__init__()
        if (
            input_length % 64
            or output_window > input_length
            or output_window <= 0
            or heads < 1
            or len(channels) != 3
        ):
            raise PaceError(
                "CNN needs input_length divisible by 64, positive central window and heads"
            )
        layers, previous = [], 4
        for width, kernel in zip(channels, (15, 7, 5), strict=True):
            layers += [
                nn.Conv1d(previous, width, kernel, padding=kernel // 2),
                nn.GELU(),
                nn.AvgPool1d(4),
            ]
            previous = width
        self.encoder = nn.Sequential(*layers)
        left = (input_length - output_window) / 2
        starts = torch.arange(input_length // 64, dtype=torch.float32) * 64
        overlap = (
            torch.minimum(starts + 64, torch.tensor(left + output_window))
            - torch.maximum(starts, torch.tensor(left))
        ).clamp(min=0)
        self.register_buffer("target_weights", overlap / overlap.sum())
        self.head = nn.Sequential(nn.Linear(previous * 2, heads), nn.Softplus())

    def forward(self, x):
        encoded = self.encoder(x)
        local = (encoded * self.target_weights).sum(dim=-1)
        context = encoded.mean(dim=-1)
        return self.head(torch.cat([local, context], dim=1))


def masked_huber(predicted, labels, scales):
    mask = torch.isfinite(labels)
    if not mask.any():
        raise PaceError("Training batch has no finite labels")
    safe = torch.where(mask, labels, torch.zeros_like(labels))
    loss = nn.functional.huber_loss(
        torch.log1p(predicted / scales), torch.log1p(safe / scales), reduction="none"
    )
    return loss[mask].mean()
