"""Bind reusable quantitative predictions to the exact prepared genome inputs."""

from __future__ import annotations

from ..errors import PaceError
from ..provenance import digest, file_hash


def genome_binding_id(cfg, units, asset):
    """Return a reproducible input binding, never an assertion of model accuracy.

    File contents, target coordinates, sample, policies and model identity are
    included. Paths and run IDs are excluded so exports remain portable.
    """
    genome = cfg["genome"]
    if genome["variant_path"] and not genome["reference_path"]:
        raise PaceError("Individual variants require genome.reference_path, including imports")
    return digest(
        {
            "binding_schema": "pace-genome-binding-1",
            "context": cfg["context"],
            "target_level": cfg["target_level"],
            "genome": {
                key: file_hash(value) if value and key.endswith("_path") else value
                for key, value in genome.items()
            },
            "model_manifest_sha256": asset["manifest_sha256"],
            "input_length": asset["input_length"],
            "output_window": asset["output_window"],
            "max_n_fraction": cfg["sequence"]["max_n_fraction"],
            "units": sorted(
                (u["element_id"], u["chrom"], u["start"], u["end"], u["anchor0"]) for u in units
            ),
        }
    )


def bind_predictions(rows, binding_id):
    """Attach the binding when predictions are produced from validated windows."""
    return [{**row, "genome_binding_id": binding_id} for row in rows]


def validate_import_binding(predictions, resolved_activity, binding_id):
    """Reject imports from a different or undocumented prepared genome input."""
    rows = list(predictions) + [
        r for r in resolved_activity if r["evidence_type"] in ("sequence_prediction", "fused")
    ]
    for row in rows:
        if row.get("genome_binding_id") != binding_id:
            raise PaceError(
                "Imported sequence activity requires a matching genome_binding_id; "
                "export with pace predict-sequence or pace run using the same reference, "
                "VCF sample, callability, ploidy, targets and model"
            )
