# Migration to the single current model

The current branch distributes one PACE scoring implementation under
`src/pace_livestock`. Formula, file schemas, tutorials and command aliases all refer
to it. Superseded region workflows, their configurations, results and model-specific
tests were removed from this branch because retaining them alongside the current
model led to incorrect public instructions.

The previous source tree is available at
[commit 8b5d12e](https://github.com/shenlinyong/PACE/tree/8b5d12e8f3ad3aa6948ee26b4ac099bde733579d).
That is an archival snapshot, not current usage or validation guidance. Git history
was preserved; historical outputs were not recomputed or relabelled.

Install the current package, prepare inputs with the [current schemas](data_dictionary.md),
and use `PACE measured`, `PACE hybrid` or `PACE genome`. `scripts/pace.py` is only a
compatibility launcher for this same installed package; it has no independent
scientific implementation. Old options are rejected rather than silently translated.

Existing raw inputs may be reusable after checking units, windows, fixed candidate
sets and evidence provenance. Old normalized scores are not interchangeable with
current results. Recompute from compatible evidence and evaluate changes with an
explicit [comparison contract](comparison.md).
