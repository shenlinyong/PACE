# PACE implementation record

Updated 2026-09-17. Maintainer: 申林用 (Linyong Shen), Northwest A&F University.

This work implements the three supplied September 2026 contracts as an installable
Python package while preserving the existing region-based scripts and their history.
The existing public repository, `https://github.com/shenlinyong/PACE`, is the upload target.

## Delivery sequence

1. Strict schemas, pure numerical kernel, candidate catalogs and provenance.
2. Measured, hybrid and genome-only pipelines with offline synthetic examples.
3. Genomic adapters, conservative variant mapping, quantitative CNN and calibration.
4. Common-denominator comparison, grouped learning and benchmarks.
5. Regression tests, wheel installation, CI, bilingual documentation and push.

## Decisions and evidence boundaries

- The supplied documents override older formula defaults: eta=0, fixed assay panel,
  no arbitrary residual mass, and missing values distinct from biological zero.
- Existing `scripts/` and `workflow/` remain a separate legacy interface. Their outputs
  must not be compared directly with canonical-grid scores.
- Only synthetic fixtures are available. No livestock weights, biological accuracy,
  or validated individual effects are claimed.
- Skills: `research-software-engineering`, with resource and documentation companions,
  from a-attia/scicomp-research-skills at commit
  `8435b16d91972c4f31b006de7bcacf1f5eb47e8e`.
- Package templates follow Scientific Python's src-layout, PEP 621 metadata and
  wheel-test conventions. No copied scientific implementation is introduced.

## Status

Implementation and local verification completed. The canonical package, real fitting/adaptation
paths, offline examples, bilingual guides and CI are present. The local full suite passed 127
tests; wheel-only offline demos and actual training/inference commands passed. See
`docs/validation.md` for precise scope and environment, and `docs/limitations.md` for research
extensions and unavailable real biological assets. GitHub upload and hosted CI verification
are the remaining delivery steps.
