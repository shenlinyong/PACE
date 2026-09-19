# PACE tutorial: from installation to a documented analysis

This tutorial uses complete offline fixtures first, then explains which files and
parameters change for a real study. A longer [Chinese manual](USER_GUIDE.zh-CN.md)
contains three fully commented project templates and track-preparation scripts.

## Start from the evidence, not from the parameter list

Use `measured` when you have qualified ATAC-seq/DNase-seq/H3K27ac and promoter-contact measurements. Use `hybrid` when those measurements are supplemented by a matching sequence model. Use `genome` when the main evidence is a reference or individual genome plus a validated sequence model and contact prior. The number of biological replicates may be one, two, three or more; each replicate is a row in `samples.tsv` and is referenced by `sample_id` in the measurement tables.

RNA-seq, histone marks other than H3K27ac, CTCF and WGBS/RRBS can be supplied as named annotations. They are not silently added to the primary score. See [MULTIOMICS.md](MULTIOMICS.md) before deciding whether a layer belongs in `observed_activity.tsv`, `features.tsv`, `expression.tsv` or `methylation.tsv`.

## 1. Install

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
conda env create -f environment.yml
conda activate pace
PACE --version
```

See [installation](INSTALLATION.md) for wget, pip/venv, optional CPU sequence
training and Docker. The default Conda/Docker setup includes `io,ml`; the offline
fixture sequence models do not need PyTorch.

## 2. Run and inspect measured evidence

```bash
PACE validate --config examples/measured/config.yaml
PACE measured --config examples/measured/config.yaml --out results/tutorial_measured
```

The complete runnable configuration is
[examples/measured/config.yaml](../examples/measured/config.yaml). It declares
synthetic species/assembly/context, one individual, fixed ATAC+H3K27ac activity,
observed contact, the candidate/promoter/unit tables and actual sample/source
metadata. Its tiny hand-calculation catalog explicitly excludes promoter units;
the main research default includes them.

```yaml
# Key measured-mode settings; merge into a full config, not a standalone file.
regime: measured
activity:
  panel: [ATAC, H3K27ac]
contact:
  mode: observed
allocation:
  eta: auto
```

For one required layer such as H3K27ac, declare `[H3K27ac]` for the whole run.
A two-layer run does not silently switch a missing element to one layer. Missing
and invalid observations remain unavailable; a real zero remains zero.

Inspect a few rows without extra software:

```bash
python - <<'PY'
import csv, gzip
from itertools import islice
with gzip.open('results/tutorial_measured/scores.tsv.gz', 'rt') as handle:
    for row in islice(csv.DictReader(handle, delimiter='\t'), 5):
        print({k: row.get(k) for k in
               ['element_id','gene_id','pace_score','A_used','Cbar','normalization_status','reason']})
PY
```

Check `gene_summary.tsv` and `qc_report.json` before treating a high score as
informative. Complete and partial backgrounds can both sum to one.

## 3. Run hybrid evidence

```bash
PACE capabilities --config examples/hybrid/config.yaml
PACE hybrid --config examples/hybrid/config.yaml --out results/tutorial_hybrid
```

The complete configuration is
[examples/hybrid/config.yaml](../examples/hybrid/config.yaml). In addition to
observations, it provides a sequence model, reference, VCF, callable regions,
chromosome ploidy and an assay-matched fusion calibrator. The observed sample and
individual genotype identify the same donor. It uses observed contact; hybrid
activity does not require contact shrinkage.

```yaml
# Hybrid-specific additions; paths here illustrate a project layout.
regime: hybrid
sequence:
  model_path: models/sequence
fusion:
  calibrator_path: models/fusion
  quality_stratum: default
genome:
  reference_path: data/reference.fa
  variant_path: data/animal1.vcf.gz
  individual_id: animal1
  sample_id: animal1
  callability_path: data/callable.bed
  ploidy_path: data/ploidy.tsv
```

For real inference install `.[sequence]` when the asset is the implemented CNN.
The calibrator combines matching assay signals in a frozen log1p space. Without
an identifiable applicable calibrator, evidence resolution prioritizes a qualified
observation and otherwise a valid prediction; it does not invent a fusion weight.
To shrink contact toward a prior, explicitly supply `contact.mode: shrinkage`,
`prior_path`, a measured/justified `reliability` and `reliability_source`.

## 4. Run genome-only evidence

```bash
PACE capabilities --config examples/genome_only/config.yaml
PACE genome --config examples/genome_only/config.yaml --out results/tutorial_genome
```

The full configuration is
[examples/genome_only/config.yaml](../examples/genome_only/config.yaml). Unlike
hybrid, it contains no observed activity/contact. Contact is an explicit prior:

```yaml
regime: genome_only
contact:
  mode: prior_only
  prior_path: models/contact
  allow_prior_fallback: true
sequence:
  model_path: models/sequence
```

It also requires the same canonical catalog and, for individual prediction, the
reference/genotype/callability/ploidy block shown above. Pure reference prediction
omits the genotype and its associated individual settings. A valid species/tissue
model is required in both cases. The supplied weights are synthetic and cannot be
used to make real livestock predictions.

Individual imports must retain matching genomic inputs and their
`genome_binding_id`; see [input preparation](input_preparation.md). Sequence-only
output means tissue-conditioned regulatory potential, not measured chromatin or
an expression effect. Unsupported structural geometry remains unavailable.

## 5. Prepare a real project

Use [input preparation](input_preparation.md) to create:

1. A fixed catalog from BED, GTF and chromosome lengths: units, promoters,
   candidates and region mappings. Include the exported `chrom_sizes.tsv` and
   candidate radius in the run contract.
2. Per-sample quantitative activity from normalized bigWig tracks. Match fixed
   windows and preserve zero versus missing values.
3. Per-sample contact from cool/mcool at one declared resolution and measurement
   scale. Do not mix normalization or balancing protocols by renaming fields.
4. Actual sample/source tables. An individual run cannot average different animals.
5. Optional [RNA, methylation and other omics](MULTIOMICS.md).

Set `execution_profile: research`, the real species, exact reference assembly and
context identifier. Validated mode requires the applicable asset reports and
checksums; selecting it does not itself validate a model. Model manifests bind
outputs to context, target level, assay, quantitative units and target window.

## 6. Key options shared and specific to modes

| Option | Meaning | Mode |
|---|---|---|
| `inputs.units/promoters/candidates` | Fixed candidate background and physical TSSs | All |
| `catalog.chrom_sizes_path` | Reference bounds for boundary-window checks | All canonical catalogs |
| `catalog.candidate_radius_bp` | Declared candidate generation distance | All; must match preparation |
| `target_level` | individual or population_mean | All |
| `inputs.samples/sources` | Real observations and their provenance | Measured/hybrid, optional annotations |
| `activity.panel` | Fixed one/two-assay activity definition | All |
| `activity.minimum_callable_fraction` | Observation usability threshold | Measured/hybrid |
| `contact.resolution/scale/normalization_id/balancing/window_id` | Complete contact measurement identity | All contact evidence |
| `contact.mode` | observed, prior_only or shrinkage | genome_only requires prior_only |
| `contact.reliability/reliability_source` | Explicit observed-contact weight and its origin | Shrinkage |
| `sequence.model_path` | Matching quantitative weights/manifest | Prediction modes |
| `fusion.calibrator_path` | Matching log-space fusion asset | Hybrid fusion |
| `genome.reference_path` | Same-assembly FASTA | Local prediction/genomic reconstruction |
| `genome.variant_path/sample_id/individual_id` | Actual individual's variants and identity | Individual prediction |
| `genome.callability_path/ploidy_path` | Reference-call evidence and chromosome copy number | Individual prediction |
| `allocation.eta` | Automatic default, or predeclared research number in [0,1] | All |
| `multiomics.mode/model_path` | Annotation or compatible independent classifier | All |

The [complete defaults and constraints](parameters.md) are generated from the
installed configuration. Direct flags, including `--contact-resolution` and
`--reference-cpg`, are shown by `PACE run --help` and the [CLI reference](cli.md).

## 7. Calibrate, compare and preserve provenance

For applicable functional data:

```bash
PACE fit-eta --config my_run.yaml --eta-labels functional_labels.tsv --out results/calibration
PACE measured --config another_animal.yaml \
  --eta-model results/calibration/eta_calibration.json --out results/another_animal
```

A fit is not automatically approved for deployment. Independent connected groups,
out-of-fold ranking, conservative selection and stability checks may legitimately
leave eta at zero. Test outcomes never tune it. See [calibration](eta_calibration.md).

Compare compatible runs using `PACE compare`, with a YAML giving left/right run
paths, rather than subtracting unmatched normalized scores. Complete and
conditional differences are separate. See [comparisons](comparison.md).

Save the exact commit, resolved configuration, manifest, QC, model identities and
source records with every research analysis. The [formula manual](FORMULA.md)
explains what the output measures; [limitations](limitations.md) lists what the
implementation and synthetic tests do not establish.
