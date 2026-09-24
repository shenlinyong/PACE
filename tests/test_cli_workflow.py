"""Flag-only command line: one-step prediction, step-by-step commands and chunked runs."""

import json
import math

import numpy as np
import pytest

from pace_livestock.cli import main
from pace_livestock.io.tables import read_table

pyBigWig = pytest.importorskip("pyBigWig")
cooler = pytest.importorskip("cooler")
pd = pytest.importorskip("pandas")

CHROMS = {"chr1": 3_000_000, "chr2": 2_000_000, "chrX": 1_500_000}
RESOLUTION = 10_000


def make_inputs(root):
    """Small three-chromosome genome with peaks, genes, two bigWigs and a Hi-C map."""
    rng = np.random.default_rng(7)
    root.mkdir(parents=True, exist_ok=True)
    peaks, gtf = [], []
    for chrom, length in CHROMS.items():
        for start in sorted(rng.choice(np.arange(20_000, length - 20_000, 700), 60, False)):
            peaks.append(f"{chrom}\t{start}\t{start + 400}\n")
        for i, tss in enumerate(sorted(rng.choice(np.arange(50_000, length - 50_000, 5000), 8))):
            strand = "+" if i % 2 else "-"
            a, b = (tss + 1, tss + 5000) if strand == "+" else (tss - 5000, tss + 1)
            gene = f"{chrom}_G{i}"
            biotype = "lncRNA" if i == 0 else "protein_coding"
            attrs = f'gene_id "{gene}"; transcript_id "{gene}.1"; gene_biotype "{biotype}";'
            gtf.append(f"{chrom}\tsrc\ttranscript\t{a}\t{b}\t.\t{strand}\t.\t{attrs}\n")
            if i == 3:  # a second, distinct TSS for one gene
                attrs = attrs.replace(".1", ".2")
                gtf.append(f"{chrom}\tsrc\ttranscript\t{a + 800}\t{b}\t.\t{strand}\t.\t{attrs}\n")
    (root / "peaks.bed").write_text("".join(peaks))
    (root / "genes.gtf").write_text("".join(gtf))
    (root / "genome.fa.fai").write_text(
        "".join(f"{c}\t{n}\t0\t60\t61\n" for c, n in CHROMS.items())
    )
    for name in ("atac", "h3k27ac"):
        bw = pyBigWig.open(str(root / f"{name}.bw"), "w")
        bw.addHeader(list(CHROMS.items()))
        for chrom, length in CHROMS.items():
            values = rng.gamma(0.5, 1.0, length // 100).tolist()
            bw.addEntries(chrom, 0, values=values, span=100, step=100)
        bw.close()
    bins, pixels, offset = [], [], 0
    for chrom, length in CHROMS.items():
        n = math.ceil(length / RESOLUTION)
        bins.append(
            pd.DataFrame(
                {
                    "chrom": chrom,
                    "start": np.arange(n) * RESOLUTION,
                    "end": np.minimum((np.arange(n) + 1) * RESOLUTION, length),
                }
            )
        )
        for d in range(0, 150):
            i = np.arange(n - d)
            counts = rng.poisson(400 * max(d, 0.5) ** -1.0, len(i))
            keep = counts > 0
            pixels.append(
                pd.DataFrame(
                    {
                        "bin1_id": offset + i[keep],
                        "bin2_id": offset + i[keep] + d,
                        "count": counts[keep],
                    }
                )
            )
        offset += n
    pixels = pd.concat(pixels).sort_values(["bin1_id", "bin2_id"])
    cooler.create_cooler(str(root / "hic.cool"), pd.concat(bins), pixels, ordered=True)
    with _quiet():
        cooler.balance_cooler(cooler.Cooler(str(root / "hic.cool")), store=True)
    return root


class _quiet:
    def __enter__(self):
        import warnings

        self.context = warnings.catch_warnings()
        self.context.__enter__()
        warnings.simplefilter("ignore")

    def __exit__(self, *exc):
        self.context.__exit__(*exc)


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    return make_inputs(tmp_path_factory.mktemp("genome"))


def call(*argv):
    return main([str(a) for a in argv])


def scores(path):
    rows = read_table(path / "scores.tsv.gz")
    return {(r["element_id"], r["gene_id"]): r for r in rows}


def context():
    return ["--species", "pig", "--assembly", "toy", "--tissue", "liver"]


def test_predict_one_step_with_hic(inputs, tmp_path, capsys):
    out = tmp_path / "pred"
    code = call(
        "predict",
        "-b",
        inputs / "peaks.bed",
        "-g",
        inputs / "genes.gtf",
        "-c",
        inputs / "genome.fa.fai",
        "--atac",
        inputs / "atac.bw",
        "--h3k27ac",
        inputs / "h3k27ac.bw",
        "--hic",
        inputs / "hic.cool",
        "-r",
        1_000_000,
        "--chunk-pairs",
        700,
        *context(),
        "-o",
        out,
    )
    assert code == 0, capsys.readouterr().err
    summary = json.loads(capsys.readouterr().out)
    assert summary["contact_mode"] == "observed"
    qc = json.loads((out / "qc_report.json").read_text())
    assert qc["execution"]["mode"] == "by_chromosome"
    assert len(qc["execution"]["chunks"]) > 1
    assert (out / "prepared/contact_prior/manifest.json").is_file()
    samples = read_table(out / "inferred_samples.tsv")
    assert {r["assay"] for r in samples} == {"ATAC", "H3K27ac", "HiC"}
    # The result folder is self-contained and can be re-scored as one unchunked run.
    code = call("run", "--config", out / "resolved_config.yaml", "-o", tmp_path / "rerun")
    assert code == 0, capsys.readouterr().err
    chunked, single = scores(out), scores(tmp_path / "rerun")
    assert chunked.keys() == single.keys()
    for key, row in single.items():
        for field in ("pace_score", "support", "Cbar", "A_used"):
            a, b = row[field], chunked[key][field]
            assert (a is None and b is None) or float(a) == pytest.approx(float(b), rel=1e-12)
    genes = {}
    for r in single.values():
        if r["pace_score"] is not None:
            genes[r["gene_id"]] = genes.get(r["gene_id"], 0) + float(r["pace_score"])
    assert genes and all(v == pytest.approx(1) for v in genes.values())


def test_predict_needs_contact_choice(inputs, tmp_path, capsys):
    base = ["predict", "-b", inputs / "peaks.bed", "-g", inputs / "genes.gtf"]
    code = call(*base, "--atac", inputs / "atac.bw", *context(), "-o", tmp_path / "x")
    assert code == 2 and "--abc-prior" in capsys.readouterr().err
    assert not (tmp_path / "x").exists()
    code = call(*base, "--hic", inputs / "hic.cool", *context(), "-o", tmp_path / "y")
    assert code == 2 and "at least one activity bigWig" in capsys.readouterr().err


def test_predict_abc_prior_and_gene_types(inputs, tmp_path, capsys):
    out = tmp_path / "abc"
    code = call(
        "predict",
        "-b",
        inputs / "peaks.bed",
        "-g",
        inputs / "genes.gtf",
        "--atac",
        inputs / "atac.bw",
        "--abc-prior",
        "--gene-types",
        "protein_coding",
        "-r",
        500_000,
        *context(),
        "-o",
        out,
    )
    assert code == 0, capsys.readouterr().err
    genes = {r["gene_id"] for r in read_table(out / "gene_summary.tsv")}
    assert genes and not any(g.endswith("_G0") for g in genes)
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["asset_manifests"]["contact_prior"]["transfer_status"].startswith("unvalidated")


def test_step_by_step_commands_without_yaml(inputs, tmp_path, capsys):
    cat = tmp_path / "catalog"
    assert (
        call(
            "catalog",
            "-b",
            inputs / "peaks.bed",
            "-g",
            inputs / "genes.gtf",
            "-c",
            inputs / "atac.bw",
            "-r",
            1_000_000,
            "-o",
            cat,
        )
        == 0
    ), capsys.readouterr().err
    for assay, name in (("ATAC", "atac"), ("H3K27ac", "h3k27ac")):
        code = call(
            "activity", "-i", inputs / f"{name}.bw", "-a", assay, "-d", cat, "-o", tmp_path / name
        )
        assert code == 0, capsys.readouterr().err
    assert (
        call(
            "merge",
            "-t",
            "observed_activity",
            "-i",
            tmp_path / "atac/observed_activity.tsv",
            tmp_path / "h3k27ac/observed_activity.tsv",
            "-o",
            tmp_path / "activity",
        )
        == 0
    )
    code = call(
        "contacts", "-i", inputs / "hic.cool", "-r", RESOLUTION, "-d", cat, "-o", tmp_path / "hic"
    )
    assert code == 0, capsys.readouterr().err
    contacts = read_table(tmp_path / "hic/observed_contacts.tsv")
    assert {r["scale"] for r in contacts} == {"balanced_contact"}
    code = call(
        "fit-prior",
        "--cooler",
        inputs / "hic.cool",
        *context(),
        "--scale",
        "balanced_contact",
        "--normalization-id",
        "cooler_weight",
        "-o",
        tmp_path / "prior",
    )
    assert code == 0, capsys.readouterr().err
    capsys.readouterr()
    run = [
        "run",
        "-d",
        cat,
        "--activity",
        tmp_path / "activity/observed_activity.tsv",
        "--contacts",
        tmp_path / "hic/observed_contacts.tsv",
        "--contact-prior",
        tmp_path / "prior",
        "--allow-prior-fallback",
        *context(),
    ]
    # No samples/sources tables and no scale label: both are taken from the inputs.
    assert call(*run, "-o", tmp_path / "single") == 0, capsys.readouterr().err
    chunked_run = [*run, "--by-chromosome", "--chunk-pairs", 1]
    assert call(*chunked_run, "-o", tmp_path / "chunked") == 0
    assert call(*chunked_run, "--threads", 2, "-o", tmp_path / "parallel") == 0
    assert (tmp_path / "parallel/scores.tsv.gz").read_bytes() == (
        tmp_path / "chunked/scores.tsv.gz"
    ).read_bytes() or read_table(tmp_path / "parallel/scores.tsv.gz") == read_table(
        tmp_path / "chunked/scores.tsv.gz"
    )
    single, chunked = scores(tmp_path / "single"), scores(tmp_path / "chunked")
    assert single.keys() == chunked.keys()
    assert all(
        (single[k]["pace_score"] is None) == (chunked[k]["pace_score"] is None)
        and (
            single[k]["pace_score"] is None
            or float(single[k]["pace_score"]) == pytest.approx(float(chunked[k]["pace_score"]))
        )
        for k in single
    )
    manifest = json.loads((tmp_path / "chunked/run_manifest.json").read_text())
    assert manifest["execution"]["chunks"] == [["chr1"], ["chr2"], ["chrX"]]
    capsys.readouterr()
    code = call(
        "compare",
        "--left",
        tmp_path / "chunked",
        "--right",
        tmp_path / "chunked",
        "-o",
        tmp_path / "self",
    )
    assert code == 0, capsys.readouterr().err
    deltas = read_table(tmp_path / "self/comparison.tsv")
    assert all(float(r["full_delta_pace"]) == 0 for r in deltas if r["reason"] == "complete")


def test_catalog_rejects_unlisted_chromosomes_with_advice(inputs, tmp_path, capsys):
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t3000000\nchr2\t2000000\n")
    base = ["catalog", "-b", inputs / "peaks.bed", "-g", inputs / "genes.gtf", "-c", sizes]
    assert call(*base, "-o", tmp_path / "strict") == 2
    assert "--skip-unlisted-chroms" in capsys.readouterr().err
    assert call(*base, "--skip-unlisted-chroms", "-o", tmp_path / "loose") == 0
    units = read_table(tmp_path / "loose/units.tsv")
    assert {r["chrom"] for r in units} == {"chr1", "chr2"}
    report = json.loads((tmp_path / "loose/preparation_report.json").read_text())
    assert report["skipped_unlisted_chromosomes"]["chromosomes"] == ["chrX"]


def test_help_lists_commands_and_hides_yaml(capsys):
    assert main([]) == 0
    text = capsys.readouterr().out
    for command in ("predict", "catalog", "activity", "contacts", "run", "compare"):
        assert f"\n  {command} " in text
    with pytest.raises(SystemExit):
        main(["run", "--help"])
    assert "--config" not in capsys.readouterr().out


def script_commands(path):
    """pace commands of an example script, with its shell variables expanded."""
    import shlex

    variables, commands, pending = {}, [], ""
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("set "):
            continue
        if line.endswith("\\"):
            pending += line[:-1] + " "
            continue
        line, pending = pending + line, ""
        name, sep, value = line.partition("=")
        if sep and name.isidentifier() and name.isupper():
            variables[name] = shlex.split(value)[0]
            continue
        for key, value in variables.items():
            line = line.replace(f"${key}", value)
        words = shlex.split(line)
        assert words[0] == "pace", line
        commands.append(words[1:])
    return commands


@pytest.mark.parametrize("example", ["contact", "measured"])
def test_example_scripts_use_only_flags(tmp_path, example, capsys):
    import os
    import shutil
    from pathlib import Path

    source = Path(__file__).resolve().parents[1] / "examples" / example
    root = tmp_path / example
    shutil.copytree(source, root)
    commands = script_commands(root / "run.sh")
    assert commands and not any("--config" in c for c in commands)
    cwd = os.getcwd()
    os.chdir(root)
    try:
        for argv in commands:
            assert main(argv) == 0, (argv, capsys.readouterr().err)
    finally:
        os.chdir(cwd)
    if example == "contact":
        calibration = json.loads((root / "calibrated/eta_calibration.json").read_text())
        assert calibration["status"] == "weak_fitted"
        assert not calibration["functional_validation"]
