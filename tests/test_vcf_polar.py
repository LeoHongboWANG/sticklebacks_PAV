"""Integration and regression tests for biological decisions and VCF integrity."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pysam
import pytest

import vcf_polar as vp

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="ALT count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Called copies">
##INFO=<ID=AF,Number=A,Type=Float,Description="ALT frequency">
##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequency">
##INFO=<ID=ANN,Number=.,Type=String,Description="Allele-dependent effect">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log likelihoods">
##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype probabilities">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phased genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="ALT dosage">
"""


def make_vcf(tmp_path, rows, samples=("s1", "s2"), name="input.vcf", header=HEADER):
    path = tmp_path / name
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if samples:
        columns += "\tFORMAT\t" + "\t".join(samples)
    path.write_text(header + columns + "\n" + "".join(row + "\n" for row in rows))
    return path


def read_records(path):
    with pysam.VariantFile(str(path)) as handle:
        return list(handle)


def run_file(tmp_path, input_path, *options, output="result.vcf"):
    path = tmp_path / output
    code = vp.main([str(input_path), "-o", str(path), *map(str, options)])
    return code, path


def test_default_annotation_preserves_reference_genotypes_and_allelic_fields(tmp_path):
    source = make_vcf(tmp_path, [
        "chr1\t1\tid1\tA\tG\t60\tPASS\tAA=G;AC=1;AN=4;AF=0.25;MAF=0.25;ANN=G|missense\tGT:AD:PL\t0|0:10,3:0,30,100\t0/1:5,5:40,0,40",
        "chr1\t2\t.\tC\tT\t.\tPASS\tAA=C\tGT:AD:PL\t1:1,8:80,0\t./1:2,4:.",
    ])
    code, output = run_file(tmp_path, source, "--report", tmp_path / "report.json")
    assert code == 0
    before, after = read_records(source), read_records(output)
    for old, new in zip(before, after):
        assert (new.contig, new.pos, new.id, new.alleles, new.qual, tuple(new.filter)) == (
            old.contig, old.pos, old.id, old.alleles, old.qual, tuple(old.filter)
        )
        for name in old.samples:
            for tag in ("GT", "AD", "PL"):
                assert new.samples[name][tag] == old.samples[name][tag]
            assert new.samples[name].phased == old.samples[name].phased
    assert after[0].info["AF"] == pytest.approx((0.25,))
    assert after[0].info["MAF"] == pytest.approx(0.25)
    assert after[0].info["ANN"] == ("G|missense",)
    assert after[0].samples["s1"]["VP_DOSAGE"] == 2
    assert after[0].samples["s2"]["VP_DOSAGE"] == 1
    assert after[0].info["VP_DAF"] == pytest.approx(0.75)
    assert after[1].samples["s1"]["VP_DOSAGE"] == 1
    assert after[1].samples["s2"]["VP_DOSAGE"] is None
    assert (after[1].info["VP_DAC"], after[1].info["VP_DAN"]) == (2, 2)
    report = json.loads((tmp_path / "report.json").read_text())
    assert report["resolved_records"] == report["output_records"] == 2
    assert report["swapped_records"] == 0


def test_outgroup_heterozygotes_are_evidence_and_no_implicit_info_fallback(tmp_path):
    names = ("fish", "o1", "o2", "o3", "o4")
    source = make_vcf(tmp_path, [
        "chr1\t1\t.\tA\tG\t.\tPASS\tAA=G\tGT:DP:GQ\t0/0:20:50\t1/1:10:50\t0/1:10:50\t0/1:10:50\t0/1:10:50",
    ], names)
    ancestors = tmp_path / "outgroups.txt"
    ancestors.write_text("o1\no2\no3\no4\n")
    code, output = run_file(tmp_path, source, "--anc-list", ancestors)
    assert code == 0
    rec = read_records(output)[0]
    assert rec.info["VP_STATUS"] == "LOW_SUPPORT"
    assert rec.info["VP_N_OUT"] == 4
    assert rec.info["VP_N_OUT_ALLELES"] == 8
    assert rec.info["VP_SUPPORT"] == pytest.approx(0.625)
    assert "AA" not in rec.info
    assert rec.info["VP_INPUT_AA"] == "G"
    assert "VP_DAF" not in rec.info
    code, output = run_file(tmp_path, source, "--anc-list", ancestors, "--fallback-info", output="fallback.vcf")
    assert code == 0
    rec = read_records(output)[0]
    assert rec.info["VP_STATUS"] == "PASS"
    assert rec.info["VP_SOURCE"] == "info_fallback"
    assert rec.info["VP_FALLBACK"] == "LOW_SUPPORT"
    assert rec.info["VP_DAC"] == 2
    assert rec.info["VP_DAN"] == 2


@pytest.mark.parametrize("out1,out2,extra,status,n", [
    ("1/1:10:50", "1/1:10:50", [], "PASS", 2),
    ("1:10:50", "1:10:50", [], "PASS", 2),
    ("1/1:2:50", "1/1:10:50", [], "PASS", 1),
    ("1/1:10:5", "1/1:10:50", ["--min-outgroups", "2"], "TOO_FEW_OUTGROUPS", 1),
    ("1/1:.:.", "1/1:.:.", [], "TOO_FEW_OUTGROUPS", 0),
    ("1/1:.:.", "1/1:.:.", ["--min-dp", "0", "--min-gq", "0"], "PASS", 2),
    ("./1:10:50", "1/1:10:50", ["--max-outgroup-missing", "0.4"], "TOO_MISSING_OUTGROUPS", 1),
    ("0/0:10:50", "1/1:10:50", [], "LOW_SUPPORT", 2),
])
def test_outgroup_quality_missingness_and_ploidy(tmp_path, out1, out2, extra, status, n):
    source = make_vcf(tmp_path, [
        f"chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t0/0:10:50\t{out1}\t{out2}"
    ], ("fish", "o1", "o2"))
    ancestors = tmp_path / "anc.txt"
    ancestors.write_text("o1\no2\n")
    code, output = run_file(tmp_path, source, "--anc-list", ancestors, *extra)
    assert code == 0
    rec = read_records(output)[0]
    assert rec.info["VP_STATUS"] == status
    assert rec.info["VP_N_OUT"] == n
    if status == "PASS":
        assert rec.info["VP_DAF"] == 1
        assert rec.info["VP_DAN"] == 2


def test_projected_fasta_handles_lowercase_missing_and_third_allele(tmp_path):
    fasta = tmp_path / "ancestor.fa"
    fasta.write_text(">chr1\nGgNC\n")
    pysam.faidx(str(fasta))
    source = make_vcf(tmp_path, [
        f"chr1\t{pos}\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0/0\t1/1"
        for pos in range(1, 6)
    ])
    code, output = run_file(tmp_path, source, "--ancestral-fasta", fasta)
    assert code == 0
    records = read_records(output)
    assert [r.info["VP_STATUS"] for r in records] == [
        "PASS", "LOW_CONFIDENCE_AA", "NO_AA", "AA_NOT_IN_ALLELES", "FASTA_MISSING",
    ]
    assert records[0].info["AA"] == "G"
    assert records[3].info["AA"] == "C"
    assert "VP_DAF" not in records[3].info
    code, output = run_file(tmp_path, source, "--ancestral-fasta", fasta, "--allow-lowercase", output="lower.vcf")
    assert code == 0
    assert read_records(output)[1].info["VP_STATUS"] == "PASS"


def test_unresolved_and_unsupported_records_retained_with_reasons(tmp_path):
    source = make_vcf(tmp_path, [
        "chr1\t1\t.\tA\tG\t.\tPASS\tAA=N\tGT\t0/0\t./.",
        "chr1\t2\t.\tA\tAT\t.\tPASS\tAA=A\tGT\t0/0\t0/1",
        "chr1\t3\t.\tA\tG,T\t.\tPASS\tAA=A\tGT\t1/2\t0/1",
        "chr1\t4\t.\tA\t.\t.\tPASS\tAA=A\tGT\t0/0\t0/0",
    ])
    report = tmp_path / "sites.tsv"
    code, output = run_file(tmp_path, source, "--drop-unresolved", "--sites-report", report)
    assert code == 0
    assert read_records(output) == []
    assert len(report.read_text().splitlines()) == 5
    code, output = run_file(tmp_path, source, output="retained.vcf")
    assert code == 0
    records = read_records(output)
    assert len(records) == 4
    assert [r.info["VP_STATUS"] for r in records] == ["NO_AA"] + ["UNSUPPORTED_VARIANT"] * 3
    assert records[2].samples["s1"]["GT"] == (1, 2)


def test_legacy_swap_updates_arrays_counts_and_phase_without_inverting_maf(tmp_path):
    source = make_vcf(tmp_path, [
        "chr1\t1\t.\tA\tG\t60\tPASS\tAA=G;AC=1;AN=4;AF=0.25;MAF=0.25;ANN=G|missense"
        "\tGT:DP:GQ:AD:PL:GL:GP:PGT:DS"
        "\t0|0:13:50:10,3:0,30,100:0,-3,-10:0.99,0.009,0.001:0|0:0.1"
        "\t0/1:10:50:5,5:40,0,40:-4,0,-4:0.01,0.98,0.01:0|1:1",
    ])
    report = tmp_path / "legacy.json"
    code, output = run_file(tmp_path, source, "--mode", "ref", "--report", report)
    assert code == 0
    rec = read_records(output)[0]
    assert rec.alleles == ("G", "A")
    assert rec.samples["s1"]["GT"] == (1, 1)
    assert rec.samples["s1"].phased
    assert rec.samples["s1"]["AD"] == (3, 10)
    assert rec.samples["s1"]["PL"] == (100, 30, 0)
    assert rec.samples["s1"]["GL"] == pytest.approx((-10, -3, 0))
    assert rec.samples["s1"]["GP"] == pytest.approx((0.001, 0.009, 0.99))
    assert rec.samples["s1"]["PGT"] == "1|1"
    assert rec.samples["s1"]["VP_DOSAGE"] == 2
    assert rec.info["AC"] == (3,)
    assert rec.info["AN"] == 4
    assert rec.info["AF"] == pytest.approx((0.75,))
    assert rec.info["MAF"] == pytest.approx(0.25)
    assert "ANN" not in rec.info
    assert "DS" not in rec.format
    assert set(rec.info["VP_DROPPED"]) == {"INFO_ANN", "FORMAT_DS"}
    counts = json.loads(report.read_text())
    assert counts["dropped_info_fields"] == {"ANN": 1}
    assert counts["dropped_format_fields"] == {"DS": 1}
    code, restored = run_file(tmp_path, output, "--mode", "alt", output="restored.vcf")
    assert code == 0
    old = read_records(source)[0]
    new = read_records(restored)[0]
    assert new.alleles == old.alleles
    for name in old.samples:
        for tag in ("GT", "AD", "PL"):
            assert new.samples[name][tag] == old.samples[name][tag]
        assert new.samples[name].phased == old.samples[name].phased


@pytest.mark.parametrize("gt,pl,expected_gt,expected_pl", [
    ("0", "0,60", (1,), (60, 0)),
    ("0/0/1", "0,10,20,30", (1, 1, 0), (30, 20, 10, 0)),
    ("./.", ".", (None, None), (None,)),
    (".|1", ".", (None, 0), (None,)),
])
def test_legacy_haploid_polyploid_and_missing_gt(tmp_path, gt, pl, expected_gt, expected_pl):
    source = make_vcf(tmp_path, [
        f"chr1\t1\t.\tA\tG\t.\tPASS\tAA=G\tGT:PL\t{gt}:{pl}",
    ], ("sample",))
    code, output = run_file(tmp_path, source, "--mode", "ref")
    assert code == 0
    rec = read_records(output)[0]
    assert rec.samples["sample"]["GT"] == expected_gt
    actual = rec.samples["sample"]["PL"]
    if expected_pl == (None,):
        assert all(v is None for v in actual)
    else:
        assert actual == expected_pl


def test_legacy_site_only_counts(tmp_path):
    source = make_vcf(tmp_path, [
        "chr1\t1\t.\tA\tG\t.\tPASS\tAA=G;AC=2;AN=10;AF=0.2;MAF=0.2",
    ], ())
    code, output = run_file(tmp_path, source, "--mode", "ref")
    assert code == 0
    rec = read_records(output)[0]
    assert rec.info["AC"] == (8,)
    assert rec.info["AF"] == pytest.approx((0.8,))
    assert rec.info["MAF"] == pytest.approx(0.2)


@pytest.mark.parametrize("extension", ["vcf.gz", "bcf"])
def test_compressed_and_bcf_outputs_have_working_indexes(tmp_path, extension):
    source = make_vcf(tmp_path, [
        "chr1\t1\t.\tA\tG\t.\tPASS\tAA=G\tGT\t0/0\t0/1",
        "chr1\t2\t.\tC\tT\t.\tPASS\tAA=C\tGT\t1/1\t0/1",
    ])
    code, output = run_file(tmp_path, source, "--index", "--threads", "2", output=f"result.{extension}")
    assert code == 0
    suffix = ".csi" if extension == "bcf" else ".tbi"
    assert Path(str(output) + suffix).exists()
    with pysam.VariantFile(str(output)) as handle:
        assert [r.pos for r in handle.fetch("chr1", 1, 2)] == [2]
    code, second = run_file(tmp_path, output, output="second.vcf")
    assert code == 0
    assert [r.pos for r in read_records(second)] == [1, 2]


def test_original_entry_point_and_standard_streams(tmp_path):
    source = make_vcf(tmp_path, [
        "chr1\t1\t.\tA\tG\t.\tPASS\tAA=G\tGT\t0|0\t0/1",
    ])
    script = Path(vp.__file__).parent / "02.deleterious_mutation" / "VCF-Polar.py"
    result = subprocess.run(
        [sys.executable, str(script), "-", "--by-chrom", "--threads", "1"],
        input=source.read_text(), text=True, capture_output=True,
        env=os.environ.copy(),
    )
    assert result.returncode == 0, result.stderr
    assert "deprecated" in result.stderr
    assert "daemonic" not in result.stderr
    assert "VP_DOSAGE" in result.stdout
    assert result.stdout.splitlines()[-1].startswith("chr1\t1\t.\tA\tG")
    assert "Processed" not in result.stdout


@pytest.mark.parametrize("options", [
    ["--consensus-threshold", "nan"], ["--consensus-threshold", "0.5"],
    ["--consensus-threshold", "1.1"], ["--min-dp", "-1"],
    ["--min-gq", "nan"], ["--max-outgroup-missing", "inf"],
    ["--threads", "0"], ["--min-outgroups", "0"],
    ["--index"], ["--source", "outgroup"], ["--source", "fasta"],
])
def test_invalid_parameters_fail_without_output(tmp_path, options):
    source = make_vcf(tmp_path, ["chr1\t1\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0/0\t0/1"])
    code, output = run_file(tmp_path, source, *options)
    assert code == 2
    assert not output.exists()


def test_sample_names_and_overlap_are_validated(tmp_path):
    source = make_vcf(tmp_path, ["chr1\t1\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0/0\t0/1"])
    names = tmp_path / "names.txt"
    for content, extra in [
        ("absent\n", []), ("s1\ns1\n", []), ("", []),
        ("s1\n", ["--ingroup-list", str(names)]),
    ]:
        names.write_text(content)
        code, output = run_file(tmp_path, source, "--anc-list", names, *extra)
        assert code == 2
        assert not output.exists()


def test_input_and_reports_cannot_overwrite_each_other(tmp_path):
    source = make_vcf(tmp_path, ["chr1\t1\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0/0\t0/1"])
    before = source.read_bytes()
    assert vp.main([str(source), "-o", str(source)]) == 2
    assert source.read_bytes() == before
    output = tmp_path / "result.vcf"
    assert vp.main([str(source), "-o", str(output), "--report", str(output)]) == 2
    assert not output.exists()


def test_index_failure_preserves_existing_output_and_cleans_temporary_files(tmp_path):
    source = make_vcf(tmp_path, [
        "chr1\t2\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0/0\t0/1",
        "chr1\t1\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0/0\t0/1",
    ])
    output = tmp_path / "result.vcf.gz"
    output.write_bytes(b"previous result")
    code, _ = run_file(tmp_path, source, "--index", output=output.name)
    assert code == 2
    assert output.read_bytes() == b"previous result"
    assert not list(tmp_path.glob(".vcf-polar-*"))


def test_empty_input_has_a_valid_empty_output(tmp_path):
    source = make_vcf(tmp_path, [])
    code, output = run_file(tmp_path, source, "--by-chrom", "--report", tmp_path / "empty.json")
    assert code == 0
    assert read_records(output) == []
    assert json.loads((tmp_path / "empty.json").read_text())["input_records"] == 0


def test_ordinary_gzip_input_from_v1_is_readable(tmp_path):
    import gzip
    source = make_vcf(tmp_path, ["chr1\t1\t.\tA\tG\t.\tPASS\tAA=G\tGT\t0/0\t0/1"])
    compressed = tmp_path / "old.vcf.gz"
    with gzip.open(compressed, "wt") as handle:
        handle.write(source.read_text())
    code, output = run_file(tmp_path, compressed)
    assert code == 0
    assert read_records(output)[0].info["VP_DAF"] == pytest.approx(0.75)


def test_external_evidence_does_not_require_existing_aa_header(tmp_path):
    source = make_vcf(
        tmp_path, ["chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t1/1:10:50"],
        ("out",), header=HEADER.replace(
            '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">\n', ""
        ),
    )
    names = tmp_path / "anc.txt"
    names.write_text("out\n")
    code, output = run_file(tmp_path, source, "--anc-list", names)
    assert code == 0
    rec = read_records(output)[0]
    assert rec.info["AA"] == "G"
    assert rec.info["VP_DAN"] == 0
    assert "VP_DAF" not in rec.info


def test_index_paths_cannot_replace_an_input_or_report(tmp_path):
    source = make_vcf(tmp_path, ["chr1\t1\t.\tA\tG\t.\tPASS\tAA=G\tGT\t0/0\t0/1"])
    output = tmp_path / "result.vcf.gz"
    assert vp.main([str(source), "-o", str(output), "--report", str(output) + ".tbi"]) == 2
    assert not output.exists()


def test_legacy_mixed_phasing_fails_without_replacing_an_existing_result(tmp_path):
    source = make_vcf(tmp_path, [
        "chr1\t1\t.\tA\tG\t.\tPASS\tAA=G\tGT\t0|1/1"
    ], ("sample",))
    output = tmp_path / "existing.vcf"
    output.write_text("previous result\n")
    code, _ = run_file(tmp_path, source, "--mode", "ref", output=output.name)
    assert code == 2
    assert output.read_text() == "previous result\n"
    code, annotated = run_file(tmp_path, source, output="annotated.vcf")
    assert code == 0
    assert annotated.read_text().splitlines()[-1].split("\t")[9].split(":")[0] == "0|1/1"


def test_bundled_example_has_documented_outgroup_decisions(tmp_path):
    examples = Path(vp.__file__).parent / "02.deleterious_mutation" / "examples"
    code, output = run_file(tmp_path, examples / "tiny.vcf", "--anc-list", examples / "outgroups.txt")
    assert code == 0
    assert [r.info["VP_STATUS"] for r in read_records(output)] == [
        "PASS", "PASS", "LOW_SUPPORT", "TOO_FEW_OUTGROUPS", "UNSUPPORTED_VARIANT", "PASS",
    ]
