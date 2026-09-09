"""
Unit tests for pypolca's mapping.py file

These check the exact command lines built for the external tools. This matters
because a malformed samtools command does not always fail loudly - e.g.
`samtools index -@ <bam>` makes samtools swallow the BAM path as the argument
to -@, leaving no input file. Samtools >=1.16 then prints its usage to stdout
and exits 0, so pypolca happily continues without ever writing the .bai index,
while older samtools exits 1 and aborts the run.
"""

from pathlib import Path

import pytest

import pypolca.utils.mapping as mapping

# samtools/bwa options that must be followed by a value
OPTS_NEEDING_A_VALUE = {"-@", "-t", "-m", "-o", "-f", "-v", "-b"}


@pytest.fixture
def captured(monkeypatch):
    """Capture the ExternalTool passed to run_tool instead of running it."""
    calls = []

    def fake_run_tool(tool, ctx=None, to_stdout=False):
        calls.append(tool)

    monkeypatch.setattr(mapping.ExternalTool, "run_tool", fake_run_tool)
    return calls


def assert_options_have_values(command):
    """Every value-taking option must be followed by a non-option token."""
    for i, token in enumerate(command):
        if token in OPTS_NEEDING_A_VALUE:
            assert i + 1 < len(command), f"{token} is the last token in {command}"
            value = command[i + 1]
            assert not value.startswith("-"), f"{token} has no value in {command}"


def test_samtools_index_passes_threads_and_bam(captured, tmp_path):
    """
    Regression test for https://github.com/gbouras13/pypolca/issues/28

    -@ must be given a thread count, and the BAM must survive as a positional
    argument. Previously this built `samtools index -@ <bam>`.
    """
    bam = tmp_path / "temp_bwa_sorted.bam"
    mapping.samtools_index(bam, 8, tmp_path / "logs")

    (tool,) = captured
    assert tool.command == ["samtools", "index", "-@", "8", str(bam)]
    assert_options_have_values(tool.command)


def test_sam_to_bam_command(captured, tmp_path):
    sam = tmp_path / "temp_bwa.sam"
    bam = tmp_path / "temp_bwa.bam"
    mapping.sam_to_bam(sam, bam, 8, tmp_path / "logs")

    (tool,) = captured
    assert tool.command == [
        "samtools",
        "view",
        "-h",
        "-@",
        "8",
        "-b",
        str(sam),
        "-o",
        str(bam),
    ]


def test_bam_to_sorted_bam_command(captured, tmp_path):
    bam = tmp_path / "temp_bwa.bam"
    sorted_bam = tmp_path / "temp_bwa_sorted.bam"
    mapping.bam_to_sorted_bam(bam, sorted_bam, 8, "1G", tmp_path / "logs")

    (tool,) = captured
    assert tool.command == [
        "samtools",
        "sort",
        "-m",
        "1G",
        "-@",
        "8",
        str(bam),
        "-o",
        str(sorted_bam),
    ]
    assert_options_have_values(tool.command)


def test_bwa_paired_command(captured, tmp_path):
    genome = tmp_path / "assembly.fasta"
    r1, r2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
    mapping.bwa_paired(r1, r2, genome, tmp_path / "out.sam", 8, tmp_path / "logs")

    (tool,) = captured
    assert tool.command == [
        "bwa",
        "mem",
        "-SP",
        "-t",
        "8",
        str(genome),
        str(r1),
        str(r2),
    ]
    assert_options_have_values(tool.command)


def test_bwa_single_command(captured, tmp_path):
    genome = tmp_path / "assembly.fasta"
    r1 = tmp_path / "r1.fastq.gz"
    mapping.bwa_single(r1, genome, tmp_path / "out.sam", 8, tmp_path / "logs")

    (tool,) = captured
    assert tool.command == ["bwa", "mem", "-SP", "-t", "8", str(genome), str(r1)]
    assert_options_have_values(tool.command)
