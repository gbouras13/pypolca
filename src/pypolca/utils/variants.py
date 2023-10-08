from pathlib import Path

from pypolca.utils.external_tools import ExternalTool


def samtools_faidx(genome: Path, logdir: Path) -> None:
    samtools_faidx = ExternalTool(
        tool="samtools",
        input="",
        output="",
        params=f" faidx {genome} ",
        logdir=logdir,
        outfile="",
    )

    # need to write to stdout
    ExternalTool.run_tool(samtools_faidx, to_stdout=False)


def run_freebayes(
    genome: Path, sorted_bam: Path, output: Path, vcf: Path, logdir: Path
) -> None:
    """
    Run FreeBayes variant calling tool on input data.

    Args:
        genome (Path): Path to the reference genome file.
        sorted_bam (Path): Path to the sorted BAM file containing aligned sequencing reads.
        output (Path): Path to the directory where the output files will be stored.
        vcf (Path): vcf Path
        logdir (Path): Path to the directory where log files will be stored.

    Returns:
        None
    """

    freebayes = ExternalTool(
        tool="freebayes",
        input="",
        output="",
        params=f"  -m 0 --min-coverage 3 -R 0 -p 1 -F 0.2 -E 0 -b {sorted_bam}  -f {genome} -v {vcf} ",
        logdir=logdir,
        outfile="",
    )

    # need to write to stdout
    ExternalTool.run_tool(freebayes, to_stdout=False)
