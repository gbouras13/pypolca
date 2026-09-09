from pathlib import Path

from pypolca.utils.external_tools import ExternalTool


def bwa_index(genome: Path, logdir: Path) -> None:
    """
    Runs bwa index

    Parameters:
        genome (Path): Path to the reference genome to index

    Returns:
        None

    """

    bwa_index = ExternalTool(
        tool="bwa",
        input=f"",
        output=f"",
        params=f"index {genome} ",
        logdir=logdir,
        outfile="",
    )

    ExternalTool.run_tool(bwa_index)


def bwa_paired(
    reads1: Path, reads2: Path, genome: Path, out_sam: Path, threads: int, logdir: Path
) -> None:
    """
    Align paired-end sequencing reads to a reference genome using bwa mem.

    Parameters:
        reads1 (Path): Path to the first set of input FASTQ files containing paired-end reads.
        reads2 (Path): Path to the second set of input FASTQ files containing paired-end reads.
        genome (Path): Path to the reference genome to align the reads against.
        out_sam (Path): Path to the output SAM where the alignment results will be stored.
        threads (int): Number of threads/cores to use for the alignment.
        logdir (Path): Path to the directory where log files will be stored.

    Returns:
        None

    This function uses bwa to perform alignment of paired-end sequencing reads ('reads1' and 'reads2')
    against a reference genome ('genome'). The aligned output is stored in a SAM (Sequence Alignment/Map)
    format file named "temp_bwa.sam" in the 'output' directory.

    The 'threads' parameter specifies the number of CPU cores to use for the alignment, and 'logdir' is the
    directory where log files for the alignment process will be stored.
    """

    bwa_paired = ExternalTool(
        tool="bwa",
        input=f"",
        output=f"",
        params=f"mem -SP -t {str(threads)} {genome} {reads1} {reads2} ",
        logdir=logdir,
        outfile=out_sam,
    )

    ExternalTool.run_tool(bwa_paired, to_stdout=True)


def bwa_single(
    reads1: Path, genome: Path, out_sam: Path, threads: int, logdir: Path
) -> None:
    """
    Align single-end sequencing reads to a reference genome using bwa mem.

    Parameters:
        reads1 (Path): Path to the first set of input FASTQ files containing paired-end reads.
        genome (Path): Path to the reference genome to align the reads against.
        out_sam (Path): Path to the output SAM where the alignment results will be stored.
        threads (int): Number of threads/cores to use for the alignment.
        logdir (Path): Path to the directory where log files will be stored.

    Returns:
        None

    This function uses bwa to perform alignment of paired-end sequencing reads ('reads1' and 'reads2')
    against a reference genome ('genome'). The aligned output is stored in a SAM (Sequence Alignment/Map)
    format file named "temp_bwa.sam" in the 'output' directory.

    The 'threads' parameter specifies the number of CPU cores to use for the alignment, and 'logdir' is the
    directory where log files for the alignment process will be stored.
    """

    bwa_single = ExternalTool(
        tool="bwa",
        input=f"",
        output=f"",
        params=f"mem -SP -t {str(threads)} {genome} {reads1} ",
        logdir=logdir,
        outfile=out_sam,
    )

    ExternalTool.run_tool(bwa_single, to_stdout=True)


def sam_to_bam(sam: Path, bam: Path, threads: int, logdir: Path) -> None:
    """converts sam to bam with samtools
    :param outdir: output directory path
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    sam_to_bam = ExternalTool(
        tool="samtools",
        input="",
        output="",
        params=f" view -h -@ {threads} -b {sam} -o {bam}",
        logdir=logdir,
        outfile="",
    )

    # need to write to stdout
    ExternalTool.run_tool(sam_to_bam, to_stdout=False)


def bam_to_sorted_bam(
    bam: Path, sorted_bam: Path, threads: int, memory_limit: str, logdir: Path
) -> None:
    """converts bam to sorted bam with samtools
    :param outdir: output directory path
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    bam_to_sorted_bam = ExternalTool(
        tool="samtools",
        input="",
        output="",
        params=f" sort -m {memory_limit} -@ {threads} {bam} -o {sorted_bam}",
        logdir=logdir,
        outfile="",
    )

    # need to write to stdout
    ExternalTool.run_tool(bam_to_sorted_bam, to_stdout=False)


def samtools_index(sorted_bam: Path, threads: int, logdir: Path) -> None:
    """indexes a sorted bam with samtools
    :param sorted_bam: sorted bam path
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    samtools_index = ExternalTool(
        tool="samtools",
        input="",
        output="",
        params=f" index -@ {threads} {sorted_bam} ",
        logdir=logdir,
        outfile="",
    )

    # need to write to stdout
    ExternalTool.run_tool(samtools_index, to_stdout=False)


def samtools_cat(bams: list, merged_bam: Path, logdir: Path) -> None:
    """concatenates bams with samtools cat

    The input BAMs are all aligned against the same bwa index, so their @SQ
    headers match and they can be concatenated without being sorted first.

    No -@ here on purpose. samtools cat only takes the long form --threads
    (and not at all in older versions), and it copies BGZF blocks without
    recompressing them, so threading buys effectively nothing.

    :param bams: list of bam paths to concatenate
    :param merged_bam: output merged bam path
    :param logdir: logdir
    :return:
    """

    bam_str = " ".join(str(bam) for bam in bams)

    samtools_cat = ExternalTool(
        tool="samtools",
        input="",
        output="",
        params=f" cat -o {merged_bam} {bam_str} ",
        logdir=logdir,
        outfile="",
    )

    # need to write to stdout
    ExternalTool.run_tool(samtools_cat, to_stdout=False)
