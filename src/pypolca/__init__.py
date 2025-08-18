#!/usr/bin/env python3
"""pypolca"""

from pathlib import Path

import click
from loguru import logger

from pypolca.utils.fix_consensus_from_vcf import fix_consensus_from_vcf
from pypolca.utils.mapping import (
    bam_to_sorted_bam,
    bwa_index,
    bwa_paired,
    bwa_single,
    sam_to_bam,
    samtools_index,
)
from pypolca.utils.report import create_report
from pypolca.utils.util import (
    begin_pypolca,
    copy_file,
    end_pypolca,
    get_version,
    print_citation,
    remove_directory,
)
from pypolca.utils.validation import (
    check_dependencies,
    check_memory_limit,
    instantiate_dirs,
    validate_fasta,
    validate_fastq,
)
from pypolca.utils.variants import run_freebayes, samtools_faidx

"""
some code adapted from tbpore https://github.com/mbhall88/tbpore
"""

log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "-a",
            "--assembly",
            help="Path to assembly contigs or scaffolds.",
            type=click.Path(),
            required=True,
        ),
        click.option(
            "-1",
            "--reads1",
            required=True,
            type=click.Path(),
            help="Path to polishing reads R1 FASTQ. Can be FASTQ or FASTQ gzipped. Required.",
        ),
        click.option(
            "-2",
            "--reads2",
            default=None,
            type=click.Path(),
            help="Path to polishing reads R2 FASTQ. Can be FASTQ or FASTQ gzipped. Optional. Only use -1 if you have single end reads.",
        ),
        click.option(
            "-t",
            "--threads",
            help="Number of threads.",
            default=1,
            show_default=True,
        ),
        click.option(
            "-o",
            "--output",
            default="output_pypolca",
            show_default=True,
            type=click.Path(),
            help="Output directory path",
        ),
        click.option(
            "-f", "--force", is_flag=True, help="Force overwrites the output directory"
        ),
        click.option(
            "--min_alt",
            help="Minimum alt allele count to make a change",
            default=2,
            type=int,
            show_default=True,
        ),
        click.option(
            "--min_ratio",
            help="Minimum alt allele to ref allele ratio to make a change",
            default=2.0,
            type=float,
            show_default=True,
        ),
        click.option(
            "--careful", is_flag=True, help="Equivalent to --min_alt 4 --min_ratio 3"
        ),
        click.option(
            "--homopolymers",
            help="Ignore all changes except for homopolymer-length changes, with homopolymers defined by this length",
            default=None,
            type=int,
            show_default=True,
        ),
        click.option(
            "-n",
            "--no_polish",
            is_flag=True,
            help="do not polish, just create vcf file, evaluate the assembly and exit",
        ),
        click.option(
            "-m",
            "--memory_limit",
            type=str,
            default="2G",
            help="Memory per thread to use in samtools sort, set to 2G or more for large genomes",
            show_default=True,
        ),
        click.option(
            "-p",
            "--prefix",
            type=str,
            default="pypolca",
            help="prefix for output files",
            show_default=True,
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
def main_cli():
    1 + 1


"""
run command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
def run(
    ctx,
    assembly,
    reads1,
    reads2,
    threads,
    output,
    force,
    min_alt,
    min_ratio,
    homopolymers,
    careful,
    no_polish,
    memory_limit,
    prefix,
    **kwargs,
):
    """Python implementation of the POLCA polisher from MaSuRCA"""

    if careful:
        min_alt, min_ratio = 4, 3
    params = {
        "--assembly": assembly,
        "--reads1": reads1,
        "--reads2": reads2,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--careful": careful,
        "--min_alt": min_alt,
        "--min_ratio": min_ratio,
        "--homopolymers": homopolymers,
        "--memory_limit": memory_limit,
        "--no_polish": no_polish,
        "--prefix": prefix,
    }

    # validates the directory  (need to before I start or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)

    temp_dir: Path = Path(output) / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    logdir: Path = Path(output) / "logs"
    logdir.mkdir(parents=True, exist_ok=True)

    # begin pypolca - initial logging etc
    start_time = begin_pypolca(params)

    # check dependencies
    check_dependencies()

    # validates fasta
    validate_fasta(assembly)

    # validating FASTQ
    # whether paired or single
    # paired_flag is True
    if reads2 is None:
        logger.warning("You have not specified -2 or --reads2.")
        logger.warning("pypolca will proceed using single end reads only.")
        paired_flag = False
        validate_fastq(reads1)
    else:
        paired_flag = True
        validate_fastq(reads1)
        validate_fastq(reads2)

    logger.info(f"Checking memory limit of {memory_limit}.")
    check_memory_limit(memory_limit)

    #########
    #
    logger.info("Creating BWA index")
    assembly_temp: Path = Path(temp_dir) / "assembly.fasta"

    copy_file(Path(assembly), assembly_temp)
    bwa_index(assembly_temp, logdir)

    logger.info("Aligning reads with BWA")
    sam: Path = temp_dir / "temp_bwa.sam"
    if paired_flag is True:
        bwa_paired(reads1, reads2, assembly_temp, sam, threads, logdir)
    else:
        bwa_single(reads1, assembly_temp, sam, threads, logdir)

    logger.info("Sorting and indexing alignment file")
    bam: Path = temp_dir / "temp_bwa.bam"
    sorted_bam: Path = temp_dir / "temp_bwa_sorted.bam"
    sam_to_bam(sam, bam, threads, logdir)
    bam_to_sorted_bam(bam, sorted_bam, threads, memory_limit, logdir)
    samtools_index(sorted_bam, logdir)

    logger.info("Calling variants.")
    samtools_faidx(assembly_temp, logdir)
    vcf: Path = Path(output) / f"{prefix}.vcf"
    run_freebayes(assembly_temp, sorted_bam, output, vcf, logdir)

    ####
    # fix_consensus_from_vcf

    out_fasta: Path = Path(output) / f"{prefix}_corrected.fasta"
    if no_polish is False:
        subs, indels = fix_consensus_from_vcf(
            assembly_temp, vcf, out_fasta, min_alt, min_ratio, homopolymers
        )
    else:
        subs, indels = 0, 0

    ######
    # report.py
    report_file: Path = Path(output) / f"{prefix}.report"
    create_report(vcf, assembly_temp, report_file, subs, indels)

    # finalise
    if no_polish is False:
        logger.info(
            f"Final report is in {report_file}; polished assembly is in {out_fasta}"
        )
    else:
        logger.info(f"Final report is in {report_file}")

    # cleanup temp
    remove_directory(temp_dir)

    # end pypolca
    end_pypolca(start_time)


@click.command()
def citation(**kwargs):
    """Print the citations for polca"""
    print_citation()


main_cli.add_command(citation)


def main():
    main_cli()


if __name__ == "__main__":
    main()
