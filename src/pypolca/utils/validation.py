import gzip
import os
import shutil
import subprocess as sp
import sys
from pathlib import Path

from Bio import SeqIO
from loguru import logger


def instantiate_dirs(output_dir: str, force: bool) -> None:
    """Checks the output directory
    :param out_dir: output directory path
    :param force: force flag
    :param logger: logger
    :return: out_dir: final output directory
    """

    # Checks the output directory
    # remove outdir on force
    logger.add(lambda _: sys.exit(1), level="ERROR")
    logger.info(f"Checking the output directory {output_dir}")
    if force is True:
        if Path(output_dir).exists():
            shutil.rmtree(output_dir)
        else:
            logger.info(
                "--force was specified even though the output directory does not already exist. Continuing."
            )
    else:
        if Path(output_dir).exists():
            logger.error(
                "Output directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory."
            )

    # instantiate outdir
    if Path(output_dir).exists() is False:
        Path(output_dir).mkdir(parents=True, exist_ok=True)


def validate_fasta(input_fasta: Path) -> None:
    """
    Validates  FASTA input - checks that the input is a FASTA
    """
    logger.info(f"Checking that the input file {input_fasta} is in FASTA format.")
    # to get extension
    with open(input_fasta, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        if any(fasta):
            logger.info(f"{input_fasta} file checked.")
        else:
            logger.error(
                f"Error: {input_fasta} file is not in the FASTA format. Please check your input file."
            )


def validate_fastq(file: Path) -> bool:
    """Checks the input fastq is really a fastq
        :param file: fastq file
    :return: zipped - Boolean whether the input fastq is gzipped.
    """

    # to get extension
    filename, file_extension = os.path.splitext(file)
    # flag for whether file is zipped
    zipped = True
    if file_extension == ".gz":
        # if gzipped
        with gzip.open(file, "rt") as handle:
            fastq = SeqIO.parse(handle, "fastq")
            if any(fastq):
                logger.info(f"FASTQ {file} checked")
            else:
                logger.error(f"Input file {file} is not in the FASTQ format.")
    else:
        zipped = False
        with open(file, "r") as handle:
            fastq = SeqIO.parse(handle, "fastq")
            if any(fastq):
                logger.info(f"FASTQ {file} checked")
            else:
                logger.error(f"Input file {file} is not in the FASTQ format.")
    return zipped


def check_memory_limit(memory_limit):
    """
    checks if memory limit is a string that ends with K, G or M.
    """

    memory_limit = memory_limit.upper()
    # -1 gets the last char
    if not any(x in memory_limit[-1] for x in ["K", "G", "M"]):
        logger.error(
            f"Memory limit format must be in format of e.g. 2GB or 2MB (GB being gigabytes of RAM and MB being megabytes of RAM)"
        )


def check_dependencies() -> None:
    """Checks the version of freebayes, bwa and samtools
    :return:
    """

    # samtools
    try:
        process = sp.Popen(["samtools", "--version"], stdout=sp.PIPE, stderr=sp.PIPE)
        samtools_out, _ = process.communicate()
        samtools_out = samtools_out.decode()
        samtools_version = samtools_out.split("\n")[0].split(" ")[1]
        message = f"Samtools v{samtools_version} found."
        logger.info(message)
    except Exception:
        logger.error("Samtools not found.")

    # freebayes
    try:
        process = sp.Popen(["freebayes", "--version"], stdout=sp.PIPE, stderr=sp.PIPE)
        freebayes_out, _ = process.communicate()
        freebayes_out = freebayes_out.decode()
        freebayes_version = freebayes_out.split("\n")[0].split(" ")[2]
        message = f"freebayes {freebayes_version} found."
        logger.info(message)
    except Exception:
        logger.error("freebayes not found.")

    # bwa

    try:
        process = sp.Popen(["bwa"], stdout=sp.PIPE, stderr=sp.PIPE)
        bwa_out, _ = process.communicate()
        # a bit weird
        bwa_out = _.decode()
        version_line = []
        for line in bwa_out.split("\n"):
            if "Version" in line:
                version_line.append(line)
        bwa_version = version_line[0].split(" ")[1]
        message = f"bwa v{bwa_version} found."
        logger.info(message)
    except Exception:
        logger.error("bwa not found")
