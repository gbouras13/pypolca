import math
from pathlib import Path

from Bio import SeqIO
from loguru import logger


def create_report(
    vcf: Path, genome: Path, report_file: Path, numsub: int, nind: int
) -> None:
    with open(vcf, "r") as file:
        lines = file.readlines()

    total_size = 0
    for record in SeqIO.parse(genome, "fasta"):
        total_size += len(record.seq)

    numerr = numsub + nind
    erate = numerr / total_size

    qual = round(100 - (erate) * 100, 2)

    qv = round(-10 * math.log10(erate + 0.0000000001), 2)

    with open(report_file, "w") as f:
        f.write("Stats BEFORE polishing:\n")
        logger.info("Stats BEFORE polishing:")
        f.write(f"Substitution Errors Found: {numsub}\n")
        logger.info(f"Substitution Errors Found: {numsub}")
        f.write(f"Insertion/Deletion Errors Found: {nind}\n")
        logger.info(f"Insertion/Deletion Errors Found: {nind}")
        f.write(f"Assembly Size: {total_size}\n")
        logger.info(f"Assembly Size: {total_size}")
        f.write(f"Consensus Quality Before Polishing: {qual}\n")
        logger.info(f"Consensus Quality Before Polishing: {qual}")
        f.write(f"Consensus QV Before Polishing: {qv}\n")
        logger.info(f"Consensus QV Before Polishing: {qv}")
