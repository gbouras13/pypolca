#!/usr/bin/env python3

import re
from pathlib import Path
from typing import Dict

from Bio import SeqIO
from Bio.Seq import Seq
from loguru import logger

from pypolca.utils.util import copy_file


def fix_consensus_from_vcf(ref_contigs: Path, vcf: Path, out_fasta: Path) -> None:
    """
    Fix errors in the consensus called in a VCF file by FreeBayes.

    Args:
        ref_contigs (Path): Path to the reference contigs in FASTA format.
        vcf (Path): Path to the VCF file containing variants.
        out_fasta (Path): Path to the output FASTA file.

    Returns:
        None
    """

    rseq = {}
    ctg, seq = "", ""
    rseq = {}

    with open(ref_contigs, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    rseq[ctg] = seq
                ctg = line[1:]
                seq = ""
            else:
                seq += line

    if seq:
        rseq[ctg] = seq

    # Initialize variables
    ctg = ""
    fixes = []
    originals = []
    offsets = []

    total_count = 0

    # Read and process VCF file
    with open(vcf) as vcf_file:
        for line in vcf_file:
            line = line.strip()
            if line.startswith("#"):
                continue

            f = line.split()
            if "," in f[4] or f[0] not in rseq:
                continue

            # if the contig is new it will instantaite a new contig - is reading line by line so needed
            if f[0] != ctg:
                if fixes:
                    if ctg not in rseq:
                        raise Exception(
                            f"sequence {ctg} not found in the input fasta file"
                        )
                    oldseq = rseq[ctg]

                    # Proceed with fixing
                    for i in range(
                        len(fixes) - 1, -1, -1
                    ):  # Going in reverse order to avoid shifting sequence due to indels
                        # First, we check if the sequence at the given offset matches the original variant
                        original_seq = oldseq[
                            offsets[i] - 1 : offsets[i] - 1 + len(originals[i])
                        ]
                        if (
                            any(c in "acgtnACGTN" for c in original_seq)
                            and not original_seq.upper() == originals[i].upper()
                        ):
                            logger.warning(
                                "WARNING! Sequence does not match the original:",
                                ctg,
                                original_seq,
                                originals[i],
                                offsets[i],
                            )
                        else:
                            # Then substitute
                            oldseq = (
                                oldseq[: offsets[i] - 1]
                                + fixes[i]
                                + oldseq[offsets[i] - 1 + len(originals[i]) :]
                            )

                    rseq[ctg] = oldseq

                fixes = []
                originals = []
                offsets = []
                ctg = f[0]

            # append if meets criteria for POLCA
            ff = f[9].split(":")
            if int(ff[5]) > 1:
                if int(ff[5]) >= 2 * int(ff[3]):
                    fixes.append(f[4])
                    originals.append(f[3])
                    offsets.append(int(f[1]))
                    total_count += 1

    # actually fix the report now
    # if fixes as previously means that this step wouldn't continue if the last contig had no changes
    # therefore, use a total count variable to achieve this.
    if total_count > 0:
        logger.info(f"POLCA has found variants. Fixing")
        # Proceed with fixing
        oldseq = rseq[ctg]
        for i in range(
            len(fixes) - 1, -1, -1
        ):  # Going in reverse order to avoid shifting sequence due to indels
            if ctg not in rseq:
                raise Exception(f"sequence {ctg} not found in the input fasta file")
            original_seq = oldseq[offsets[i] - 1 : offsets[i] - 1 + len(originals[i])]
            if (
                any(c in "acgtnACGTN" for c in original_seq)
                and not original_seq.upper() == originals[i].upper()
            ):
                logger.warning(
                    "WARNING! Sequence does not match the original:",
                    ctg,
                    original_seq,
                    originals[i],
                    offsets[i],
                )
            else:
                oldseq = (
                    oldseq[: offsets[i] - 1]
                    + fixes[i]
                    + oldseq[offsets[i] - 1 + len(originals[i]) :]
                )

        rseq[ctg] = oldseq

        records = []
        for contig, sequence in rseq.items():
            record = SeqIO.SeqRecord(seq=Seq(sequence), id=contig, description="")
            records.append(record)

        with open(out_fasta, "w") as output_file:
            SeqIO.write(records, output_file, "fasta")

    else:
        logger.warning(f"POLCA has found 0 variants.")
        logger.warning(
            f"The corrected FASTA {out_fasta} will be the same as the input FASTA {ref_contigs}."
        )
        # copy the reference to the output
        copy_file(ref_contigs, out_fasta)
