#!/usr/bin/env python3

import re
from pathlib import Path
from typing import Dict, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from loguru import logger

from pypolca.utils.util import copy_file


def fix_consensus_from_vcf(
    ref_contigs: Path,
    vcf: Path,
    out_fasta: Path,
    min_alt: int,
    min_ratio: float,
) -> Tuple[int, int]:
    """
    Fix errors in the consensus called in a VCF file by FreeBayes.

    Args:
        ref_contigs (Path): Path to the reference contigs in FASTA format.
        vcf (Path): Path to the VCF file containing variants.
        out_fasta (Path): Path to the output FASTA file.

    Returns:
        total_subs (int): the total number of substitution changes made
        total_indels (int): the total number of indel changes made
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
                ctg = line[1:].split()[0]
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
    total_subs, total_indels = 0, 0

    # Read and process VCF file
    with open(vcf) as vcf_file:
        for line in vcf_file:
            line = line.strip()
            if line.startswith("#"):
                continue

            f = line.split()
            # this line skips the multiple allele cases
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
            ref_count, alt_count = int(ff[3]), int(ff[5])
            ref_seq, alt_seq = f[3], f[4]
            if alt_count >= min_alt:
                if alt_count >= min_ratio * ref_count:
                    fixes.append(alt_seq)
                    originals.append(ref_seq)
                    offsets.append(int(f[1]))
                    subs, indels = edit_distance(ref_seq, alt_seq)
                    total_subs += subs
                    total_indels += indels
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

    return total_subs, total_indels


def edit_distance(s1, s2):
    """
    Calculate the global edit distance between two strings, providing separate counts for
    substitutions and indels.
    """
    s1, s2 = s1.upper(), s2.upper()
    dp = [[0 for n in range(len(s2) + 1)] for m in range(len(s1) + 1)]
    for i in range(len(s1) + 1):
        dp[i][0] = i
    for j in range(len(s2) + 1):
        dp[0][j] = j
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            dp[i][j] = min(dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + cost)
    subs, indels = 0, 0
    i, j = len(s1), len(s2)
    while i > 0 or j > 0:
        if (
            i > 0
            and j > 0
            and s1[i - 1] != s2[j - 1]
            and dp[i][j] == dp[i - 1][j - 1] + 1
        ):
            subs += 1
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + 1:
            indels += 1
            i -= 1
        elif j > 0 and dp[i][j] == dp[i][j - 1] + 1:
            indels += 1
            j -= 1
        else:
            i -= 1
            j -= 1
    return subs, indels
