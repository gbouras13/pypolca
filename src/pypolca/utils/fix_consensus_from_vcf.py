#!/usr/bin/env python3

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

    # Read reference sequences into memory
    with open(ref_contigs) as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            rseq[record.id] = str(record.seq)

    # Initialize variables
    ctg = ""
    fixes = []
    originals = []
    offsets = []
    fixed_sequences = {}

    # Read and process VCF file
    with open(vcf) as vcf_file:
        for line in vcf_file:
            line = line.strip()
            if line.startswith("#"):
                continue
            fields = line.split()
            if "," in fields[4] or fields[0] not in rseq:
                continue
            if fields[0] != ctg:
                if fixes:
                    if ctg not in rseq:
                        raise Exception(
                            "Sequence {} not found in the input fasta file".format(ctg)
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
                            set(original_seq).intersection("acgtnACGTN")
                            and not originals[i].upper() == original_seq.upper()
                        ):
                            logger.warning(
                                "WARNING! Sequence does not match the original {} {} {} {}".format(
                                    ctg, original_seq, originals[i], offsets[i]
                                )
                            )
                        else:
                            # Then, substitute
                            oldseq = (
                                oldseq[: offsets[i] - 1]
                                + fixes[i]
                                + oldseq[offsets[i] - 1 + len(originals[i]) :]
                            )

                    fixed_sequences[ctg] = oldseq

                fixes = []
                originals = []
                offsets = []
                ctg = fields[0]

            ff = fields[9].split(":")
            if int(ff[5]) > 1 and int(ff[5]) >= 2 * int(ff[3]):
                fixes.append(fields[4])
                originals.append(fields[3])
                offsets.append(int(fields[1]))

    if fixes:
        logger.info(f"POLCA has found variants. Fixing")
        # Proceed with fixing for the last contig
        if ctg not in rseq:
            raise Exception("Sequence {} not found in the input fasta file".format(ctg))
        oldseq = rseq[ctg]

        for i in range(
            len(fixes) - 1, -1, -1
        ):  # Going in reverse order to avoid shifting sequence due to indels
            # First, we check if the sequence at the given offset matches the original variant
            original_seq = oldseq[offsets[i] - 1 : offsets[i] - 1 + len(originals[i])]
            if (
                set(original_seq).intersection("acgtnACGTN")
                and not originals[i].upper() == original_seq.upper()
            ):
                logger.warning(
                    "WARNING! Sequence does not match the original {} {} {} {}".format(
                        ctg, original_seq, originals[i], offsets[i]
                    )
                )
            else:
                # Then, substitute
                oldseq = (
                    oldseq[: offsets[i] - 1]
                    + fixes[i]
                    + oldseq[offsets[i] - 1 + len(originals[i]) :]
                )

        fixed_sequences[ctg] = oldseq

        records = []
        for contig, sequence in fixed_sequences.items():
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
