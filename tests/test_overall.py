"""
Unit tests for dnaapler overall

Usage: pytest .

"""

# import
import os
import shutil

# import functions
import subprocess
import unittest
from pathlib import Path

import pytest

test_data = Path("tests/test_data")

def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")

# define threads 
threads: int = 4

def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None



def test_C347(tmp_dir):
    """test C347"""
    input_fasta: Path = f"{test_data}/chromosome.fasta"
    r1: Path = f"{test_data}/C347_R1.fastq.gz"
    r2: Path = f"{test_data}/C347_R2.fastq.gz"
    outdir: Path = "output_dir"
    cmd = f"pypolca run -a {input_fasta} -1 {r1} -2 {r2} -t {threads} -o {outdir} -f"
    exec_command(cmd)
    remove_directory(outdir)

def test_C347_single(tmp_dir):
    """test C347 single"""
    input_fasta: Path = f"{test_data}/chromosome.fasta"
    r1: Path = f"{test_data}/C347_R1.fastq.gz"
    outdir: Path = "output_dir"
    cmd = f"pypolca run -a {input_fasta} -1 {r1} -t {threads} -o {outdir} -f"
    exec_command(cmd)
    remove_directory(outdir)

def test_citation():
    """test citation"""
    cmd = f"dnaapler citation"
    exec_command(cmd)


class TestExits(unittest.TestCase):
    """Tests of End to End common failures"""

    def test_FASTQ_as_FASTA(self):
        """test FASTQ as FASTA to fail"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{test_data}/chromosome.fasta"
            r1: Path = f"{test_data}/C347_R1.fastq.gz"
            outdir: Path = "output_dir"
            cmd = f"pypolca run -a {r1} -1 {r1} -t {threads} -o {outdir} -f"
            exec_command(cmd)
            remove_directory(outdir)

   
