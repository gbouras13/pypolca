[![CI](https://github.com/gbouras13/pypolca/actions/workflows/ci.yaml/badge.svg)](https://github.com/gbouras13/pypolca/actions/workflows/ci.yaml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/pypolca/badges/version.svg)](https://anaconda.org/bioconda/pypolca)
[![Bioconda Downloads](https://img.shields.io/conda/dn/bioconda/pypolca)](https://img.shields.io/conda/dn/bioconda/pypolca)
[![PyPI version](https://badge.fury.io/py/pypolca.svg)](https://badge.fury.io/py/pypolca)
[![Downloads](https://static.pepy.tech/badge/pypolca)](https://pepy.tech/project/pypolca)


# pypolca

`pypolca` is a Standalone Python re-imnplementation of the POLCA polisher from the [MaSuRCA genome assembly and analysis toolkit](https://github.com/alekseyzimin/masurca).

## Quick Start

```
# creates conda environment with dnaapler
conda create -n pypolca_env polca

# activates conda environment
conda activate pypolca_env

# runs pypolca
pypolca run -a <genome> -1 <R1 short reads file> -2 <R2 short reads file> -t <threads> -o <output directory> 
```

## Table of Contents
- [pypolca](#pypolca)
  - [Quick Start](#quick-start)
  - [Table of Contents](#table-of-contents)
  - [Description](#description)
  - [Installation](#installation)
    - [Conda](#conda)
    - [Pip](#pip)
    - [Source](#source)
  - [Usage](#usage)

## Description

`pypolca` is a python reimplenetation of the POLCA polisher from the [MaSuRCA genome assembly and analysis toolkit](https://github.com/alekseyzimin/masurca).

It was written for a number of reasons:

* MaSuRCA is only available on Linux, not for MacOS.
* The original `polca.sh` script from MaSuRCA was difficult to use because you could not specify an output directory. Additionally, due to its shell implementation, both FASTQ read files needed to be input together as a string
* To use `polca.sh`, you need to install the entire MaSuRCA assembly toolkit.
* POLCA is recommended for long-read only bacterial only polishing (see [this paper](https://doi.org/10.1371/journal.pcbi.1010905)) and I wanted to include it for MacOS and Linux in my assembly tool [hybracter](https://github.com/gbouras13/hybracter).

Note: I do not guarantee `pypolca` will give identical results to POLCA implemented in MaSuRCA. This is because of the different versions of [freebayes](https://github.com/freebayes/freebayes) used.

Note if you really want to replicate POLCA, the latest versions of MaSuRCA uses freebayes `v1.3.1-dirty`.

## Installation

Installation from conda is recommended as this will install all non-python dependencies.

### Conda

`pypolca` is available on bioconda.

```
conda install -c bioconda pypolca
```

### Pip

You can also install the Python components of `pypolca` with pip.

```
pip install pypolca
```

### Source

Alternatively, the development version of `pypolca` can be installed manually via github.

```
git clone https://github.com/gbouras13/pypolca.git
cd pypolca
pip install -e .
pypolca -h
```

You will then need to install the external dependencies separately, which can be found in `build/environment.yaml`

* [bwa](https://github.com/lh3/bwa) >=0.7.17
* [Samtools](https://github.com/samtools/samtools) >=1.18
* [freebayes](https://github.com/freebayes/freebayes) >=1.3.1



## Usage

```
Usage: pypolca [OPTIONS] COMMAND [ARGS]...

Options:
  -h, --help     Show this message and exit.
  -V, --version  Show the version and exit.

Commands:
  citation  Print the citations for polca
  run       Python implementation of the POLCA polisher from MaSuRCA
```

```
Usage: pypolca run [OPTIONS]

  Python implementation of the POLCA polisher from MaSuRCA

Options:
  -h, --help               Show this message and exit.
  -V, --version            Show the version and exit.
  -a, --assembly PATH      Path to assembly contigs or scaffolds.  [required]
  -1, --reads1 PATH        Path to polishing reads R1 FASTQ. Can be FASTQ or
                           FASTQ gzipped. Required.  [required]
  -2, --reads2 PATH        Path to polishing reads R2 FASTQ. Can be FASTQ or
                           FASTQ gzipped. Optional. Only use -1 if you have
                           single end reads.
  -t, --threads INTEGER    Number of threads.  [default: 1]
  -o, --output PATH        Output directory path  [default: output_polca]
  -f, --force              Force overwrites the output directory
  -n, --no_polish          do not polish, just create vcf file, evaluate the
                           assembly and exit
  -m, --memory_limit TEXT  Memory per thread to use in samtools sort, set to
                           2G or more for large genomes  [default: 2G]
  -p, --prefix TEXT        prefix for output files  [default: polca]
```

The polished output FASTA will be `{prefix}_corrected.fasta` in the specified output directory and the POLCA report will be the textfile `{prefix}.report`


