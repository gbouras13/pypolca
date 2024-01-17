[![CI](https://github.com/gbouras13/pypolca/actions/workflows/ci.yaml/badge.svg)](https://github.com/gbouras13/pypolca/actions/workflows/ci.yaml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![DOI](https://zenodo.org/badge/700722839.svg)](https://zenodo.org/badge/latestdoi/700722839)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/pypolca/badges/version.svg)](https://anaconda.org/bioconda/pypolca)
[![Bioconda Downloads](https://img.shields.io/conda/dn/bioconda/pypolca)](https://img.shields.io/conda/dn/bioconda/pypolca)
[![PyPI version](https://badge.fury.io/py/pypolca.svg)](https://badge.fury.io/py/pypolca)
[![Downloads](https://static.pepy.tech/badge/pypolca)](https://pepy.tech/project/pypolca)


# pypolca

`pypolca` is a Standalone Python re-implementation of the POLCA polisher from the [MaSuRCA genome assembly and analysis toolkit](https://github.com/alekseyzimin/masurca).

## Quick Start

```
# creates conda environment with pypolca 
conda create -n pypolca_env pypolca

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
    - [Note of Caution for Large (e.g. Eukaryotic) Genomes](#note-of-caution-for-large-eg-eukaryotic-genomes)
  - [Installation](#installation)
    - [Conda](#conda)
    - [Pip](#pip)
    - [Source](#source)
  - [Usage](#usage)
- [Benchmarking](#benchmarking)
- [Citation](#citation)

## Description

`pypolca` is a python reimplenetation of the POLCA polisher from the [MaSuRCA genome assembly and analysis toolkit](https://github.com/alekseyzimin/masurca) that was made for inclusion into the hybrid bacterial genome assembly tool [hybracter](https://github.com/gbouras13/hybracter).

It was written for a number of reasons:

* MaSuRCA is only available on Linux, not for MacOS.
* The original `polca.sh` script from MaSuRCA was difficult to use because you could not specify an output directory. Additionally, due to its shell implementation, both FASTQ read files needed to be input together as a string
* To use `polca.sh`, you need to install the entire MaSuRCA assembly toolkit.
* POLCA is recommended for long-read first bacterial assembly polishing (see [this paper](https://doi.org/10.1371/journal.pcbi.1010905)) and I wanted to include it for MacOS and Linux in my assembly tool [hybracter](https://github.com/gbouras13/hybracter).

Note: I neither guarantee nor desire that `pypolca` will give identical results to POLCA implemented in MaSuRCA. This is because of the different versions of [freebayes](https://github.com/freebayes/freebayes) and Samtools that might be used as a dependency. 

In testing, `pypolca` v0.2.0 (running Freebayes v1.3.6 and Samtools v1.18) was extremely similar, but not identical to POLCA (running Freebayes v1.3.1-dirty and Samtools v0.1.20). Please see [benchmarking](benchmarking.md) for more details.

I have decided to use the newest versions of freebayes and Samtools possible rather than the version installed with MaSuRCA, for ease of maintenance and particularly because the version of Samtools used is a major version behind and the CLI has changed. 

### Note of Caution for Large (e.g. Eukaryotic) Genomes

* I have implemeted `pypolca` predominantly for the use-case of polishing long-read bacterial genome assemblies with short reads. Therefore, I decided not to implement the batched multiprocessing of freebayes included in POLCA, because it was a lot of work for no benefit for most bacterial genomes. 
* However, this is certainly not true for larger genomes such as eukaryotic organisms. `pypolca` should be a lot slower than POLCA for such organisms if you run both with more than 1 thread. 
* I do not intend to implement multiprocessing but if someone wants to feel free to make a PR.

## Installation

Installation from conda is recommended as this will install all non-python dependencies.

### Conda

`pypolca` will soon be available on bioconda.

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

If you have install `pypolca` with pip or from source, you will then need to install the external dependencies separately, which can be found in `build/environment.yaml`

* [bwa](https://github.com/lh3/bwa) >=0.7.17
* [Samtools](https://github.com/samtools/samtools) >=1.18
* [freebayes](https://github.com/freebayes/freebayes) >=1.3.1,<1.3.7

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
  --min_alt INTEGER        Minimum alt allele count to make a change
                           [default: 2]
  --min_ratio FLOAT        Minimum alt allele to ref allele ratio to make a
                           change  [default: 2.0]
  --careful                Equivalent to --min_alt 4 --min_ratio 3
  -n, --no_polish          do not polish, just create vcf file, evaluate the
                           assembly and exit
  -m, --memory_limit TEXT  Memory per thread to use in samtools sort, set to
                           2G or more for large genomes  [default: 2G]
  -p, --prefix TEXT        prefix for output files  [default: polca]
```

The polished output FASTA will be `{prefix}_corrected.fasta` in the specified output directory and the POLCA report will be the textfile `{prefix}.report`

# Benchmarking

Please see [benchmarking](benchmarking.md) for more details. As can be seen, `pypolca` v0.2.0 was extremely similar, but not identical to POLCA.


# Citation

Please cite `pypolca` in your paper using:

Bouras G, Wick RR (2023) pypolca: Standalone Python reimplementation of the genome polishing tool POLCA. https://github.com/gbouras13/pypolca. 

Zimin AV, Salzberg SL (2020) The genome polishing tool POLCA makes fast and accurate corrections in genome assemblies. PLoS Comput Biol 16(6): e1007981. https://doi.org/10.1371/journal.pcbi.1007981.

