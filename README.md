[![CI](https://github.com/gbouras13/pypolca/actions/workflows/ci.yaml/badge.svg)](https://github.com/gbouras13/pypolca/actions/workflows/ci.yaml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![DOI](https://zenodo.org/badge/700722839.svg)](https://zenodo.org/badge/latestdoi/700722839)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/pypolca/badges/version.svg)](https://anaconda.org/bioconda/pypolca)
[![Bioconda Downloads](https://img.shields.io/conda/dn/bioconda/pypolca)](https://img.shields.io/conda/dn/bioconda/pypolca)
[![PyPI version](https://badge.fury.io/py/pypolca.svg)](https://badge.fury.io/py/pypolca)
[![Downloads](https://static.pepy.tech/badge/pypolca)](https://pepy.tech/project/pypolca)


# pypolca

`pypolca` is a Standalone Python re-implementation of the POLCA polisher from the [MaSuRCA genome assembly and analysis toolkit](https://github.com/alekseyzimin/masurca).

`pypolca` also adds a `--careful` flag that we have shown is (generally) the best performing short-read bacterial isolate genome polishing option across a range of depths and samples.

The `‑‑careful` flag requires at least 4 reads to support the alternative allele, and a minimum ratio of 3 (i.e. at least three times as many reads must support the alterative allele compared to the reference allele). This prevents almost all false positives at low depths without sacrificing error removal sensitivity.

**Therefore, we recommend you always use `pypolca` with `--careful`**

## Manuscript

For more information about `pypolca` and `pypolca --careful`, along with [Polypolish](https://github.com/rrwick/Polypolish)'s new `--careful` option and our analysis of short-read polishing methods for near-perfect Nanopore assemblies, please read our manuscript introducing `pypolca`:

Bouras G, Judd LM, Edwards RA, Vreugde S, Stinear TP, Wick RR. How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. Microbial Genomics. 2024. doi: [https://doi.org/10.1099/mgen.0.001254](https://doi.org/10.1099/mgen.0.001254).

## Quick Start

```
# creates conda environment with pypolca 
conda create -n pypolca_env pypolca

# activates conda environment
conda activate pypolca_env

# runs pypolca with --careful
pypolca run -a <genome> -1 <R1 short reads file> -2 <R2 short reads file> -t <threads> -o <output directory> --careful
```

## Table of Contents
- [pypolca](#pypolca)
  - [Manuscript](#manuscript)
  - [Quick Start](#quick-start)
  - [Table of Contents](#table-of-contents)
  - [Description](#description)
    - [Note of Caution for Large (e.g. Eukaryotic) Genomes](#note-of-caution-for-large-eg-eukaryotic-genomes)
  - [Installation](#installation)
    - [Conda](#conda)
    - [Pip](#pip)
    - [Container](#container)
    - [Source](#source)
  - [Usage](#usage)
  - [Homopolymer-only mode](#homopolymer-only-mode)
  - [Citation](#citation)

## Description

`pypolca` is a python reimplenetation of the POLCA polisher from the [MaSuRCA genome assembly and analysis toolkit](https://github.com/alekseyzimin/masurca) that was made for inclusion into the hybrid bacterial genome assembly tool [hybracter](https://github.com/gbouras13/hybracter).

It was written for a number of reasons:

* MaSuRCA is only available on Linux, not for MacOS.
* The original `polca.sh` script from MaSuRCA was difficult to use because you could not specify an output directory. Additionally, due to its shell implementation, both FASTQ read files needed to be input together as a string
* To use `polca.sh`, you need to install the entire MaSuRCA assembly toolkit.
* POLCA is recommended for long-read first bacterial assembly polishing (see [this paper](https://doi.org/10.1371/journal.pcbi.1010905)) and I wanted to include it for MacOS and Linux in [hybracter](https://github.com/gbouras13/hybracter).

Note: I neither guarantee nor desire that `pypolca` will give identical results to POLCA implemented in MaSuRCA. This is because of the different versions of [freebayes](https://github.com/freebayes/freebayes) and Samtools that might be used as a dependency. 

In testing, `pypolca` v0.2.0 (running Freebayes v1.3.6 and Samtools v1.18) was extremely similar, but not identical to POLCA (running Freebayes v1.3.1-dirty and Samtools v0.1.20). Please see [benchmarking](benchmarking.md) for more details.

I have decided to use the newest versions of freebayes and Samtools possible rather than the version installed with MaSuRCA, for ease of maintenance and particularly because the version of Samtools used is a major version behind and the CLI has changed. 

Note: as of 21 March 2025, I bumped the dependency of freebayes to v1.3.9. Testing shows identical results to v1.3.6 and it is compatible with Apple Silicon MacOSX builds.(As an aside, v1.3.7 conda package appears to be broken)

### Note of Caution for Large (e.g. Eukaryotic) Genomes

* I have implemeted `pypolca` predominantly for the use-case of polishing long-read bacterial genome assemblies with short reads. Therefore, I decided not to implement the batched multiprocessing of freebayes included in POLCA, because it was a lot of work for no benefit for most bacterial genomes. 
* However, this is certainly not true for larger genomes such as eukaryotic organisms. `pypolca` should be a lot slower than POLCA for such organisms if you run both with more than 1 thread. 
* I do not intend to implement multiprocessing but if someone wants to feel free to make a PR.

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

### Container

If you have Docker/Singularity/Apptainer installed, you can use the [biocontainers container](https://quay.io/repository/biocontainers/pypolca?tab=tags) (yes, every bioconda package has one!)

For example to install `pypolca v0.3.1` with Singularity

```
IMAGE_DIR="<the directory you want the .sif file to be in >"
singularity pull --dir $IMAGE_DIR docker://quay.io/biocontainers/pypolca:0.3.1--pyhdfd78af_0

containerImage="$IMAGE_DIR/pypolca_0.3.1--pyhdfd78af_0.sif"
singularity exec $containerImage pypolca run -h
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
* [freebayes](https://github.com/freebayes/freebayes) >=1.3.9

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
  -o, --output PATH        Output directory path  [default: output_pypolca]
  -f, --force              Force overwrites the output directory
  --min_alt INTEGER        Minimum alt allele count to make a change
                           [default: 2]
  --min_ratio FLOAT        Minimum alt allele to ref allele ratio to make a
                           change  [default: 2.0]
  --careful                Equivalent to --min_alt 4 --min_ratio 3
  --homopolymers INTEGER   Ignore all changes except for homopolymer-length
                           changes, with homopolymers defined by this length
  -n, --no_polish          do not polish, just create vcf file, evaluate the
                           assembly and exit
  -m, --memory_limit TEXT  Memory per thread to use in samtools sort, set to
                           2G or more for large genomes  [default: 2G]
  -p, --prefix TEXT        prefix for output files  [default: pypolca]
```

The polished output FASTA will be `{prefix}_corrected.fasta` in the specified output directory and the POLCA report will be the textfile `{prefix}.report`

## Homopolymer-only mode

The `--homopolymers` option tells Pypolca to ignore all changes except those changing the length of homopolymers of at least the given length.

For example, with `--homopolymers 6`:
* `ACGTA` → `ACATA` ❌ Not applied (change is not in a homopolymer).
* `ACGTTTTTTTCAA` → `ACGTTTTTTTTCAA` ✅ Applied (homopolymer ≥6 bp).
* `GCTAAAATCG` → `GCTAAAAATCG` ❌ Not applied (homopolymer <6 bp).

This mode is useful when polishing an ONT-only assembly without having short reads from the same sample. You can instead use reads from a closely related sample. The underlying assumptions are:
1. ONT-only assemblies often contain errors in homopolymer length, especially for long homopolymers.
2. Short reads will be more accurate than ONT reads for long homopolymers.
3. Homopolymer runs are generally conserved between closely related genomes.

Example command, where 'sample A' is your ONT-only assembly and 'sample B' is a closely related genome with short reads:
```bash
pypolca run -a sample_A_draft.fasta -1 sample_B_1.fastq.gz -2 sample_B_2.fastq.gz -t 16 --careful --homopolymers 6
```

If you don’t have reads from a related sample but do have a reference genome, you can simulate reads with [wgsim](https://github.com/lh3/wgsim):
```bash
pair_count=$(seqtk size sample_B.fasta | awk '{s+=$2} END{print int(s/4)}')  # ~50× depth
wgsim -e 0.0 -r 0.0 -N "$pair_count" -1 100 -2 100 sample_B.fasta temp_1.fastq temp_2.fastq
pypolca run -a sample_A_draft.fasta -1 temp_1.fastq -2 temp_2.fastq -t 16 --careful --homopolymers 6
rm temp_1.fastq temp_2.fastq
```

We recommend a threshold of 6, since modern ONT-only assemblies are typically accurate for shorter homopolymers.

## Citation

Please cite `pypolca` in your paper using both:

Bouras G, Judd LM, Edwards RA, Vreugde S, Stinear TP, Wick RR. How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. Microbial Genomics. 2024. doi: [https://doi.org/10.1099/mgen.0.001254](https://doi.org/10.1099/mgen.0.001254).

Zimin AV, Salzberg SL (2020) The genome polishing tool POLCA makes fast and accurate corrections in genome assemblies. PLoS Comput Biol 16(6): e1007981. [https://doi.org/10.1371/journal.pcbi.1007981](https://doi.org/10.1371/journal.pcbi.1007981).

