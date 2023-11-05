=======
History
=======

0.2.0 (2023-11-05)
------------------

* Fixes bug where `pypolca` would always warn it found 0 variants.
* Also would result in the polishing not working. 
* Please upgrade if you have used v0.1.1!
* Adds benchmarking of `pypolca` v0.2.0 vs POLCA more thoroughly. 

0.1.1 (2023-10-09)
------------------

* Patch release to fix the release.yaml uploads to pypi.

0.1.0 (2023-10-09)
------------------

* Initial release

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