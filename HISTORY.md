History
=======

0.5.0 (2026-09-09)
------------------

* Adds `-s`/`--reads_se` to provide extra single end reads (such as the unpaired reads kept by fastp) alongside `-1` and `-2`. Implements [#8](https://github.com/gbouras13/pypolca/issues/8) - thanks @oschwengers for the request and for your patience.
* **This is a quality-of-life change only.** It exists so you don't need to throw away the unpaired reads your QC tool kept, and for parity with POLCA, which accepts them via `polca.sh -a assembly.fna -r 'R1.fastq R2.fastq SE.fastq'`. Please do not expect it to improve your polished assembly - these files are usually well under 1% of the library.
* Implementation: the single end reads are aligned separately with `bwa mem` and their alignments concatenated onto the paired alignments with `samtools cat` before sorting and variant calling.
* One caveat worth knowing. Unpaired reads kept by QC tools tend to be lower quality than the paired reads they came from, because they are exactly the reads that only just survived filtering. On the `C347` sample bundled with the tests, the single end reads were 0.86% of the library but had 22% of bases below Q20, against 7% for the paired reads. At the default thresholds they added 41 extra substitutions (97 to 138) by tipping marginal positions over `--min_alt 2`. With `--careful` they changed almost nothing (73 substitutions either way). If you use `-s`, `--careful` is recommended.

0.4.2 (2026-09-09)
------------------

* `pypolca run` now warns if you have not specified `--careful`, pointing at the benchmarking in [Bouras et al. (2024)](https://doi.org/10.1099/mgen.0.001254). The README has always recommended `--careful`; this just makes the recommendation visible at runtime.
* The warning is skipped if you have set `--min_alt` or `--min_ratio` yourself, since that is a deliberate choice - telling someone running `--min_alt 10` to loosen to 4 would be bad advice.
* Nothing about the default behaviour has changed. `--careful` is still opt-in.

0.4.1 (2026-09-09)
------------------

* Fixes [#28](https://github.com/gbouras13/pypolca/issues/28) - `samtools index` failing with return code 1.
* The `samtools index` command was built as `samtools index -@ <bam>`, with the thread count missing after `-@`. Samtools therefore consumed the BAM path as the argument to `-@`, leaving it with no input file.
    * With Samtools v1.4-v1.15, this exited with code 1 and aborted the run - the error reported in the issue.
    * With Samtools <=v1.3.1, `-@` was not a valid option for `samtools index` at all, giving `index: invalid option -- '@'` and exiting 1.
    * With Samtools >=v1.16, it printed the usage message to stdout and exited 0, so `pypolca` continued silently without ever writing the `.bai` index. Polishing results were unaffected, because `freebayes` is run over the whole BAM and does not need an index.
* Thanks @SuPrSh, @jyap21 and @stevebaeyen for reporting.

0.4.0 (2025-08-19)
------------------

* Adds `--homopolymers` as an option for using closely related (but not identical) short-read data for polishing only homopolymers
* Thanks @rrwick for implementing this!
* Specifically, the `--homopolymers` option tells Pypolca to ignore all changes except those changing the length of homopolymers of at least the given length.
* For example, with `--homopolymers 6`:
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

0.3.1 (2024-02-02)
------------------

* Adds `--careful` to log file to make it clear if it has been run

0.3.0 (2024-01-17)
------------------

* Includes a number of changes by @rrwick
* Fixes errors with counts in the report
* Fixes issue where pypolca would crash if there were spaces in the FASTA header
* Adds `--min_alt`, `--min_ratio` and `--careful` parameters
* Using `--careful` is recommended for low read depth (<25x coverage)

0.2.1 (2023-12-19)
------------------

* Fixes bug in generating the report. If there were multiple alleles at a certain site (on a certain line in the vcf), this would crash pypolca. Exception handling has been added.

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
