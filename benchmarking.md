# Benchmarking

I benchmarked `pypolca` v0.2.0 against POLCA (in MaSuRCA v4.1.0) using a variety of isolates used to benchmark [hybracter](https://github.com/gbouras13/hybracter).

Overall, there were 18 isolates: 12 from this [preprint](https://www.biorxiv.org/content/10.1101/2023.09.25.559359v1) by Lerminiaux et al, JKD6159 from Wick et al [here](https://doi.org/10.1128/mra.01129-22) and 5 ATCC strain FASTQs provided by Ryan Wick and Lousie Judd from [here](https://rrwick.github.io/2023/05/05/ont-only-accuracy-with-r10.4.1.html).

All output and this script can be found at the following link on Zenodo [here](https://zenodo.org/doi/10.5281/zenodo.10072192). The FASTQs will be available separately (soon) with the forthcoming hybracter manuscript and benchmarking.

Benchmarking was conducted on an Intel® Core™ i7-10700K CPU @ 3.80 GHz on a machine running Ubuntu 20.04.6 LTS.

The full methodology in how to download the FASTQs and FASTAs and subsample the FASTQs for the following isolates can be found at the hybracter benchmarking repository https://github.com/gbouras13/hybracter_benchmarking. They were all subsampled to 100x genome size. 

The assemblies polished were intermediate assemblies from [Flye](https://github.com/fenderglass/Flye) v2.9.2 done within hybracter.

The dependencies were:

* `pypolca`: Freebayes v1.3.6, bwa v 0.7.17-r1188 and Samtools v1.18
* `POLCA`: Freebayes v1.3.1-dirty, bwa v 0.7.17-r1188 and Samtools v0.1.20

# Results & Conclusion

* There were 16/18 identical assemblies.
* ATCC_33560 had 2 SNPs between `pypolca` and POLCA genomes.
* Lerminiaux_isolateI had 2 SNPs and 1 Indel between `pypolca` and POLCA genomes.
* I also ran those 2 isolates with `pypolca` specifying Freebayes v1.3.1-dirty, with identical results to v1.3.6.
* I also attempted to install Samtools v0.1.20 with `pypolca`, but this was not available on bioconda. While I could attempt a source install etc, I decided I wasn't going to for wider purposes (e.g. including in hybracter) as bioconda integration is essential.
  * I also tried to install Samtools v1.18 with POLCA, but this throws an error and won't work without modifying the polca.sh script - as the parameters of  `samtools sort` have changed (line 136 [here](https://github.com/alekseyzimin/PacBio/blob/1efa908ed3127226aed40448f8d6edd20f798e53/src_reconcile/polca.sh#L136C1-L136C1)).

* Therefore, the likely difference is most likely to be explained by the different Samtools versions.
* Given `pypolca` runs (by default) a newer version of Freebayes and Samtools and that there are very few differences compared to POLCA, I am deciding to move forward with the `pypolca` implementation.
* I have decided to use the newest versions of freebayes and Samtools possible rather than the version installed with MaSuRCA, for ease of maintenance and particularly because the version of Samtools used is a major version behind and the CLI has changed. 

### Setting up the environments

```
mamba create -n pypolcaENV pypolca==0.2.0
mamba create -n polcaENV masurca==4.1.0
```

### Running the commands

* To run pypolca starting inside the `pypolca_vs_POLCA_benchmarking` directory:

```
mkdir -p pypolca_outputs

conda activate pypolcaENV

samples=(
  ATCC_10708
  ATCC_17802
  ATCC_25922
  ATCC_33560
  ATCC_BAA_679
  JKD6159
  Lerminiaux_isolateA
  Lerminiaux_isolateB
  Lerminiaux_isolateC
  Lerminiaux_isolateD
  Lerminiaux_isolateE
  Lerminiaux_isolateF
  Lerminiaux_isolateG
  Lerminiaux_isolateH
  Lerminiaux_isolateI
  Lerminiaux_isolateJ
  Lerminiaux_isolateK
  Lerminiaux_isolateL
)

for sample in "${samples[@]}"; do

pypolca run -a flye_assemblies/${sample}_flye.fasta -1 all_sr_fastqs/${sample}_100x_1.fastq.gz -2 all_sr_fastqs/${sample}_100x_2.fastq.gz -t 8 -o pypolca_outputs/${sample} -p ${sample} -f
done 

conda deactivate

```

* To run `pypolca` on ATCC_33560 and  Lerminiaux_isolateI with Freebayes v1.3.1 and Samtools

```

mamba create -n pypolcafreebayes1.3.1 freebayes==1.3.1 samtools==0.1.20 bwa pypolca==0.2.0
conda activate pypolcafreebayes1.3.1

samples=(
  ATCC_33560
  Lerminiaux_isolateI
)

for sample in "${samples[@]}"; do

pypolca run -a flye_assemblies/${sample}_flye.fasta -1 all_sr_fastqs/${sample}_100x_1.fastq.gz -2 all_sr_fastqs/${sample}_100x_2.fastq.gz -t 8 -o pypolca_outputs/${sample}_freebayes1.3.1 -p ${sample} -f
done 


```

* To run POLCA

```

conda activate polcaENV

mkdir POLCA_outputs

polca.sh -a flye_assemblies/JKD6159_flye.fasta -t 8 -r 'all_sr_fastqs/JKD6159_100x_1.fastq.gz all_sr_fastqs/JKD6159_100x_2.fastq.gz'
mv JKD6159* POLCA_outputs
polca.sh -a flye_assemblies/ATCC_17802_flye.fasta -t 8 -r 'all_sr_fastqs/ATCC_17802_100x_1.fastq.gz all_sr_fastqs/ATCC_17802_100x_2.fastq.gz'
mv ATCC_17802* POLCA_outputs
polca.sh -a flye_assemblies/ATCC_25922_flye.fasta -t 8 -r 'all_sr_fastqs/ATCC_25922_100x_1.fastq.gz all_sr_fastqs/ATCC_25922_100x_2.fastq.gz'
mv ATCC_25922* POLCA_outputs
polca.sh -a flye_assemblies/ATCC_33560_flye.fasta -t 8 -r 'all_sr_fastqs/ATCC_33560_100x_1.fastq.gz all_sr_fastqs/ATCC_33560_100x_2.fastq.gz'
mv ATCC_33560* POLCA_outputs
polca.sh -a flye_assemblies/ATCC_10708_flye.fasta -t 8 -r 'all_sr_fastqs/ATCC_10708_100x_1.fastq.gz all_sr_fastqs/ATCC_10708_100x_2.fastq.gz'
mv ATCC_10708* POLCA_outputs
polca.sh -a flye_assemblies/ATCC_BAA_679_flye.fasta -t 8 -r 'all_sr_fastqs/ATCC_BAA_679_100x_1.fastq.gz all_sr_fastqs/ATCC_BAA_679_100x_2.fastq.gz'
mv ATCC_BAA_679* POLCA_outputs
polca.sh -a flye_assemblies/Lerminiaux_isolateA_flye.fasta -t 8 -r 'all_sr_fastqs/Lerminiaux_isolateA_100x_1.fastq.gz all_sr_fastqs/Lerminiaux_isolateA_100x_2.fastq.gz'
mv Lerminiaux_isolateA* POLCA_outputs
polca.sh -a flye_assemblies/Lerminiaux_isolateB_flye.fasta -t 8 -r 'all_sr_fastqs/Lerminiaux_isolateB_100x_1.fastq.gz all_sr_fastqs/Lerminiaux_isolateB_100x_2.fastq.gz'
mv Lerminiaux_isolateB* POLCA_outputs
polca.sh -a flye_assemblies/Lerminiaux_isolateC_flye.fasta -t 8 -r 'all_sr_fastqs/Lerminiaux_isolateC_100x_1.fastq.gz all_sr_fastqs/Lerminiaux_isolateC_100x_2.fastq.gz'
mv Lerminiaux_isolateC* POLCA_outputs
polca.sh -a flye_assemblies/Lerminiaux_isolateD_flye.fasta -t 8 -r 'all_sr_fastqs/Lerminiaux_isolateD_100x_1.fastq.gz all_sr_fastqs/Lerminiaux_isolateD_100x_2.fastq.gz'
mv Lerminiaux_isolateD* POLCA_outputs
polca.sh -a flye_assemblies/Lerminiaux_isolateE_flye.fasta -t 8 -r 'all_sr_fastqs/Lerminiaux_isolateE_100x_1.fastq.gz all_sr_fastqs/Lerminiaux_isolateE_100x_2.fastq.gz'
mv Lerminiaux_isolateE* POLCA_outputs
polca.sh -a flye_assemblies/Lerminiaux_isolateF_flye.fasta -t 8 -r 'all_sr_fastqs/Lerminiaux_isolateF_100x_1.fastq.gz all_sr_fastqs/Lerminiaux_isolateF_100x_2.fastq.gz'
mv Lerminiaux_isolateF* POLCA_outputs
polca.sh -a flye_assemblies/Lerminiaux_isolateG_flye.fasta -t 8 -r 'all_sr_fastqs/Lerminiaux_isolateG_100x_1.fastq.gz all_sr_fastqs/Lerminiaux_isolateG_100x_2.fastq.gz'
mv Lerminiaux_isolateG* POLCA_outputs
polca.sh -a flye_assemblies/Lerminiaux_isolateH_flye.fasta -t 8 -r 'all_sr_fastqs/Lerminiaux_isolateH_100x_1.fastq.gz all_sr_fastqs/Lerminiaux_isolateH_100x_2.fastq.gz'
mv Lerminiaux_isolateH* POLCA_outputs
polca.sh -a flye_assemblies/Lerminiaux_isolateI_flye.fasta -t 8 -r 'all_sr_fastqs/Lerminiaux_isolateI_100x_1.fastq.gz all_sr_fastqs/Lerminiaux_isolateI_100x_2.fastq.gz'
mv Lerminiaux_isolateI* POLCA_outputs
polca.sh -a flye_assemblies/Lerminiaux_isolateJ_flye.fasta -t 8 -r 'all_sr_fastqs/Lerminiaux_isolateJ_100x_1.fastq.gz all_sr_fastqs/Lerminiaux_isolateJ_100x_2.fastq.gz'
mv Lerminiaux_isolateJ* POLCA_outputs
polca.sh -a flye_assemblies/Lerminiaux_isolateK_flye.fasta -t 8 -r 'all_sr_fastqs/Lerminiaux_isolateK_100x_1.fastq.gz all_sr_fastqs/Lerminiaux_isolateK_100x_2.fastq.gz'
mv Lerminiaux_isolateK* POLCA_outputs
polca.sh -a flye_assemblies/Lerminiaux_isolateL_flye.fasta -t 8 -r 'all_sr_fastqs/Lerminiaux_isolateL_100x_1.fastq.gz all_sr_fastqs/Lerminiaux_isolateL_100x_2.fastq.gz'
mv Lerminiaux_isolateL* POLCA_outputs

```

* To run [dnadiff](https://github.com/marbl/MUMmer3/blob/master/docs/dnadiff.README) to compare the outputs:

```
mamba create -n mummer mummer
conda activate mummer

mkdir -p dnadiff_results

samples=(
  ATCC_10708
  ATCC_17802
  ATCC_25922
  ATCC_33560
  ATCC_BAA_679
  JKD6159
  Lerminiaux_isolateA
  Lerminiaux_isolateB
  Lerminiaux_isolateC
  Lerminiaux_isolateD
  Lerminiaux_isolateE
  Lerminiaux_isolateF
  Lerminiaux_isolateG
  Lerminiaux_isolateH
  Lerminiaux_isolateI
  Lerminiaux_isolateJ
  Lerminiaux_isolateK
  Lerminiaux_isolateL
)


for sample in "${samples[@]}"; do

dnadiff -p dnadiff_results/$sample POLCA_outputs/${sample}_flye.fasta.PolcaCorrected.fa pypolca_outputs/$sample/${sample}_corrected.fasta

done

samples=(
  ATCC_33560
  Lerminiaux_isolateI
)

for sample in "${samples[@]}"; do

dnadiff -p dnadiff_results/${sample}_freebayes1.3.1 POLCA_outputs/${sample}_flye.fasta.PolcaCorrected.fa pypolca_outputs/${sample}_freebayes1.3.1/${sample}_corrected.fasta

done

```