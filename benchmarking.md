# Benchmarking

I benchmarked `pypolca` against POLCA using a variety of isolates used to benchmark [hybracter](https://github.com/gbouras13/hybracter).

All output can be found at the following link on Zenodo [here](_).

Benchmarking was conducted on an Intel® Core™ i7-10700K CPU @ 3.80 GHz on a machine running Ubuntu 20.04.6 LTS.

The full methodology in how to download the FASTQs and FASTAs and subsample the FASTQs for the following isolates can be found at the hybracter benchmarking repository [here](_). They were all subsampled to 100x genome size. 

The assemblies polished were intermediate assemblies from [Flye](https://github.com/fenderglass/Flye) v2.9.2 done within hybracter.

### Setting up the environments

```
mamba create -n pypolcaENV pypolca==0.2.0
mamba create -n polcaENV masurca==4.1.0
```

### running the commands

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







```
mamba create -n mummer mummer
conda activate mummer


dnadiff -p <> <> <>


```