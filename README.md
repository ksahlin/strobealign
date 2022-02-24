strobealign
==============

Strobealign is a fast short-read aligner. It achieves the speedup by using a dynamic seed size obtained from syncmer-thinned strobemers. Strobealign is multithreaded, implements alignment (SAM) and mapping (PAF), and benchmarked for SE and PE reads of lengths between 100-300bp. A preprint describing **v0.4** is available [here](https://doi.org/10.1101/2021.06.18.449070).

**Current version is 0.6.1**. See the performance of v0.6.1 [here](https://github.com/ksahlin/StrobeAlign/edit/main/README.md#v061-performance).

v0.6.1 implements:
1. Runtime bugfix introduced in v0.6.

v0.6 implements:
1. Crucial bugfix to v0.5: Rare but occasional alignments to very long reference regions.
2. Identifying symmetrical hash collisions and testing reverse orientation. This leads to a slightly increased alignment accuracy over previous versions, particularly for shorter read lengths.
3. Fixes reporting of template len field in SAM output if deletion in alignment.

v0.5 implements:
1. Several improvements for downstream SNP and INDEL calling. SNV and small indel calling benchmark below.
2. Option to report secondary alignments. 
3. Base level SW alignment parameters are now parameters to strobealign. 
4. And more.. (See release notes)


INSTALLATION
----------------

You can acquire precompiled binaries for Linux and Mac OSx from the [release page](https://github.com/ksahlin/StrobeAlign/releases) compiled with `-O3 -mavx2`. 

It has been [reported](https://github.com/ksahlin/StrobeAlign/issues/6) that `strobealign` is even faster if compliled with flag `-march=skylake-avx512` for avx512 supported processors.

If you want to compile from the source, you need to have a newer `g++` and [zlib](https://zlib.net/) installed. Then do the following:

```
git clone https://github.com/ksahlin/StrobeAlign
cd StrobeAlign
# Needs a newer g++ version. Tested with version 8 and upwards.
g++ -std=c++14 main.cpp source/index.cpp source/xxhash.c source/ksw2_extz2_sse.c source/ssw_cpp.cpp source/ssw.c -lz -fopenmp -o strobealign -O3 -mavx2
```

### Zlib linking error

If you have `zlib` installed, and the `zlib.h` file is in folder `/path/to/zlib/include` and the `libz.so` file in `/path/to/zlib/lib` but you get 

```
main.cpp:12:10: fatal error: zlib.h: No such file or directory
 #include <zlib.h>
          ^~~~~~~~
compilation terminated.
```

add `-I/path/to/zlib/include -L/path/to/zlib/lib` to the compilation, that is

```
g++ -std=c++14 -I/path/to/zlib/include -L/path/to/zlib/lib main.cpp source/index.cpp source/xxhash.c source/ksw2_extz2_sse.c source/ssw_cpp.cpp source/ssw.c -lz -fopenmp -o strobealign -O3 -mavx2
``` 


USAGE
-------

### Alignment

Strobealign comes with a parameter `-r read_length` that sets suitable seed parameters for the rough read length. Specifically, it sets parameters `-k`, `-l` and `-u`. If not specified, it defaults to 150. The value of `r` does not have to match the exact read length.

For alignment to SAM file:

```
strobealign -r <read_length> ref.fa reads.fa > output.sam
```

To report secondary alignments, set parameter `-N [INT]` for maximum of `[INT]` secondary alignments. 

### Mapping

For mapping to PAF file (option -x):

```
strobealign -r <read_length> -x ref.fa reads.fa > output.sam
```


V0.6.1 PERFORMANCE
---------------

## Mapping accuracy and runtime

Below shows the accuracy (panel A) runtime (panel B) and %-aligned reads (panel C) for the SIM3 and REPEATS dataset in the preprint using strobealign v0.6.1. On all but the 2x100 dataset, it has comparable or higher accuracy than BWA MEM while being substantially faster. On the 2x100 dataset it has the second highest accuracy after BWA MEM while being substantially faster.

![v0 6 1_sim_and_repeats_experiment 001](https://user-images.githubusercontent.com/1714667/155538019-a5be1a72-9d14-4a4d-831c-e5de32d79c34.jpeg)
Fig 1. Accuracy (panel A) runtime (panel B) and %-aligned reads (panel C) for the SIM3 dataset

FIG2 TO BE GENERATED

Fig 2. Accuracy (panel A) runtime (panel B) and %-aligned reads (panel C) for the REPEATS dataset

## Variant calling benchmark
A small SNV and INDEL calling benchmark with strobealign v0.6 is provided below. We used `bcftools` to call SNPs and indels on a simulated repetitive genome based on alignments from strobealign, BWA-MEM, and minimap2. The genome is a 16.8Mbp sequence consisting of 500 concatenated copies of a 40kbp sequence which is mutated through substitutions (5%) and removing segments of size 1bp-1kbp (0.5%) along the oringinal 20Mbp string. 

Then, 2 million paired-end reads (lengths 100, 150, 200, 250, 300) from a related genome with high variation rate: 0.5% SNVs and 0.5% INDELs. The challange is to find the right location of reads in the repetitive genome to predict the SNVs and INDELs in the related genome. In the genome where the reads are simulated from there is about 78k SNVs and INDELS, respectively. The precision (P), recall (R), and F-score are computed from these numbers. Results in table below. 

In the experiments strobealign is in general the fastest tool, has the highest SNV precision, and *highest precision, recall, and F-score* for indels. 

There are frequent indels in this dataset (every 200th bases on average) requiring calls to base level alignments for most reads. Between 65-85% of strobealign's runtime is spent on base level alignments with third-party SSW alignment module. The longer the reads the higher % of time is spent on base level alignment. Speed improvements to base-level alignment libraries will greatly reduce runtime on this dataset.


| Read length  | Tool        | SNVs (P) | SNVs (R) | SNVs (F-score) | Indels (P) | Indels (R) | Indels (F-score) | Alignment time (s) |
| :---         | :---        |      ---: |       ---:  |       ---: |       ---:  |       ---: |       ---: |       ---: |
| 100 | strobealign | **97.9** | 93.5 | **95.6** | **55.6** | **41.1** | **47.2** | **424** |
| &nbsp; | minimap2 | 91.4 | 94.3 | 92.8 | 55.2 | 39.1 | 45.8 | 605 |
| &nbsp; | bwa_mem | 93.7 | **95.9** | 94.8 | 55.3 | 30.0 | 38.9 | 1020 |
| &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
| 150 | strobealign | **96.6** | 92.7 | 94.6 | **55.2** | **46.2** | **50.3** |  **350** |
| &nbsp; | minimap2 | 89.8 | 94.6 | 92.1 | 54.9 | 44.8 | 49.3 | 902 |
| &nbsp; | bwa_mem | 96.0 | **96.0** | **96.0** | 55.0 | 39.6 | 46.1 | 1010|
| &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
| 200 | strobealign | **97.4** | 94.1 | 95.7 | **55.3** | **45.8** | **50.1** | **487**|
| &nbsp; | minimap2 | 88.1 | **96.7** | 92.2 | 55.0 | 44.7 | 49.3 | 1290 |
| &nbsp; | bwa_mem | 95.2 | 96.5 | **95.8** | 55.1 | 42.3 | 47.8 | 1263 |
| &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
| 250 | strobealign | **96.4** | 93.3 | 94.8 | **55.1** | **45.0** | **49.6** | **697** |
| &nbsp; | minimap2 | 87.7 | 94.8 | 91.1 | 54.9 | 43.8 | 48.7 | 998 |
| &nbsp; | bwa_mem | 94.3 | **96.2** | **95.2** | **55.1** | 42.3 | 47.8 | 1593 |
| &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
| 300 | strobealign | **95.7** | 92.7 | 94.1 | **55.1** | **44.5** | **49.2** | **1005** |
| &nbsp; | minimap2 | 88.2 | 94.3 | 91.2 | 54.8 | 43.4 | 48.4 | 1046 |
| &nbsp; | bwa_mem | 93.7 | **96.4** | **95.0** | 54.9 | 42.0 | 47.6 | 1988 |


For the results, we ran

```
bcftools mpileup -O z --fasta-ref ref aligned.bam > aligned.vcf.gz
bcftools call -v -c -O v aligned.vcf.gz > aligned.variants.vcf.gz

# Split into SNP and INDELS
grep -v -E -e "INDEL;" aligned.variants.vcf.gz > aligned.variants.SNV.vcf
grep "#"  aligned.variants.vcf.gz > aligned.variants.INDEL.vcf
grep -E -e "INDEL;" aligned.variants.vcf.gz >> aligned.variants.INDEL.vcf

for type in SNV INDEL
do
	bcftools sort -Oz aligned.variants.$type.vcf.gz -o aligned.variants.sorted.$type.vcf.gz
	bcftools index aligned.variants.sorted.$type.vcf.gz
	bcftools isec --nfiles 2 -O u true_variants.sorted.$type.vcf.gz  aligned.variants.sorted.$type.vcf -p out_$type
done
```


CREDITS
----------------

Kristoffer Sahlin. Flexible seed size enables ultra-fast and accurate read alignment. bioRxiv, 2021. doi:10.1101/2021.06.18.449070. Preprint available [here](https://doi.org/10.1101/2021.06.18.449070).


VERSION INFO
---------------

See [release page](https://github.com/ksahlin/StrobeAlign/releases)


LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/uLTRA/blob/master/LICENCE.txt).

