strobealign
==============

Strobealign is a fast short-read aligner. It achieves the speedup by using a dynamic seed size obtained from syncmer-thinned strobemers. Strobealign is multithreaded, implements alignment (SAM) and mapping (PAF), and benchmarked for SE and PE reads of lengths between 100-300bp. A preprint describing **v0.4** is available [here](https://doi.org/10.1101/2021.06.18.449070).

**Current version is 0.5**, implements:
1. Several improvements for downstream SNP andf INDEL calling. See benchmark below.
2. Option to report secondary alignments and more. 
3. Option to set base level alignment parameters. 
4. And more (See release notes)


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

### Zlib linking

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

VARIANT CALLING BENCHMARK
---------------

A small SNV and INDEL calling benchmark is provided below. The experiment evaluates SNP and short INDEL calling for a simulated repetitive genome. The genome consist of 500 concatenated copies of a 50kbp sequence which is then mutated through INDELS (0.5%) and substitutions (5%) between the copies. 

Then, 2 million paired-end reads (lengths 100, 150, 200, 250, 300) from a related genome with 0.5% SNV rate and 0.5% INDEL rate. The challange is to find the right location for each read pair to predict the SNVs in the simulated reads (similar but not identical to the REPEATS example given in the [preprint](https://doi.org/10.1101/2021.06.18.449070)). In the dataset, there is a total of 78,623 SNVs and 78,015 INDELS in the genome where the reads are simulated from. The precision (P), recall (R), and F-score are computed from these numbers. Results in table below 

| Read length  | Tool        | SNVs (P) | SNVs (R) | SNVs (F-score) | Indels (P) | Indels (R) | Indels (F-score) | Alignment time (s) |
| :---         | :---        |      ---: |       ---:  |       ---: |       ---:  |       ---: |       ---: |       ---: |
| 100          | strobealign |  **99.6**  |  93.5       |  **96.4**  |  **55.6**  | **41.1**  | **47.3**   | **455** |
| 100          | minimap2    | 92.6       |  94.3       | 93.4       |  55.2      | 39.1      |  45.8      | 605     |
| 100          | bwa mem     | 93.9       |  **95.9**   |  94.9      |  55.3      | 30.0      | 38.9       | 1020    |

<!-- | 150          | strobealign |   | | |
| 150          | minimap2    | 
| 150          | bwa mem     |  -->

For the results, we ran

```
bcftools mpileup -O z --fasta-ref ref aligned.bam > aligned.vcf.gz
bcftools call -v -c -O v aligned.vcf.gz > aligned.variants.vcf.gz

# Split into SNP and INDELS
grep "#"  aligned.variants.vcf.gz > aligned.variants.SNV.vcf
grep -E "\t[ACGT]\t[ACGT]\t" aligned.variants.vcf.gz >> aligned.variants.SNV.vcf
grep -v -E "\t[ACGT]\t[ACGT]\t" aligned.variants.vcf.gz > aligned.variants.INDEL.vcf

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

