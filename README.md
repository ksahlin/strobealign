strobealign
===========

Strobealign is a fast **short-read aligner**. It achieves the speedup by using a dynamic seed size obtained from [syncmer](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7869670/)-thinned [strobemers](https://github.com/ksahlin/strobemers#what-is-a-strobemer).

Strobealign is multithreaded, aligns single-end and paired-end reads, and outputs mapped reads either in SAM format (default) or PAF format (extension-free mapping). A paper describing version 0.7.1 is available [here](https://doi.org/10.1186/s13059-022-02831-7). Strobealign is benchmarked for read lengths between 100 and 500 bp. 

See [INSTALLATION](https://github.com/ksahlin/strobealign#installation) and [USAGE](https://github.com/ksahlin/strobealign#usage) to install and run strobealign. See [v07 PERFORMANCE](https://github.com/ksahlin/strobealign#v07-performance) for the accuracy and runtime performance of strobealign.

INSTALLATION
----------------

### Conda
Strobealign can be installed through conda. Simply run

```
conda create -n strobealign strobealign
```

### Binaries
You can acquire precompiled binaries for Linux and macOS from the [release page](https://github.com/ksahlin/strobealign/releases) compiled with `-O3 -mavx2`.

It has been [reported](https://github.com/ksahlin/strobealign/issues/6) that strobealign is even faster if compiled with flag `-march=skylake-avx512` for processors that support AVX512 instructions.

### From source
If you want to compile from the source, you need to have CMake, a recent `g++` (tested with version 8) and [zlib](https://zlib.net/) installed.
Then do the following:

```
git clone https://github.com/ksahlin/strobealign
cd strobealign
cmake -B build
make -j -C build
```

### Development installation

You can add `-DCMAKE_BUILD_TYPE=RelWithDebInfo` to the
`cmake` options to get debug symbols.

Run `make` with `VERBOSE=1` to get more logging output.


USAGE
-----

```
strobealign ref.fa reads.fq > output.sam                # Single-end reads
strobealign ref.fa reads1.fq reads2.fq > output.sam     # Paired-end reads

strobealign -x ref.fa reads.fq > output.paf             # Single-end reads mapping only (PAF)
strobealign -x ref.fa reads1.fq reads2.fq > output.paf  # Paired-end reads mapping only (PAF)
```

To report secondary alignments, set parameter `-N [INT]` for maximum of `[INT]` secondary alignments. 


Index files
-----------

Every time you run strobealign, the program by default computes an index of the
provided reference FASTA. This is done because the indexing parameters vary
according to read length. This is relatively fast: Depending on CPU, indexing a
human-sized genome takes 2-5 minutes, but when you work with multiple libraries
that all (roughly) have the same read length, you can save time by letting
strobealign precompute an index file with the `-i` option.

Since strobealign needs to know the read length, either provide it with the
reads file as if you wanted to map them:

    strobalign -i ref.fa reads.1.fastq.gz reads.2.fastq.gz

Or set the read length explicitly with `-r`:

    strobealign -i ref.fa -r 150

This creates a file named `ref.fa.sti` containing the strobemer index. To use
it, provide option `--use-index` when doing the actual mapping:

    strobealign --use-index ref.fa reads.1.fastq.gz reads.2.fastq.gz | samtools ...


Changelog
---------

See [Changelog](CHANGES.md).


Video
-----

For an introduction, see the [RECOMB-Seq video from 2022: “Flexible seed size enables ultra-fast and accurate read alignment”](https://www.youtube.com/watch?v=cn32telW63w) (12 minutes). For a more detailed presentation of the underlying seeding mechanism in strobealign (strobemers) see [“Efficient sequence similarity searches with strobemers”](https://www.youtube.com/watch?v=DS4tURz1Wio) (73 minutes).


V0.7 PERFORMANCE
----------------

We have in below three sections investigated accuracy and runtime metrics for v0.7 on SIM3 and REPEATS datasets included in the preprint, as well as performance of SNV and small indel calling for additional simulated and biological (GIAB) datasets.

For the biological SNV and indel experiments, we used GIAB datasets (HG004; Mother) with 2x150bp reads (subsampled to ~26x coverage) and 2x250bp reads (~17x coverage).

## Mapping accuracy and runtime

Below shows the accuracy (panel A) runtime (panel B) and %-aligned reads (panel C) for the SIM3 (Fig 1) and REPEATS (Fig 2) datasets in the preprint using strobealign v0.7. On all but the 2x100 datasets, strobealign has comparable or higher accuracy than BWA-MEM while being substantially faster. On the 2x100 datasets, strobealign has the second highest accuracy after BWA-MEM on SIM3 while being substantially faster, and comparable accuracy to minimap2 and BWA-MEM on the REPEATS dataset while being twice as fast.

![v0 6 1_sim3 001 jpeg 001](https://user-images.githubusercontent.com/1714667/155574282-35bd370c-e7f5-4e59-896d-465473e6a71f.jpeg)
Figure 1. Accuracy (panel A) runtime (panel B) and %-aligned reads (panel C) for the SIM3 dataset

![v0 6 1_repeats_experiment 001](https://user-images.githubusercontent.com/1714667/155572146-9d51b822-e9bc-4dda-8703-71ea0306330b.jpeg)
Figure 2. Accuracy (panel A) runtime (panel B) and %-aligned reads (panel C) for the REPEATS dataset

## Variant calling benchmark (simulated REPEATS)
A small SNV and INDEL calling benchmark with strobealign v0.7 is provided below. We used `bcftools` to call SNPs and indels on a simulated repetitive genome based on alignments from strobealign, BWA-MEM, and minimap2 (using one core). The genome is a 16.8Mbp sequence consisting of 500 concatenated copies of a 40kbp sequence which is mutated through substitutions (5%) and removing segments of size 1bp-1kbp (0.5%) along the oringinal 20Mbp string.

Then, 2 million paired-end reads (lengths 100, 150, 200, 250, 300) from a related genome with high variation rate: 0.5% SNVs and 0.5% INDELs. The challange is to find the right location of reads in the repetitive genome to predict the SNVs and INDELs in the related genome. In the genome where the reads are simulated from there is about 78k SNVs and INDELS, respectively. Locations of true SNVs and INDELs and provided by the read simulator. The precision (P), recall (R), and [F-score](https://en.wikipedia.org/wiki/F-score) are computed based on the true variants (for details see section [Variant calling benchmark method](https://github.com/ksahlin/strobealign#variant-calling-benchmark-method)). Results in table below.

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


## Variant calling benchmark (simulated SIM3)

We simulated 2x150 and 2x250 reads at 30x coverage from a human genome with SNV and indel rate according to the SIM3 genome (described in the preprint). We aligned the reads to hg38 without alternative haplotypes as proposed [here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use). We used 16 cores for all aligners.

Results are shown for SNVs and indels separately in Figure 3. For SNVs, predictions with strobealign as the aligner have an [F-score](https://en.wikipedia.org/wiki/F-score) on par with most other aligners. BWA has the best performance on this dataset. However, indel predictions have both the highest recall and precision using strobealign. Minimap2 is the close second best aligner for calling indels on this dataset, having only 0.1% lower recall and precision to strobealign.


![sv_calling_sim 001](https://user-images.githubusercontent.com/1714667/156549014-66d7b015-877e-48c2-a3b2-131e7a4a2db0.jpeg)
Figure 3. Recall precision and F-score for the aligners on 2x150 and 2x250 datasets from SIM3.

## Variant calling benchmark (GIAB)

We used Illumina paired-end reads from the [GIAB datasets](https://github.com/genome-in-a-bottle/giab_data_indexes) HG004 (Mother) with the 2x150bp reads (subsampled to ~26x coverage; using only the reads in `140818_D00360_0047_BHA66FADXX/Project_RM8392`) and 2x250bp reads (~17x coverage). We aligned the reads to hg38 without alternative haplotypes as proposed [here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use). We used 16 cores for all aligners. We obtain the "true" SNVs and INDELs from the GIAB gold standard predictions formed from several sequencing technologies. They are provided [here](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/).

Results are shown for SNVs and indels separately in Figure 4. For SNVs, predictions with strobealign as the aligner have the highest [F-score](https://en.wikipedia.org/wiki/F-score) of all benchmarked aligners on both datasets. Strobealign's alignments yield the highest precision at the cost of a slightly lower recall. As for indels, predictions have a low recall, precision, and F-score with all aligners. This may be because we benchmarked against all gold standard SVs for HG004 that were not SNVs (see below method for the evaluation). Overall, predictions using Bowtie2 are the most desirable on these datasets.


![sv_calling 001](https://user-images.githubusercontent.com/1714667/156128690-266ccfad-14c0-48a7-8fe7-de794d98cc65.jpeg)
Figure 4. Recall precision and F-score for the aligners on 2x150 and 2x250 datasets from HG004.

## Runtime 

For the four larger datasets above we show the runtime of aligners using 16 threads in Figure 5. The two SIM3 datasets are denoted SIM150 and SIM250, and the two GIAB datasets are denoted BIO150 and BIO250. Urmap was excluded from the timing benchmark because we can only get it to run with 1 core on our server as reported [here](https://github.com/rcedgar/urmap/issues/8). Strobealign is the fastest aligner across datasets. While urmap could be faster (based on the singlethreaded benchmarks), strobealign has substaintially better accuracy and downstream SV calling statistics (as seen in previous sections).

![runtime_sv](https://user-images.githubusercontent.com/1714667/161267752-8ff8d0ca-512d-46bc-a5f7-c3548dc603a3.png)
Figure 5. Runtime of aligners using 16 threads on two simulaed and two biological datasets of about 20-30x coverage of a human genome.


## Variant calling benchmark method

For the results, we ran

```
bcftools mpileup -O z --fasta-ref ref aligned.bam > aligned.vcf.gz
bcftools call -v -c -O v aligned.vcf.gz > aligned.variants.vcf.gz

# Split into SNP and INDELS
grep -v -E -e "INDEL;" aligned.variants.vcf.gz > aligned.variants.SNV.vcf
grep "#"  aligned.variants.vcf.gz > aligned.variants.INDEL.vcf
grep -E -e "INDEL;" aligned.variants.vcf.gz >> aligned.variants.INDEL.vcf

# Separate GIAB SNVs and INDELS
shell('zgrep "#" true.variants.vcf > true.variants.SNV.vcf')
shell('zgrep -P  "\t[ACGT]\t[ACGT]\t" true.variants.vcf >> true.variants.SNV.vcf')
shell('zgrep -v -P  "\t[ACGT]\t[ACGT]\t" true.variants.vcf > true.variants.INDEL.vcf')

for type in SNV INDEL
do
	bcftools sort -Oz aligned.variants.$type.vcf.gz -o aligned.variants.sorted.$type.vcf.gz
	bcftools index aligned.variants.sorted.$type.vcf.gz
	bcftools isec --nfiles 2 -O u true_variants.sorted.$type.vcf.gz  aligned.variants.sorted.$type.vcf -p out_$type
done
```


CREDITS
----------------

Sahlin, K. Strobealign: flexible seed size enables ultra-fast and accurate read alignment. Genome Biol 23, 260 (2022). [https://doi.org/10.1186/s13059-022-02831-7](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02831-7)


VERSION INFO
---------------

See [release page](https://github.com/ksahlin/strobealign/releases)


LICENCE
----------------

MIT license, see [LICENSE](https://github.com/ksahlin/strobealign/blob/main/LICENSE).

