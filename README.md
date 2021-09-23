StrobeAlign
==============

Short-read aligner using syncmer-thinned strobemers. The implementation is currently proof-of-concept. The default parameter setting is tailored for Illumina reads of length 150-500. For longer reads, we need chaining of seeds instead of the current approach, which is merging of overlapping seeds as described in the [preprint](https://doi.org/10.1101/2021.06.18.449070).


INSTALLATION
----------------

```
git clone https://github.com/ksahlin/StrobeAlign
cd StrobeAlign
g++ -std=c++14 main.cpp source/index.cpp source/ksw2_extz2_sse.c -lz -fopenmp -o StrobeAlign -O3 -mavx2
```


USAGE
-------

For alignment to SAM file:

```
StrobeAlign [-k 22 -s 18 -f 0.0002] -o <output.sam> ref.fa reads.fa 
```

For mapping to PAF file (option -x):

```
StrobeAlign [-k 22 -s 18 -f 0.0002] -x -o <output.sam> ref.fa reads.fa 
```

TODO
-------

1. Add option to separate build index and perform alignment in separate steps.
2. Implement multi-threading.
3. Allow reads in fastq format
4. (Eventually) Consider mate inference for Illumina Paired-End mapping.


CREDITS
----------------

Kristoffer Sahlin. Faster short-read mapping with strobemer seeds in syncmer space. bioRxiv, 2021. doi:10.1101/2021.06.18.449070. Preprint available [here](https://doi.org/10.1101/2021.06.18.449070).


LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/uLTRA/blob/master/LICENCE.txt).


