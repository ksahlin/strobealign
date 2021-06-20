StrobeAlign
==============

Short-read aligner using syncmer-thinned strobemers. This implementation is currently proof-of-concept. The default parameter setting is tailored for Illumina reads of length 150-500. For longer reads, we need chaining of seeds instead of the current merging of overlapping seeds as described in the [preprint]().


INSTALLATION
----------------

```
git clone https://github.com/ksahlin/StrobeAlign
cd StrobeAlign
g++ -std=c++11 main.cpp source/index.cpp source/ksw2_extz2_sse.c -o StrobeMap -O3 -mavx2
```


USAGE
-------

```
StrobeAlign -o <output.sam> ref.fa reads.fa [-k 22 -s 18 -f 0.0002]
```


TODO
-------

1. Add option to separate build index and perform alignment in separate steps.
2. Implement multi-threading.
3. Allow reads in fastq format
4. (Eventually) Consider mate inference for Illumina Paired-End mapping.


<!-- CREDITS
----------------


1. Kristoffer Sahlin 2021. "" [preprint available here]().
 -->


LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/uLTRA/blob/master/LICENCE.txt).


