StrobeAlign
==============

Short-read aligner using syncmer-thinned strobemers. Strobealign implements PE and SE Illumina read alignment/mapping and is multithreaded. The default parameter setting is tailored for Illumina reads of lengths about 150-500nt. 

For reads shorter than 150nt, a lower value of parameter `-k` is recommentded (e.g. 15-17). 

StrobeAlign is currently not recommentded for long reads (>500nt) as significant implementation changes is needed. For long reads we need a different extention algorithm (chaining of seeds instead of the current approach described in the [preprint](https://doi.org/10.1101/2021.06.18.449070)) and split-mapping funcitionality.


INSTALLATION
----------------

You can acquire precompiled binaries for Linux and Mac OSx from [here](https://github.com/ksahlin/StrobeAlign/tree/main/bin). For example, for linux, simply do

```
wget https://github.com/ksahlin/StrobeAlign/tree/main/bin/Linux/StrobeAlign-v0.0.3.1
chmod +x StrobeAlign-v0.0.3.1
./StrobeAlign-v0.0.3.1  # test program
```

If you want to compile from the source, you need to have a newer `g++` and [zlib]](https://zlib.net/) installed. Then do the following:

```
git clone https://github.com/ksahlin/StrobeAlign
cd StrobeAlign
# Needs a newer g++ version. Tested with version 8 and upwards.
g++ -std=c++14 main.cpp source/index.cpp source/ksw2_extz2_sse.c -lz -fopenmp -o StrobeAlign -O3 -mavx2
```

## Common installation from source errors

If you have `zlib` installed, and the `zlib.h` file is in folder `/path/to/zlib/include` and the `libz.so` file in `/path/to/zlib/lib` but you get 

```
main.cpp:12:10: fatal error: zlib.h: No such file or directory
 #include <zlib.h>
          ^~~~~~~~
compilation terminated.
```

add `-I/path/to/zlib/include -L/path/to/zlib/lib` to the compilation, that is

```
g++ -std=c++14 -I/path/to/zlib/include -L/path/to/zlib/lib main.cpp source/index.cpp source/ksw2_extz2_sse.c -lz -fopenmp -o StrobeAlign -O3 -mavx2
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
2. Consider some smart PE-alignment mode (rescue mates etc).


CREDITS
----------------

Kristoffer Sahlin. Faster short-read mapping with strobemer seeds in syncmer space. bioRxiv, 2021. doi:10.1101/2021.06.18.449070. Preprint available [here](https://doi.org/10.1101/2021.06.18.449070).

VERSION INFO
---------------

### Version 0.0.3.1

1. Bugfix. Takes care of segmentation fault bug in paired-end mapping mode (-x) when none of the reads have NAMs.

### Version 0.0.3

1. Implements a paired-end alignment mode.
2. Implements a rescue mode both in SE and PE alignment modes (described in preprint v2).
3. Changed to symmetrical strobemer hashvalues due to inversions (described in preprint v2).


### Version 0.0.2

1. Implements multi-threading.
2. Allow reads in fast[a/q] format and gzipped files through [kseqpp library](https://github.com/cartoonist/kseqpp).

### Version 0.0.1

The aligner used for the experiments presented in the preprint (v1) on bioRxiv. Only single threaded alignment and aligns reads as single reads (no PE mapping).

LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/uLTRA/blob/master/LICENCE.txt).


