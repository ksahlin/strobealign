StrobeAlign
==============

Strobealign is a fast short-read aligner. It achieves the speedup by using a dynamic seed size obtained from syncmer-thinned strobemers. Strobealign is multithreaded, implements alignment (SAM) and mapping (PAF), and has high accuracy for reads of lengths between 100-300bp and insert sizes up to roughly 1000bp. A somewhat outdated preprint describing v0.0.3 is available [here](https://doi.org/10.1101/2021.06.18.449070).

Results for version 0.2 below when alinging PE reads simulated at various variation rates (SIM1-3) to hg38. Solid lines (alignment) are what matters in practice, mapping mode included for completion.

<img width="2374" alt="accuracy" src="https://user-images.githubusercontent.com/1714667/147755184-f92fdf90-250f-4768-88f9-9ff6c180de2f.png">

<img width="2318" alt="runtime" src="https://user-images.githubusercontent.com/1714667/147755453-5b230acf-fb1a-4b2e-b487-53ebb084e273.png">

INSTALLATION
----------------

You can acquire precompiled binaries for Linux and Mac OSx from the [release page](https://github.com/ksahlin/StrobeAlign/releases) compiled with `-O3 -mavx2`. 

It has been [reported](https://github.com/ksahlin/StrobeAlign/issues/6) that `strobealign` is even faster if compliled with flag `-march=skylake-avx512` for avx512 supported processors.

If you want to compile from the source, you need to have a newer `g++` and [zlib](https://zlib.net/) installed. Then do the following:

```
git clone https://github.com/ksahlin/StrobeAlign
cd StrobeAlign
# Needs a newer g++ version. Tested with version 8 and upwards.
g++ -std=c++14 main.cpp source/index.cpp source/xxhash.c source/ksw2_extz2_sse.c source/ssw_cpp.cpp source/ssw.c -lz -fopenmp -o StrobeAlign -O3 -mavx2
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
g++ -std=c++14 -I/path/to/zlib/include -L/path/to/zlib/lib main.cpp source/index.cpp source/xxhash.c source/ksw2_extz2_sse.c source/ssw_cpp.cpp source/ssw.c -lz -fopenmp -o StrobeAlign -O3 -mavx2
``` 


USAGE
-------

Strobealign v0.1 and up comes with a parameter `-r read_length` that sets suitable seed parameters for the rough read length. Specifically, it sets parameters `-k`, `-l` and `-u`. If not specified, it defaults to 150. The value of `r` does not have to match the exact read length.

For alignment to SAM file:

```
StrobeAlign -r <read_length> -o <output.sam> ref.fa reads.fa 
```

For mapping to PAF file (option -x):

```
StrobeAlign -r <read_length> -x -o <output.sam> ref.fa reads.fa 
```

TODO
-------

1. Add option to separate build index and perform alignment in separate steps.


CREDITS
----------------

Kristoffer Sahlin. Faster short-read mapping with strobemer seeds in syncmer space. bioRxiv, 2021. doi:10.1101/2021.06.18.449070. Preprint available [here](https://doi.org/10.1101/2021.06.18.449070).

VERSION INFO
---------------

### Version 0.2.1
1. Inroduced a max seed size contraint when sampling seeds, only active in few regions where syncmers are sparsely sampled.
2. Parameter `-r` can now take any integer value. 

### Version 0.2

Important bugfix [1](https://github.com/ksahlin/StrobeAlign/commit/39c8c45afd6d9b35ea55da5744ae12b810fc8086) and added ssw for rescue alignment [2](https://github.com/ksahlin/StrobeAlign/commit/8282043129f2f8fcdab9106c6e5f3c5777e4220c) since ksw is only for extension. These fixes improve both accuracy and speed in paired-end alignment mode further.

### Version 0.1

Major update to algorithm. See release page.

### Version 0.0.3.2

1. Takes care of negative alignment coordinate bug.
2. Minimal value for repetitive seed filtering implemented. Previously, the top fraction of `-f (0.0002)` seeds was filtered regardless of how repetitive they were. Assume `-f` filtered everything above `X` occurrences. The new version filters seeds with occurrences over `max(X, 30)`. This threshold is usually not active for hg38 as `X>40` for hg38. 

### Version 0.0.3.1

1. Bugfix. Takes care of segmentation fault bug in paired-end mapping mode (-x) when none of the reads have NAMs.

### Version 0.0.3

1. Implements a paired-end alignment mode.
2. Implements a rescue mode both in SE and PE alignment modes (described in preprint v2).
3. Changed to symmetrical strobemer hash values due to inversions (described in preprint v2).


### Version 0.0.2

1. Implements multi-threading.
2. Allow reads in fast[a/q] format and gzipped files through [kseqpp library](https://github.com/cartoonist/kseqpp).

### Version 0.0.1

The aligner used for the experiments presented in the first preprint (v1) on bioRxiv. Only single-threaded alignment and aligns reads as single reads (no PE mapping).

LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/uLTRA/blob/master/LICENCE.txt).


