# StrobeAlign
Aligns short reads using strobemers

# INSTALLATION

```
git clone https://github.com/ksahlin/StrobeAlign
cd StrobeAlign
g++ -std=c++11 main.cpp source/index.cpp source/edlib.cpp source/ksw2_extz2_sse.c -o StrobeAlign -O3 -mavx2
```


# USAGE

```
StrobeAlign -o <output.sam> ref.fa reads.fa [-k 22 -s 18 -f 0.0002]
```