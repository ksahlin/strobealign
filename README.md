# StrobeAlign
Aligns short reads using strobemers

# INSTALLATION

```
git clone https://github.com/ksahlin/StrobeAlign
cd StrobeAlign
g++ -std=c++11 main.cpp source/index.cpp source/edlib.cpp source/ksw2_extz2_sse.c -o StrobeMap -O3 -mavx2
```
