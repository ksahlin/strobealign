cmake_minimum_required(VERSION 3.16)

project(StrobeAlign VERSION 0.5)

find_package(ZLIB)

find_package(OpenMP)

add_executable(StrobeAlign
   main.cpp
   source/index.cpp
   source/ksw2_extz2_sse.c
   source/xxhash.c
   source/ssw_cpp.cpp
   source/ssw.c
)

target_link_libraries(StrobeAlign PUBLIC ZLIB::ZLIB OpenMP::OpenMP_CXX)

target_compile_features(StrobeAlign PUBLIC cxx_std_14)