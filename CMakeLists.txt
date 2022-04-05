cmake_minimum_required(VERSION 3.16)

project(strobealign VERSION 0.7)

option(ENABLE_AVX "Enable AVX2 support" ON)

find_package(ZLIB)
find_package(Threads)
find_package(OpenMP)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_FLAGS "-O3 ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "-O3 ${CMAKE_CXX_FLAGS}")

set(SOURCES
main.cpp
source/index.cpp
source/ksw2_extz2_sse.c
source/xxhash.c
source/ssw_cpp.cpp
source/ssw.c
source/pc.cpp
source/aln.cpp
)

add_executable(strobealign ${SOURCES})

target_link_libraries(strobealign PUBLIC ZLIB::ZLIB Threads::Threads)

install(TARGETS strobealign DESTINATION bin)

IF(ENABLE_AVX)
  target_compile_options(strobealign PUBLIC "-mavx2")
ENDIF()


