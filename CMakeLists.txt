cmake_minimum_required(VERSION 3.16)

project(strobealign VERSION 0.7.1)

configure_file(
    "${PROJECT_SOURCE_DIR}/src/version.hpp.in"
    "${PROJECT_BINARY_DIR}/version.hpp"
)
option(ENABLE_AVX "Enable AVX2 support" OFF)

find_package(ZLIB)
find_package(Threads)
find_package(OpenMP)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_FLAGS "-O3 ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "-O3 ${CMAKE_CXX_FLAGS}")

set(SOURCES
src/main.cpp
src/cmdline.cpp
src/index.cpp
src/pc.cpp
src/aln.cpp
ext/xxhash.c
ext/ssw_cpp.cpp
ext/ssw.c
)

add_executable(strobealign ${SOURCES})

target_link_libraries(strobealign PUBLIC ZLIB::ZLIB Threads::Threads)
target_include_directories(strobealign PUBLIC src/ ext/ ${PROJECT_BINARY_DIR})
install(TARGETS strobealign DESTINATION bin)

IF(ENABLE_AVX)
  target_compile_options(strobealign PUBLIC "-mavx2")
ENDIF()


