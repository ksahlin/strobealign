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

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

add_compile_options(-Wall -Wextra)

add_library(salib STATIC ${SOURCES}
  src/refs.cpp
  src/fastq.cpp
  src/cmdline.cpp
  src/index.cpp
  src/sam.cpp
  src/pc.cpp
  src/aln.cpp
  src/randstrobes.cpp
  src/readlen.cpp
  ext/xxhash.c
  ext/ssw_cpp.cpp
  ext/ssw.c
)
target_include_directories(salib PUBLIC src/ ext/ ${PROJECT_BINARY_DIR})
target_link_libraries(salib PUBLIC ZLIB::ZLIB Threads::Threads)
IF(ENABLE_AVX)
  target_compile_options(salib PUBLIC "-mavx2")
ENDIF()

add_executable(strobealign src/main.cpp)
target_link_libraries(strobealign PUBLIC salib)
install(TARGETS strobealign DESTINATION bin)


add_executable(test-strobealign tests/tests.cpp)
target_link_libraries(test-strobealign salib)
target_include_directories(test-strobealign PUBLIC src/ ext/ ${PROJECT_BINARY_DIR})

# Currently unused code for which we only want to ensure that it still compiles
add_executable(unused-code src/unused.cpp)
target_link_libraries(unused-code salib)
target_include_directories(unused-code PUBLIC src/ ext/ ${PROJECT_BINARY_DIR})
