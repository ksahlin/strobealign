cmake_minimum_required(VERSION 3.16)

project(strobealign VERSION 0.12.0)
include(FetchContent)

option(ENABLE_AVX "Enable AVX2 support" OFF)
option(PYTHON_BINDINGS "Build Python bindings" OFF)

find_package(ZLIB)
find_package(Threads)
find_package(OpenMP)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: RelWithDebInfo Debug Release" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "RelWithDebInfo" "Debug" "Release")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
add_compile_options(-Wall -Wextra -Werror=maybe-uninitialized)

add_subdirectory(ext/zstr)

# Obtain version from Git or fall back to PROJECT_VERSION if not building
# from a Git repository
add_custom_target(version
  ${CMAKE_COMMAND}
  -D CONFIGIN="${PROJECT_SOURCE_DIR}/src/version.hpp.in"
  -D CONFIGOUT="${PROJECT_BINARY_DIR}/version.hpp"
  -D DEFAULT_VERSION="${PROJECT_VERSION}"
  -P ${CMAKE_SOURCE_DIR}/GitVersion.cmake
)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/buildconfig.hpp.in"
  "${PROJECT_BINARY_DIR}/buildconfig.hpp"
)

add_library(salib STATIC ${SOURCES}
  src/refs.cpp
  src/fastq.cpp
  src/cmdline.cpp
  src/index.cpp
  src/indexparameters.cpp
  src/sam.cpp
  src/paf.cpp
  src/pc.cpp
  src/aln.cpp
  src/cigar.cpp
  src/aligner.cpp
  src/nam.cpp
  src/randstrobes.cpp
  src/readlen.cpp
  src/version.cpp
  src/io.cpp
  src/insertsizedistribution.cpp
  ext/xxhash.c
  ext/ssw/ssw_cpp.cpp
  ext/ssw/ssw.c
)
target_include_directories(salib PUBLIC src/ ext/ ${PROJECT_BINARY_DIR})
target_link_libraries(salib PUBLIC ZLIB::ZLIB Threads::Threads zstr::zstr)
IF(ENABLE_AVX)
  target_compile_options(salib PUBLIC "-mavx2")
ENDIF()
add_dependencies(salib version)

add_executable(strobealign src/main.cpp)
target_link_libraries(strobealign PUBLIC salib)
if(NOT PYTHON_BINDINGS)
  install(TARGETS strobealign DESTINATION bin)
endif()
add_executable(test-strobealign
  tests/tests.cpp
  tests/test_input.cpp
  tests/test_refs.cpp
  tests/test_sam.cpp
  tests/test_aligner.cpp
  tests/test_cigar.cpp
)
target_link_libraries(test-strobealign salib)
target_include_directories(test-strobealign PUBLIC src/ ext/ ${PROJECT_BINARY_DIR})

add_executable(dumpstrobes src/dumpstrobes.cpp)
target_link_libraries(dumpstrobes salib)
target_include_directories(dumpstrobes PUBLIC src/ ext/ ${PROJECT_BINARY_DIR})

# Currently unused code for which we only want to ensure that it still compiles
add_executable(unused-code src/unused.cpp)
target_link_libraries(unused-code salib)
target_include_directories(unused-code PUBLIC src/ ext/ ${PROJECT_BINARY_DIR})

if(PYTHON_BINDINGS)
  add_subdirectory(src/python)
endif()
