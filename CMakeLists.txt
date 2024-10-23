cmake_minimum_required(VERSION 3.16)

project(strobealign VERSION 0.14.0)
include(FetchContent)
include(ExternalProject)

option(ENABLE_AVX "Enable AVX2 support" OFF)
option(PYTHON_BINDINGS "Build Python bindings" OFF)
option(TRACE "Highly verbose debugging output" OFF)

find_package(ZLIB)
find_package(Threads)
find_package(OpenMP)
find_package(PkgConfig REQUIRED)

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

set(ISAL "system" CACHE STRING "Where to find ISA-L: system (uses pkg-config), download (at build time)")
set_property(CACHE ISAL PROPERTY STRINGS "system" "download")

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
  src/iowrap.cpp
  ext/xxhash.c
  ext/ssw/ssw_cpp.cpp
  ext/ssw/ssw.c
)
if(ISAL STREQUAL "system")
  pkg_check_modules(ISAL IMPORTED_TARGET GLOBAL libisal>=2.30.0)
  if (NOT ISAL_FOUND)
    message(FATAL_ERROR "libisal (ISA-L) could not be found. You have two options:\n"
        "- Install the library using your package manager (try 'sudo apt install "
        "libisal-dev' or 'brew install isa-l')\n"
        "- or add '-DISAL=download' to "
        "the CMake options in order to download ISA-L at build time."
      )
  endif()
  add_library(ISAL ALIAS PkgConfig::ISAL)
elseif(ISAL STREQUAL "download")
  find_program(NASM NAMES nasm)
  if (NASM STREQUAL NASM-NOTFOUND)
    message(FATAL_ERROR "NASM is required to build ISA-L. Try 'sudo apt "
        "install nasm' or 'brew install nasm'"
    )
  endif()
  set(ISAL_ROOT ${CMAKE_BINARY_DIR}/ISAL)
  ExternalProject_Add(
    isal_external
    GIT_REPOSITORY https://github.com/intel/isa-l
    GIT_TAG v2.30.0
    BUILD_COMMAND make -f Makefile.unx
    CONFIGURE_COMMAND ""
    INSTALL_DIR ${ISAL_ROOT}
    INSTALL_COMMAND make -f Makefile.unx install prefix=<INSTALL_DIR>
    BUILD_IN_SOURCE ON
  )
  target_include_directories(salib PUBLIC "${ISAL_ROOT}/include")
  add_library(ISAL STATIC IMPORTED)
  set_target_properties(ISAL PROPERTIES IMPORTED_LOCATION ${ISAL_ROOT}/lib/libisal.a)
  add_dependencies(salib isal_external)
#elseif(ISAL STREQUAL "off")
#  message(FATAL_ERROR "Building without ISA-L support is currently not supported")
endif()

target_include_directories(salib PUBLIC src/ ext/ ${PROJECT_BINARY_DIR})
target_link_libraries(salib PUBLIC ZLIB::ZLIB Threads::Threads zstr::zstr ISAL)
IF(ENABLE_AVX)
  target_compile_options(salib PUBLIC "-mavx2")
ENDIF()
if (TRACE)
   target_compile_definitions(salib PUBLIC "TRACE")
endif()

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
  tests/test_randstrobes.cpp
)
target_link_libraries(test-strobealign salib)
target_include_directories(test-strobealign PUBLIC src/ ext/ ${PROJECT_BINARY_DIR})

add_executable(dumpstrobes src/dumpstrobes.cpp)
target_link_libraries(dumpstrobes salib)
target_include_directories(dumpstrobes PUBLIC src/ ext/ ${PROJECT_BINARY_DIR})

if(PYTHON_BINDINGS)
  add_subdirectory(src/python)
endif()
