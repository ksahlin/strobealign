cmake_minimum_required(VERSION 3.19)
project(StrobeAlign)

set(CMAKE_CXX_STANDARD 14)

add_executable(StrobeAlign main.cpp source/index.cpp source/index.hpp source/robin_hood.h source/edlib.cpp source/edlib.h source/ksw2.h)