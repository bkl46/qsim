cmake_minimum_required(VERSION 3.16)
project(qsim VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(QSIM_BUILD_TESTS "Build tests with Catch2" OFF)
option(QSIM_BUILD_BENCHMARKS "Build benchmarks with Google Benchmark" OFF)
option(QSIM_BUILD_BINDINGS "Build Python bindings with pybind11" OFF)

if(QSIM_BUILD_BINDINGS)
    find_package(pybind11 REQUIRED)
endif()

if(QSIM_BUILD_TESTS)
    find_package(Catch2 REQUIRED)
endif()

if(QSIM_BUILD_BENCHMARKS)
    find_package(benchmark REQUIRED)
endif()


# goto source
add_subdirectory(src)



if(QSIM_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(QSIM_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if(QSIM_BUILD_BINDINGS)
    add_subdirectory(bindings)
endif()
