cmake_minimum_required(VERSION 3.16)
project(qsim VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
        "Choose the type of build: Debug Release RelWithDebInfo MinSizeRel"
        FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")



#

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Debug flags
    set(CMAKE_CXX_FLAGS_DEBUG
        "-O0 -g -Wall -Wextra -Wpedantic -DDEBUG"
        CACHE STRING "Debug build flags" FORCE)

    # Release flags - full optimization
    set(CMAKE_CXX_FLAGS_RELEASE
        "-O3 -march=native -mtune=native -DNDEBUG"
        CACHE STRING "Release build flags" FORCE)

    # Release with debug info
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
        "-O3 -march=native -mtune=native -g -DNDEBUG"
        CACHE STRING "RelWithDebInfo build flags" FORCE)

    # Minimum size release
    set(CMAKE_CXX_FLAGS_MINSIZEREL
        "-Os -DNDEBUG"
        CACHE STRING "MinSizeRel build flags" FORCE)

    # Enable Link Time Optimization for Release builds
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)

elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # Debug flags
    set(CMAKE_CXX_FLAGS_DEBUG
        "/Od /Zi /MDd /DDEBUG"
        CACHE STRING "Debug build flags" FORCE)

    # Release flags
    set(CMAKE_CXX_FLAGS_RELEASE
        "/O2 /MD /DNDEBUG /arch:AVX2"
        CACHE STRING "Release build flags" FORCE)

    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
        "/O2 /MD /Zi /DNDEBUG /arch:AVX2"
        CACHE STRING "RelWithDebInfo build flags" FORCE)

    # Enable LTO for Release
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

#




option(QSIM_BUILD_TESTS "Build tests with Catch2" ON)
option(QSIM_BUILD_BENCHMARKS "Build benchmarks with Google Benchmark" OFF)
option(QSIM_BUILD_BINDINGS "Build Python bindings with pybind11" ON)

if(QSIM_BUILD_BINDINGS)
    find_package(pybind11 REQUIRED)
endif()

if(QSIM_BUILD_TESTS)



    # get Catch2
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.4.0
    )
    FetchContent_MakeAvailable(Catch2)




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
