# Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

cmake_minimum_required(VERSION 3.0)

set(clang_path "/usr/bin/clang++")

if(CMAKE_BUILD_TYPE MATCHES debug AND EXISTS "${clang_path}")
    set(CMAKE_CXX_COMPILER "${clang_path}")
else()
    set(CMAKE_CXX_COMPILER "/usr/bin/g++")
endif()

message("Using ${CMAKE_CXX_COMPILER}")

project(su2-hmc CXX C)

add_definitions("-march=native")
add_definitions("-std=c++11")
#add_definitions("-Werror")
add_definitions("-Wall")
#add_definitions(-DBOOST_FILESYSTEM_NO_DEPRECATED)

###############################################################################
#                                Compile flags                                #
###############################################################################

if(CMAKE_BUILD_TYPE MATCHES release)
    add_definitions("-Werror")
endif()

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-fdiagnostics-color has_diagnostics_color)
if(has_diagnostics_color)
    add_definitions("-fdiagnostics-color=always")
endif()

check_cxx_compiler_flag(-Wdocumentation has_W_documentation)
if(has_W_documentation)
    add_definitions("-Wdocumentation")
endif()

###############################################################################
#                                  Libraries                                  #
###############################################################################

find_package(OpenMP)
if(OPENMP_FOUND AND CMAKE_BUILD_TYPE MATCHES release)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    add_definitions("-Wno-unknown-pragmas")
endif()

find_package(Boost REQUIRED COMPONENTS)

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

###############################################################################
#                                 Executables                                 #
###############################################################################

add_executable(su2-hmc

    configuration.cpp
    hybrid-monte-carlo.cpp
    main.cpp
    pauli-matrices.cpp
    sanity-checks.cpp
    su2matrix.cpp

    )

target_link_libraries(su2-hmc ${Boost_LIBS})

add_subdirectory(tests)
