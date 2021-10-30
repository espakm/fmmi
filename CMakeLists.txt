cmake_minimum_required(VERSION 3.13)

project(fmmi LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug")

set(DMATRIX_REC_TEST_MAX_MATRIX_SIZE "16" CACHE STRING
    "Maximum matrix size for FMMI tests with dynamic sized matrices")
set_property(CACHE DMATRIX_REC_TEST_MAX_MATRIX_SIZE PROPERTY STRINGS
    "2" "4" "8" "16" "32" "64" "128" "256" "512" "1024" "2048" "4096")

set(SMATRIX_REC_TEST_MAX_MATRIX_SIZE "16" CACHE STRING
    "Maximum matrix size for FMMI tests with static sized matrices")
set_property(CACHE SMATRIX_REC_TEST_MAX_MATRIX_SIZE PROPERTY STRINGS
    "2" "4" "8" "16" "32" "64" "128" "256" "512" "1024" "2048" "4096")

add_subdirectory(src)
add_subdirectory(test)
