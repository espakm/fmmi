include_directories(${CMAKE_SOURCE_DIR}/3rdparty/catch2)

set(TEST_SOURCES
    fmmi/common_test.cpp
    fmmi/dmatrix_test.cpp
    fmmi/smatrix_test.cpp
    fmmi/dmatrix_rec_test.cpp
    fmmi/smatrix_rec_test.cpp
    test.cpp
)

add_executable(test ${TEST_SOURCES})
target_link_libraries(test fmmi)
target_compile_definitions(test PRIVATE
    DMATRIX_REC_TEST_MAX_MATRIX_SIZE=${DMATRIX_REC_TEST_MAX_MATRIX_SIZE}
    SMATRIX_REC_TEST_MAX_MATRIX_SIZE=${SMATRIX_REC_TEST_MAX_MATRIX_SIZE}
)
