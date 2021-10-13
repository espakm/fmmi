include_directories(${CMAKE_SOURCE_DIR}/3rdparty/catch2)

set(TEST_SOURCES
    fmmi/dmatrix_test.cpp
    fmmi/smatrix_test.cpp
    fmmi/dmatrix_fmmi_test.cpp
    fmmi/smatrix_fmmi_test.cpp
    test.cpp
)

add_executable(test ${TEST_SOURCES})
target_link_libraries(test fmmi)
target_compile_definitions(test PRIVATE
    DMATRIX_FMMI_TEST_MAX_MATRIX_SIZE=${DMATRIX_FMMI_TEST_MAX_MATRIX_SIZE}
    SMATRIX_FMMI_TEST_MAX_MATRIX_SIZE=${SMATRIX_FMMI_TEST_MAX_MATRIX_SIZE}
)
