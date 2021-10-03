include_directories(${CMAKE_SOURCE_DIR}/3rdparty/catch2)

set(TEST_SOURCES
    test.cpp
    fmmi/smatrix_test.cpp
    fmmi/dmatrix_test.cpp
    fmmi/fmmi_test.cpp
)

add_executable(test ${TEST_SOURCES})
target_link_libraries(test fmmi)
