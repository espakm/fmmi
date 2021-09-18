include_directories(${CMAKE_SOURCE_DIR}/3rdparty/catch2)

set(TEST_SOURCES
    test.cpp
    fmmi/fmmi_test.cpp
    fmmi/matrix_test.cpp
)

add_executable(test ${TEST_SOURCES})
target_link_libraries(test fmmi)
