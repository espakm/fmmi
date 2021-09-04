include_directories(${CMAKE_SOURCE_DIR}/3rdparty/doctest)

set(TEST_SOURCES
    doctest.cpp
    fmmi/fmmi_test.cpp
    fmmi/matrix_test.cpp
)

add_executable(doctest ${TEST_SOURCES})
