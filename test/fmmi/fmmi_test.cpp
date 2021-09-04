#include <doctest.h>

#include "fmmi/fmmi.hpp"

using namespace fmmi;

TEST_CASE("fmmi test")
{
    matrix<2, 2> mx1{
        -2.0, 3.0,
        -7.0, 9.0,
    };

    matrix<2, 2> mx2{
        5.0, -8.0,
        4.0, -6.0,
    };

    matrix<2, 2> mx3;

    mul(mx1, mx2, mx3);

    matrix<2, 2> mx4;

    mul_fast(mx1, mx2, mx4);

    CHECK_EQ(mx3, mx4);
}
