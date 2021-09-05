#include <doctest.h>

#include "fmmi/fmmi.hpp"

using namespace fmmi;

TEST_CASE("fmmi mul_fast")
{
    matrix<2, 2> mx1{
        -2.0, 3.0,
        -7.0, 9.0,
    };

    matrix<2, 2> mx2{
        5.0, -8.0,
        4.0, -6.0,
    };

    SUBCASE("2x2 times 2x2")
    {
        matrix<2, 2> mx3;
        mul(mx1, mx2, mx3);

        matrix<2, 2> mx4;
        mul_fast(mx1, mx2, mx4);

        CHECK_EQ(mx3, mx4);
    }

    matrix<4, 4> mx5{
        1.0, 2.0, 3.0, 4.0,
        7.0, 8.0, 9.0, 10.0,
        13.0, 14.0, 15.0, 16.0,
        19.0, 20.0, 21.0, 22.0,
    };

    matrix<4, 4> mx6{
        3.0, 4.0, 5.0, 6.0,
        9.0, 10.0, 11.0, 12.0,
        15.0, 16.0, 17.0, 18.0,
        21.0, 22.0, 23.0, 24.0,
    };

    SUBCASE("4x4 times 4x4")
    {
        matrix<4, 4> mx3;
        mul(mx5, mx6, mx3);

        matrix<4, 4> mx4;
        mul_fast(mx5, mx6, mx4);

        CHECK_EQ(mx3, mx4);
    }
}
