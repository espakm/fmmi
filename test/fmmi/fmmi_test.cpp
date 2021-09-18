#include <catch.hpp>

#include "fmmi/fmmi.hpp"

using namespace fmmi;

TEST_CASE("fmmi mul_fast")
{
    matrix<4, 6> mx1{
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0,
    };
    matrix<4, 6> mx2{
        -1.0, -2.0, -3.0, -4.0, -5.0, -6.0,
        -7.0, -8.0, -9.0, -10.0, -11.0, -12.0,
        -13.0, -14.0, -15.0, -16.0, -17.0, -18.0,
        -19.0, -20.0, -21.0, -22.0, -23.0, -24.0,
    };

    SECTION("2x2 times 2x2")
    {
        const auto a = mx1.partition<0, 0, 2, 2>();
        const auto b = mx2.partition<0, 0, 2, 2>();

        matrix<2, 2> c;
        mul(a, b, c);

        matrix<2, 2> d;
        mul_fast(a, b, d);

        CHECK(c == d);
    }

    SECTION("4x4 times 4x4")
    {
        const auto a = mx1.partition<0, 0, 4, 4>();
        const auto b = mx2.partition<0, 0, 4, 4>();

        matrix<4, 4> c;
        mul(a, b, c);

        matrix<4, 4> d;
        mul_fast(a, b, d);

        CHECK(c == d);
    }

    SECTION("3x3 times 3x3")
    {
        const auto a = mx1.partition<0, 0, 3, 3>();
        const auto b = mx2.partition<0, 0, 3, 3>();

        matrix<3, 3> c;
        mul(a, b, c);

        matrix<3, 3> d;
        mul_fast(a, b, d);

        CHECK(c == d);
    }
}
