#include <doctest.h>

#include "fmmi/matrix.hpp"

using namespace fmmi;

TEST_CASE("matrix::matrix()")
{
    matrix<5, 4> mx{
        1.0, 2.0, 3.0, 4.0,
        5.0, 6.0, 7.0, 8.0,
        9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0,
        17.0, 18.0, 19.0, 20.0,
    };

    CHECK_EQ(mx(0, 0), 1.0);
    CHECK_EQ(mx(0, 1), 2.0);
    CHECK_EQ(mx(0, 2), 3.0);
    CHECK_EQ(mx(0, 3), 4.0);
    CHECK_EQ(mx(1, 0), 5.0);
    CHECK_EQ(mx(1, 1), 6.0);
    CHECK_EQ(mx(1, 2), 7.0);
    CHECK_EQ(mx(1, 3), 8.0);
    CHECK_EQ(mx(2, 0), 9.0);
    CHECK_EQ(mx(2, 1), 10.0);
    CHECK_EQ(mx(2, 2), 11.0);
    CHECK_EQ(mx(2, 3), 12.0);
    CHECK_EQ(mx(3, 0), 13.0);
    CHECK_EQ(mx(3, 1), 14.0);
    CHECK_EQ(mx(3, 2), 15.0);
    CHECK_EQ(mx(3, 3), 16.0);
    CHECK_EQ(mx(4, 0), 17.0);
    CHECK_EQ(mx(4, 1), 18.0);
    CHECK_EQ(mx(4, 2), 19.0);
    CHECK_EQ(mx(4, 3), 20.0);
}


TEST_CASE("matrix::operator==()")
{
    matrix<2, 3> mx1{
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
    };

    matrix<2, 3> mx2{
        1.0, 2.0, 3.0,
        5.0, 6.0, 7.0,
    };

    CHECK_EQ(mx1, mx1);

    CHECK_NE(mx1, mx2);

    mx2(1, 0) = 4.0;
    mx2(1, 1) = 5.0;
    mx2(1, 2) = 6.0;

    CHECK_EQ(mx1, mx2);
}


TEST_CASE("add()")
{
    matrix<2, 3> mx1{
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
    };

    matrix<2, 3> mx2{
        10.0, 11.0, 12.0,
        20.0, 21.0, 22.0,
    };

    matrix<2, 3> mx3;

    add(mx1, mx2, mx3);

    CHECK_EQ(mx3(0, 0), 11.0);
    CHECK_EQ(mx3(0, 1), 13.0);
    CHECK_EQ(mx3(0, 2), 15.0);
    CHECK_EQ(mx3(1, 0), 24.0);
    CHECK_EQ(mx3(1, 1), 26.0);
    CHECK_EQ(mx3(1, 2), 28.0);
}


TEST_CASE("mul()")
{
    matrix<2, 3> mx1{
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
    };

    matrix<3, 2> mx2{
        10.0, 11.0,
        20.0, 21.0,
        30.0, 31.0,
    };

    matrix<2, 2> mx3;

    mul(mx1, mx2, mx3);

    CHECK_EQ(mx3(0, 0), 140.0);
    CHECK_EQ(mx3(0, 1), 146.0);
    CHECK_EQ(mx3(1, 0), 320.0);
    CHECK_EQ(mx3(1, 1), 335.0);
}


TEST_CASE("matrix::partition()")
{
    matrix<5, 4> mx{
        1.0, 2.0, 3.0, 4.0,
        5.0, 6.0, 7.0, 8.0,
        9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0,
        17.0, 18.0, 19.0, 20.0,
    };

    matrix<2, 3> mx11{
        1.0, 2.0, 3.0,
        5.0, 6.0, 7.0,
    };

    matrix<2, 1> mx12{
        4.0,
        8.0,
    };

    matrix<3, 1> mx21{
        9.0,
        13.0,
        17.0,
    };

    matrix<3, 3> mx22{
        10.0, 11.0, 12.0,
        14.0, 15.0, 16.0,
        18.0, 19.0, 20.0,
    };

    auto part11 = mx.partition<0, 0, 2, 3>();
    auto part12 = mx.partition<0, 3, 2, 1>();
    auto part21 = mx.partition<2, 0, 3, 1>();
    auto part22 = mx.partition<2, 1, 3, 3>();

    CHECK_EQ(part11, mx11);
    CHECK_EQ(part12, mx12);
    CHECK_EQ(part21, mx21);
    CHECK_EQ(part22, mx22);
}
