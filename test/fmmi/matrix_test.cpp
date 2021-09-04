#include <doctest.h>

#include "fmmi/matrix.hpp"

TEST_CASE("matrix test")
{
    matrix<4, 4> mx{
        1.0, 2.0, 3.0, 4.0,
        5.0, 6.0, 7.0, 8.0,
        9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0,
    };

    CHECK_EQ(mx(0, 0), 1.0);
    CHECK_EQ(mx(0, 1), 5.0);
    CHECK_EQ(mx(0, 2), 9.0);
    CHECK_EQ(mx(0, 3), 13.0);
    CHECK_EQ(mx(1, 0), 2.0);
    CHECK_EQ(mx(1, 1), 6.0);
    CHECK_EQ(mx(1, 2), 10.0);
    CHECK_EQ(mx(1, 3), 14.0);
}
