#include <catch.hpp>

#include "fmmi/matrix.hpp"

using namespace fmmi;

TEST_CASE("log_2")
{
    CHECK(log_2(1) == 0);
    CHECK(log_2(2) == 1);
    CHECK(log_2(3) == 1);
    CHECK(log_2(4) == 2);
    CHECK(log_2(5) == 2);
    CHECK(log_2(6) == 2);
    CHECK(log_2(7) == 2);
    CHECK(log_2(8) == 3);
    CHECK(log_2(9) == 3);
}


TEST_CASE("exp_2")
{
    CHECK(exp_2(0) == 1);
    CHECK(exp_2(1) == 2);
    CHECK(exp_2(2) == 4);
    CHECK(exp_2(3) == 8);
    CHECK(exp_2(4) == 16);
    CHECK(exp_2(5) == 32);
}


TEST_CASE("padded_size")
{
    CHECK(padded_size(1) == 1);
    CHECK(padded_size(2) == 2);
    CHECK(padded_size(3) == 4);
    CHECK(padded_size(4) == 4);
    CHECK(padded_size(5) == 8);
    CHECK(padded_size(6) == 8);
    CHECK(padded_size(7) == 8);
    CHECK(padded_size(8) == 8);
    CHECK(padded_size(9) == 16);

    CHECK(padded_size(3, 4) == 4 * 4);
    CHECK(padded_size(3, 5) == 8 * 8);
    CHECK(padded_size(5, 3) == 8 * 8);
}


TEST_CASE("matrix::matrix()")
{
    matrix<5, 4> mx{
        1.0, 2.0, 3.0, 4.0,
        5.0, 6.0, 7.0, 8.0,
        9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0,
        17.0, 18.0, 19.0, 20.0,
    };

    CHECK(mx(0, 0) == 1.0);
    CHECK(mx(0, 1) == 2.0);
    CHECK(mx(0, 2) == 3.0);
    CHECK(mx(0, 3) == 4.0);
    CHECK(mx(1, 0) == 5.0);
    CHECK(mx(1, 1) == 6.0);
    CHECK(mx(1, 2) == 7.0);
    CHECK(mx(1, 3) == 8.0);
    CHECK(mx(2, 0) == 9.0);
    CHECK(mx(2, 1) == 10.0);
    CHECK(mx(2, 2) == 11.0);
    CHECK(mx(2, 3) == 12.0);
    CHECK(mx(3, 0) == 13.0);
    CHECK(mx(3, 1) == 14.0);
    CHECK(mx(3, 2) == 15.0);
    CHECK(mx(3, 3) == 16.0);
    CHECK(mx(4, 0) == 17.0);
    CHECK(mx(4, 1) == 18.0);
    CHECK(mx(4, 2) == 19.0);
    CHECK(mx(4, 3) == 20.0);
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

    CHECK(mx1 == mx1);

    CHECK(mx1 != mx2);

    mx2(1, 0) = 4.0;
    mx2(1, 1) = 5.0;
    mx2(1, 2) = 6.0;

    CHECK(mx1 == mx2);
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

    CHECK(mx3(0, 0) == 11.0);
    CHECK(mx3(0, 1) == 13.0);
    CHECK(mx3(0, 2) == 15.0);
    CHECK(mx3(1, 0) == 24.0);
    CHECK(mx3(1, 1) == 26.0);
    CHECK(mx3(1, 2) == 28.0);
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

    CHECK(mx3(0, 0) == 140.0);
    CHECK(mx3(0, 1) == 146.0);
    CHECK(mx3(1, 0) == 320.0);
    CHECK(mx3(1, 1) == 335.0);
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

    CHECK(part11 == mx11);
    CHECK(part12 == mx12);
    CHECK(part21 == mx21);
    CHECK(part22 == mx22);
}


TEST_CASE("matrix partition addition")
{
    matrix<4, 6> mx1{
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0,
    };

    auto mx1_part11 = mx1.partition<0, 0, 2, 3>();
    auto mx1_part12 = mx1.partition<0, 3, 2, 3>();
    auto mx1_part21 = mx1.partition<2, 0, 2, 3>();
    auto mx1_part22 = mx1.partition<2, 3, 2, 3>();

    matrix<4, 6> mx2{
        -1.0, -2.0, -3.0, -4.0, -5.0, -6.0,
        -7.0, -8.0, -9.0, -10.0, -11.0, -12.0,
        -13.0, -14.0, -15.0, -16.0, -17.0, -18.0,
        -19.0, -20.0, -21.0, -22.0, -23.0, -24.0,
    };

    auto mx2_part11 = mx2.partition<0, 0, 2, 3>();
    auto mx2_part12 = mx2.partition<0, 3, 2, 3>();
    auto mx2_part21 = mx2.partition<2, 0, 2, 3>();
    auto mx2_part22 = mx2.partition<2, 3, 2, 3>();

    matrix<4, 6> mx3;

    auto mx3_part11 = mx3.partition<0, 0, 2, 3>();
    auto mx3_part12 = mx3.partition<0, 3, 2, 3>();
    auto mx3_part21 = mx3.partition<2, 0, 2, 3>();
    auto mx3_part22 = mx3.partition<2, 3, 2, 3>();

    add(mx1, mx2, mx3);

    matrix<4, 6> mx4;

    auto mx4_part11 = mx4.partition<0, 0, 2, 3>();
    auto mx4_part12 = mx4.partition<0, 3, 2, 3>();
    auto mx4_part21 = mx4.partition<2, 0, 2, 3>();
    auto mx4_part22 = mx4.partition<2, 3, 2, 3>();

    add(mx1_part11, mx2_part11, mx4_part11);
    add(mx1_part12, mx2_part12, mx4_part12);
    add(mx1_part21, mx2_part21, mx4_part21);
    add(mx1_part22, mx2_part22, mx4_part22);

    CHECK(mx3_part11 == mx4_part11);
    CHECK(mx3_part12 == mx4_part12);
    CHECK(mx3_part21 == mx4_part21);
    CHECK(mx3_part22 == mx4_part22);

    CHECK(mx3 == mx4);
}


TEST_CASE("matrix partition multiplication")
{
    matrix<4, 6> mx1{
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0,
    };

    auto mx1_part11 = mx1.partition<0, 0, 2, 3>();
    auto mx1_part12 = mx1.partition<0, 3, 2, 3>();
    auto mx1_part21 = mx1.partition<2, 0, 2, 3>();
    auto mx1_part22 = mx1.partition<2, 3, 2, 3>();

    matrix<6, 4> mx2{
        -1.0, -2.0, -3.0, -4.0,
        -5.0, -6.0, -7.0, -8.0,
        -9.0, -10.0, -11.0, -12.0,
        -13.0, -14.0, -15.0, -16.0,
        -17.0, -18.0, -19.0, -20.0,
        -21.0, -22.0, -23.0, -24.0,
    };

    auto mx2_part11 = mx2.partition<0, 0, 3, 2>();
    auto mx2_part12 = mx2.partition<0, 2, 3, 2>();
    auto mx2_part21 = mx2.partition<3, 0, 3, 2>();
    auto mx2_part22 = mx2.partition<3, 2, 3, 2>();

    matrix<4, 4> mx3;

    auto mx3_part11 = mx3.partition<0, 0, 2, 2>();
    auto mx3_part12 = mx3.partition<0, 2, 2, 2>();
    auto mx3_part21 = mx3.partition<2, 0, 2, 2>();
    auto mx3_part22 = mx3.partition<2, 2, 2, 2>();

    mul(mx1, mx2, mx3);

    matrix<4, 4> mx4;

    auto mx4_part11 = mx4.partition<0, 0, 2, 2>();
    auto mx4_part12 = mx4.partition<0, 2, 2, 2>();
    auto mx4_part21 = mx4.partition<2, 0, 2, 2>();
    auto mx4_part22 = mx4.partition<2, 2, 2, 2>();

    matrix<2, 2> mx_tmp1;
    matrix<2, 2> mx_tmp2;

    mul(mx1_part11, mx2_part11, mx_tmp1);
    mul(mx1_part12, mx2_part21, mx_tmp2);
    add(mx_tmp1, mx_tmp2, mx4_part11);

    mul(mx1_part11, mx2_part12, mx_tmp1);
    mul(mx1_part12, mx2_part22, mx_tmp2);
    add(mx_tmp1, mx_tmp2, mx4_part12);

    mul(mx1_part21, mx2_part11, mx_tmp1);
    mul(mx1_part22, mx2_part21, mx_tmp2);
    add(mx_tmp1, mx_tmp2, mx4_part21);

    mul(mx1_part21, mx2_part12, mx_tmp1);
    mul(mx1_part22, mx2_part22, mx_tmp2);
    add(mx_tmp1, mx_tmp2, mx4_part22);

    CHECK(mx3_part11 == mx4_part11);
    CHECK(mx3_part12 == mx4_part12);
    CHECK(mx3_part21 == mx4_part21);
    CHECK(mx3_part22 == mx4_part22);

    CHECK(mx3 == mx4);
}
