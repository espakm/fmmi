#include <catch.hpp>
#include <fmmi/smatrix.hpp>

using namespace fmmi;

TEST_CASE("smatrix ctor", "[smatrix][ctor]")
{
    smatrix<double, 5, 4> mx{
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


TEST_CASE("smatrix ==", "[smatrix][operator==]")
{
    smatrix<double, 2, 3> mx1{
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
    };

    smatrix<double, 2, 3> mx2{
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


TEST_CASE("add smatrix", "[add][smatrix]")
{
    smatrix<double, 2, 3> mx1{
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
    };

    smatrix<double, 2, 3> mx2{
        10.0, 11.0, 12.0,
        20.0, 21.0, 22.0,
    };

    smatrix<double, 2, 3> mx3;

    add(mx1, mx2, mx3);

    CHECK(mx3(0, 0) == 11.0);
    CHECK(mx3(0, 1) == 13.0);
    CHECK(mx3(0, 2) == 15.0);
    CHECK(mx3(1, 0) == 24.0);
    CHECK(mx3(1, 1) == 26.0);
    CHECK(mx3(1, 2) == 28.0);
}


TEST_CASE("mul smatrix", "[mul][smatrix]")
{
    smatrix<double, 2, 3> mx1{
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
    };

    smatrix<double, 3, 2> mx2{
        10.0, 11.0,
        20.0, 21.0,
        30.0, 31.0,
    };

    smatrix<double, 2, 2> mx3;

    mul(mx1, mx2, mx3);

    CHECK(mx3(0, 0) == 140.0);
    CHECK(mx3(0, 1) == 146.0);
    CHECK(mx3(1, 0) == 320.0);
    CHECK(mx3(1, 1) == 335.0);
}


TEST_CASE("smatrix partition", "[smatrix][partition]")
{
    smatrix<double, 5, 4> mx{
        1.0, 2.0, 3.0, 4.0,
        5.0, 6.0, 7.0, 8.0,
        9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0,
        17.0, 18.0, 19.0, 20.0,
    };

    smatrix<double, 2, 3> mx11{
        1.0, 2.0, 3.0,
        5.0, 6.0, 7.0,
    };

    smatrix<double, 2, 1> mx12{
        4.0,
        8.0,
    };

    smatrix<double, 3, 1> mx21{
        9.0,
        13.0,
        17.0,
    };

    smatrix<double, 3, 3> mx22{
        10.0, 11.0, 12.0,
        14.0, 15.0, 16.0,
        18.0, 19.0, 20.0,
    };

    const auto& part11 = mx.partition<2, 3, 0, 0>();
    const auto& part12 = mx.partition<2, 1, 0, 3>();
    const auto& part21 = mx.partition<3, 1, 2, 0>();
    const auto& part22 = mx.partition<3, 3, 2, 1>();

    CHECK(part11 == mx11);
    CHECK(part12 == mx12);
    CHECK(part21 == mx21);
    CHECK(part22 == mx22);
}


TEST_CASE("add smatrix partition", "[add][smatrix][partition]")
{
    smatrix<double, 4, 6> mx1{
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0,
    };

    const auto& mx1_part11 = mx1.partition<2, 3, 0, 0>();
    const auto& mx1_part12 = mx1.partition<2, 3, 0, 3>();
    const auto& mx1_part21 = mx1.partition<2, 3, 2, 0>();
    const auto& mx1_part22 = mx1.partition<2, 3, 2, 3>();

    smatrix<double, 4, 6> mx2{
        -1.0, -2.0, -3.0, -4.0, -5.0, -6.0,
        -7.0, -8.0, -9.0, -10.0, -11.0, -12.0,
        -13.0, -14.0, -15.0, -16.0, -17.0, -18.0,
        -19.0, -20.0, -21.0, -22.0, -23.0, -24.0,
    };

    const auto& mx2_part11 = mx2.partition<2, 3, 0, 0>();
    const auto& mx2_part12 = mx2.partition<2, 3, 0, 3>();
    const auto& mx2_part21 = mx2.partition<2, 3, 2, 0>();
    const auto& mx2_part22 = mx2.partition<2, 3, 2, 3>();

    smatrix<double, 4, 6> mx3;

    auto& mx3_part11 = mx3.partition<2, 3, 0, 0>();
    auto& mx3_part12 = mx3.partition<2, 3, 0, 3>();
    auto& mx3_part21 = mx3.partition<2, 3, 2, 0>();
    auto& mx3_part22 = mx3.partition<2, 3, 2, 3>();

    add(mx1, mx2, mx3);

    smatrix<double, 4, 6> mx4;

    auto& mx4_part11 = mx4.partition<2, 3, 0, 0>();
    auto& mx4_part12 = mx4.partition<2, 3, 0, 3>();
    auto& mx4_part21 = mx4.partition<2, 3, 2, 0>();
    auto& mx4_part22 = mx4.partition<2, 3, 2, 3>();

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


TEST_CASE("mul smatrix partition", "[mul][smatrix][partition]")
{
    smatrix<double, 4, 6> mx1{
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0,
    };

    const auto& mx1_part11 = mx1.partition<2, 3, 0, 0>();
    const auto& mx1_part12 = mx1.partition<2, 3, 0, 3>();
    const auto& mx1_part21 = mx1.partition<2, 3, 2, 0>();
    const auto& mx1_part22 = mx1.partition<2, 3, 2, 3>();

    smatrix<double, 6, 4> mx2{
        -1.0, -2.0, -3.0, -4.0,
        -5.0, -6.0, -7.0, -8.0,
        -9.0, -10.0, -11.0, -12.0,
        -13.0, -14.0, -15.0, -16.0,
        -17.0, -18.0, -19.0, -20.0,
        -21.0, -22.0, -23.0, -24.0,
    };

    const auto& mx2_part11 = mx2.partition<3, 2, 0, 0>();
    const auto& mx2_part12 = mx2.partition<3, 2, 0, 2>();
    const auto& mx2_part21 = mx2.partition<3, 2, 3, 0>();
    const auto& mx2_part22 = mx2.partition<3, 2, 3, 2>();

    smatrix<double, 4, 4> mx3;

    auto& mx3_part11 = mx3.partition<2, 2, 0, 0>();
    auto& mx3_part12 = mx3.partition<2, 2, 0, 2>();
    auto& mx3_part21 = mx3.partition<2, 2, 2, 0>();
    auto& mx3_part22 = mx3.partition<2, 2, 2, 2>();

    mul(mx1, mx2, mx3);

    smatrix<double, 4, 4> mx4;

    auto& mx4_part11 = mx4.partition<2, 2, 0, 0>();
    auto& mx4_part12 = mx4.partition<2, 2, 0, 2>();
    auto& mx4_part21 = mx4.partition<2, 2, 2, 0>();
    auto& mx4_part22 = mx4.partition<2, 2, 2, 2>();

    smatrix<double, 2, 2> mx_tmp1;
    smatrix<double, 2, 2> mx_tmp2;

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


TEST_CASE("inv smatrix", "[inv][smatrix]")
{
    f32smx_t<1, 1> a{
        5,
    };

    f32smx_t<1, 1> expected_a_inv{
        0.2,
    };

    f32smx_t<1, 1> identity1x1 = f32smx_t<1, 1>::identity();
    f32smx_t<1, 1> a_expected_a_inv;
    mul(a, expected_a_inv, a_expected_a_inv);
    CHECK(a_expected_a_inv.equals(identity1x1, 0.001));

    f32smx_t<1, 1> a_inv;
    inv(a, a_inv);
    CHECK(a_inv.equals(expected_a_inv));

    f32smx_t<1, 1> a_a_inv;
    mul(a, a_inv, a_a_inv);
    CHECK(a_a_inv.equals(identity1x1));

    f32smx_t<1, 1> a0{
        0,
    };
    f32smx_t<1, 1> a0_inv;
    CHECK_THROWS_WITH(inv(a0, a0_inv), "Matrix is not invertible.");

    f32smx_t<2, 2> b{
        5, 6,
        2, 2,
    };

    f32smx_t<2, 2> expected_b_inv{
        -1, 3,
        1, -2.5,
    };

    f32smx_t<2, 2> identity2x2 = f32smx_t<2, 2>::identity();
    f32smx_t<2, 2> b_expected_b_inv;
    mul(b, expected_b_inv, b_expected_b_inv);
    CHECK(b_expected_b_inv.equals(identity2x2, 1e-7));

    f32smx_t<2, 2> b_inv;
    inv(b, b_inv);
    CHECK(b_inv.equals(expected_b_inv, 1e-6));

    f32smx_t<2, 2> b_b_inv;
    mul(b, b_inv, b_b_inv);
    CHECK(b_b_inv.equals(identity2x2, 1e-6));

    f32smx_t<2, 2> b0{
        2, 4,
        4, 8,
    };
    f32smx_t<2, 2> b0_inv;
    CHECK_THROWS_WITH(inv(b0, b0_inv), "Matrix is not invertible.");

    f32smx_t<4, 4> c{
        5, 6, 6, 8,
        2, 2, 2, 8,
        6, 6, 2, 8,
        2, 3, 6, 7,
    };

    f32smx_t<4, 4> expected_c_inv{
        -17, -9, 12, 16,
        17, 8.75, -11.75, -16,
        -4, -2.25, 2.75, 4,
        1, 0.75, -0.75, -1,
    };

    f32smx_t<4, 4> identity4x4 = f32smx_t<4, 4>::identity();
    f32smx_t<4, 4> c_expected_c_inv;
    mul(c, expected_c_inv, c_expected_c_inv);
    CHECK(c_expected_c_inv.equals(identity4x4, 1e-6));

    f32smx_t<4, 4> c_inv;
    inv(c, c_inv);
    CHECK(c_inv.equals(expected_c_inv, 1e-4));

    f32smx_t<4, 4> c_c_inv;
    mul(c, c_inv, c_c_inv);
    CHECK(c_c_inv.equals(identity4x4, 1e-4));

    f32smx_t<4, 4> c0{
        5, 6, 6, 8,
        2, 2, 2, 8,
        6, 6, 2, 8,
        -13, -14, -10, -24,
    };
    f32smx_t<4, 4> c0_inv;
    CHECK_THROWS_WITH(inv(c0, c0_inv), "Matrix is not invertible.");

    f32smx_t<3, 3> d{
        0, -3, -2,
        1, -4, -2,
        -3, 4, 1,
    };

    f32smx_t<3, 3> expected_d_inv{
        4, -5, -2,
        5, -6, -2,
        -8, 9, 3,
    };

    f32smx_t<3, 3> identity3x3 = f32smx_t<3, 3>::identity();
    f32smx_t<3, 3> d_expected_d_inv;
    mul(d, expected_d_inv, d_expected_d_inv);
    CHECK(d_expected_d_inv.equals(identity3x3, 1e-7));

    f32smx_t<3, 3> d_inv;
    inv(d, d_inv);
    CHECK(d_inv.equals(expected_d_inv, 1e-5));

    f32smx_t<3, 3> d_d_inv;
    mul(d, d_inv, d_d_inv);
    CHECK(d_d_inv.equals(identity3x3, 1e-6));
}
