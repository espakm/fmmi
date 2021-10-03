#include <catch.hpp>
#include <fmmi/dmatrix.hpp>

using namespace fmmi;

TEST_CASE("dmatrix ctor", "[dmatrix][ctor]")
{
    dmatrix<double> mx(5, 4, {
        1.0, 2.0, 3.0, 4.0,
        5.0, 6.0, 7.0, 8.0,
        9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0,
        17.0, 18.0, 19.0, 20.0,
    });

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


TEST_CASE("dmatrix ==", "[dmatrix][==]")
{
    dmatrix<double> mx1(2, 3, {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
    });

    dmatrix<double> mx2(2, 3, {
        1.0, 2.0, 3.0,
        5.0, 6.0, 7.0,
    });

    CHECK(mx1 == mx1);

    CHECK(mx1 != mx2);

    mx2(1, 0) = 4.0;
    mx2(1, 1) = 5.0;
    mx2(1, 2) = 6.0;

    CHECK(mx1 == mx2);
}


TEST_CASE("add dmatrix", "[add][dmatrix]")
{
    dmatrix<double> mx1(2, 3, {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
    });

    dmatrix<double> mx2(2, 3, {
        10.0, 11.0, 12.0,
        20.0, 21.0, 22.0,
    });

    dmatrix<double> mx3(2, 3);

    add(mx1, mx2, mx3);

    CHECK(mx3(0, 0) == 11.0);
    CHECK(mx3(0, 1) == 13.0);
    CHECK(mx3(0, 2) == 15.0);
    CHECK(mx3(1, 0) == 24.0);
    CHECK(mx3(1, 1) == 26.0);
    CHECK(mx3(1, 2) == 28.0);
}


TEST_CASE("mul dmatrix", "[mul][dmatrix]")
{
    dmatrix<double> mx1(2, 3, {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
    });

    dmatrix<double> mx2(3, 2, {
        10.0, 11.0,
        20.0, 21.0,
        30.0, 31.0,
    });

    dmatrix<double> mx3(2, 2);

    mul(mx1, mx2, mx3);

    CHECK(mx3(0, 0) == 140.0);
    CHECK(mx3(0, 1) == 146.0);
    CHECK(mx3(1, 0) == 320.0);
    CHECK(mx3(1, 1) == 335.0);
}


TEST_CASE("dmatrix partition", "[dmatrix][partition]")
{
    dmatrix<double> mx(5, 4, {
        1.0, 2.0, 3.0, 4.0,
        5.0, 6.0, 7.0, 8.0,
        9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0,
        17.0, 18.0, 19.0, 20.0,
    });

    dmatrix<double> mx11(2, 3, {
        1.0, 2.0, 3.0,
        5.0, 6.0, 7.0,
    });

    dmatrix<double> mx12(2, 1, {
        4.0,
        8.0,
    });

    dmatrix<double> mx21(3, 1, {
        9.0,
        13.0,
        17.0,
    });

    dmatrix<double> mx22(3, 3, {
        10.0, 11.0, 12.0,
        14.0, 15.0, 16.0,
        18.0, 19.0, 20.0,
    });

    const auto& part11 = mx.partition(2, 3, 0, 0);
    const auto& part12 = mx.partition(2, 1, 0, 3);
    const auto& part21 = mx.partition(3, 1, 2, 0);
    const auto& part22 = mx.partition(3, 3, 2, 1);

    CHECK(part11 == mx11);
    CHECK(part12 == mx12);
    CHECK(part21 == mx21);
    CHECK(part22 == mx22);
}


TEST_CASE("add dmatrix partition", "[add][dmatrix][partition]")
{
    dmatrix<double> mx1(4, 6, {
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0,
    });

    const auto& mx1_part11 = mx1.partition(2, 3, 0, 0);
    const auto& mx1_part12 = mx1.partition(2, 3, 0, 3);
    const auto& mx1_part21 = mx1.partition(2, 3, 2, 0);
    const auto& mx1_part22 = mx1.partition(2, 3, 2, 3);

    dmatrix<double> mx2(4, 6, {
        -1.0, -2.0, -3.0, -4.0, -5.0, -6.0,
        -7.0, -8.0, -9.0, -10.0, -11.0, -12.0,
        -13.0, -14.0, -15.0, -16.0, -17.0, -18.0,
        -19.0, -20.0, -21.0, -22.0, -23.0, -24.0,
    });

    const auto& mx2_part11 = mx2.partition(2, 3, 0, 0);
    const auto& mx2_part12 = mx2.partition(2, 3, 0, 3);
    const auto& mx2_part21 = mx2.partition(2, 3, 2, 0);
    const auto& mx2_part22 = mx2.partition(2, 3, 2, 3);

    dmatrix<double> mx3(4, 6);

    auto mx3_part11 = mx3.partition(2, 3, 0, 0);
    auto mx3_part12 = mx3.partition(2, 3, 0, 3);
    auto mx3_part21 = mx3.partition(2, 3, 2, 0);
    auto mx3_part22 = mx3.partition(2, 3, 2, 3);

    add(mx1, mx2, mx3);

    dmatrix<double> mx4(4, 6);

    auto mx4_part11 = mx4.partition(2, 3, 0, 0);
    auto mx4_part12 = mx4.partition(2, 3, 0, 3);
    auto mx4_part21 = mx4.partition(2, 3, 2, 0);
    auto mx4_part22 = mx4.partition(2, 3, 2, 3);

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


TEST_CASE("mul dmatrix partition", "[mul][dmatrix][partition]")
{
    dmatrix<double> mx1(4, 6, {
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0,
    });

    const auto& mx1_part11 = mx1.partition(2, 3, 0, 0);
    const auto& mx1_part12 = mx1.partition(2, 3, 0, 3);
    const auto& mx1_part21 = mx1.partition(2, 3, 2, 0);
    const auto& mx1_part22 = mx1.partition(2, 3, 2, 3);

    dmatrix<double> mx2(6, 4, {
        -1.0, -2.0, -3.0, -4.0,
        -5.0, -6.0, -7.0, -8.0,
        -9.0, -10.0, -11.0, -12.0,
        -13.0, -14.0, -15.0, -16.0,
        -17.0, -18.0, -19.0, -20.0,
        -21.0, -22.0, -23.0, -24.0,
    });

    const auto& mx2_part11 = mx2.partition(3, 2, 0, 0);
    const auto& mx2_part12 = mx2.partition(3, 2, 0, 2);
    const auto& mx2_part21 = mx2.partition(3, 2, 3, 0);
    const auto& mx2_part22 = mx2.partition(3, 2, 3, 2);

    dmatrix<double> mx3(4, 4);

    auto mx3_part11 = mx3.partition(2, 2, 0, 0);
    auto mx3_part12 = mx3.partition(2, 2, 0, 2);
    auto mx3_part21 = mx3.partition(2, 2, 2, 0);
    auto mx3_part22 = mx3.partition(2, 2, 2, 2);

    mul(mx1, mx2, mx3);

    dmatrix<double> mx4(4, 4);

    auto mx4_part11 = mx4.partition(2, 2, 0, 0);
    auto mx4_part12 = mx4.partition(2, 2, 0, 2);
    auto mx4_part21 = mx4.partition(2, 2, 2, 0);
    auto mx4_part22 = mx4.partition(2, 2, 2, 2);

    dmatrix<double> mx_tmp1(2, 2);
    dmatrix<double> mx_tmp2(2, 2);

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


TEST_CASE("inv dmatrix", "[inv][dmatrix]")
{
    f32dmx_t a(1, 1, {
        5,
    });

    f32dmx_t expected_a_inv(1, 1, {
        0.2,
    });

    f32dmx_t identity1x1 = f32dmx_t::identity(1, 1);
    f32dmx_t a_expected_a_inv(1, 1);
    mul(a, expected_a_inv, a_expected_a_inv);
    CHECK(a_expected_a_inv.equals(identity1x1, 0.001));

    f32dmx_t a_inv(1, 1);
    inv(a, a_inv);
    CHECK(a_inv.equals(expected_a_inv));

    f32dmx_t a_a_inv(1, 1);
    mul(a, a_inv, a_a_inv);
    CHECK(a_a_inv.equals(identity1x1));

    f32dmx_t b(2, 2, {
        5, 6,
        2, 2,
    });

    f32dmx_t expected_b_inv(2, 2, {
        -1, 3,
        1, -2.5,
    });

    f32dmx_t identity2x2 = f32dmx_t::identity(2, 2);
    f32dmx_t b_expected_b_inv(2, 2);
    mul(b, expected_b_inv, b_expected_b_inv);
    CHECK(b_expected_b_inv.equals(identity2x2, 1e-7));

    f32dmx_t b_inv(2, 2);
    inv(b, b_inv);
    CHECK(b_inv.equals(expected_b_inv, 1e-7));

    f32dmx_t b_b_inv(2, 2);
    mul(b, b_inv, b_b_inv);
    CHECK(b_b_inv.equals(identity2x2, 1e-7));

    f32dmx_t c(4, 4, {
        5, 6, 6, 8,
        2, 2, 2, 8,
        6, 6, 2, 8,
        2, 3, 6, 7,
    });

    f32dmx_t expected_c_inv(4, 4, {
        -17, -9, 12, 16,
        17, 8.75, -11.75, -16,
        -4, -2.25, 2.75, 4,
        1, 0.75, -0.75, -1,
    });

    f32dmx_t identity4x4 = f32dmx_t::identity(4, 4);
    f32dmx_t c_expected_c_inv(4, 4);
    mul(c, expected_c_inv, c_expected_c_inv);
    CHECK(c_expected_c_inv.equals(identity4x4, 1e-7));

    f32dmx_t c_inv(4, 4);
    inv(c, c_inv);
    CHECK(c_inv.equals(expected_c_inv, 1e-7));

    f32dmx_t c_c_inv(4, 4);
    mul(c, c_inv, c_c_inv);
    CHECK(c_c_inv.equals(identity4x4, 1e-7));

    f32dmx_t d(3, 3, {
        0, -3, -2,
        1, -4, -2,
        -3, 4, 1,
    });

    f32dmx_t expected_d_inv(3, 3, {
        4, -5, -2,
        5, -6, -2,
        -8, 9, 3,
    });

    f32dmx_t identity3x3 = f32dmx_t::identity(3, 3);
    f32dmx_t d_expected_d_inv(3, 3);
    mul(d, expected_d_inv, d_expected_d_inv);
    CHECK(d_expected_d_inv.equals(identity3x3, 1e-7));

    f32dmx_t d_inv(3, 3);
    inv(d, d_inv);
    CHECK(d_inv.equals(expected_d_inv, 1e-7));

    f32dmx_t d_d_inv(3, 3);
    mul(d, d_inv, d_d_inv);
    CHECK(d_d_inv.equals(identity3x3, 1e-7));
}
