#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch.hpp>

#include <iomanip>
#include <iostream>
#include <random>

#include <fmmi/fmmi.hpp>

using namespace fmmi;

static const matrix<double, 7, 7> mx1{
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
    8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0,
    15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0,
    22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0,
    29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0,
    36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0,
    43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0,
};

static const matrix<double, 7, 7> mx2{
    -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0,
    -8.0, -9.0, -10.0, -11.0, -12.0, -13.0, -14.0,
    -15.0, -16.0, -17.0, -18.0, -19.0, -20.0, -21.0,
    -22.0, -23.0, -24.0, -25.0, -26.0, -27.0, -28.0,
    -29.0, -30.0, -31.0, -32.0, -33.0, -34.0, -35.0,
    -36.0, -37.0, -38.0, -39.0, -40.0, -41.0, -42.0,
    -43.0, -44.0, -45.0, -46.0, -47.0, -48.0, -49.0,
};

static matrix<double, 7, 7> mx3;
static matrix<double, 7, 7> mx4;


TEMPLATE_TEST_CASE_SIG("fmmi mul_fast up to 7x7", "[fmmi][mul_fast][equals]",
        ((uint16_t m, uint16_t n, uint16_t p), m, n, p),
        (1, 1, 1), (1, 1, 2), (1, 1, 3), (1, 1, 4), (1, 1, 5), (1, 1, 6), (1, 1, 7),
        (1, 2, 1), (1, 2, 2), (1, 2, 3), (1, 2, 4), (1, 2, 5), (1, 2, 6), (1, 2, 7),
        (1, 3, 1), (1, 3, 2), (1, 3, 3), (1, 3, 4), (1, 3, 5), (1, 3, 6), (1, 3, 7),
        (1, 4, 1), (1, 4, 2), (1, 4, 3), (1, 4, 4), (1, 4, 5), (1, 4, 6), (1, 4, 7),
        (1, 5, 1), (1, 5, 2), (1, 5, 3), (1, 5, 4), (1, 5, 5), (1, 5, 6), (1, 5, 7),
        (1, 6, 1), (1, 6, 2), (1, 6, 3), (1, 6, 4), (1, 6, 5), (1, 6, 6), (1, 6, 7),
        (1, 7, 1), (1, 7, 2), (1, 7, 3), (1, 7, 4), (1, 7, 5), (1, 7, 6), (1, 7, 7),
        (2, 1, 1), (2, 1, 2), (2, 1, 3), (2, 1, 4), (2, 1, 5), (2, 1, 6), (2, 1, 7),
        (2, 2, 1), (2, 2, 2), (2, 2, 3), (2, 2, 4), (2, 2, 5), (2, 2, 6), (2, 2, 7),
        (2, 3, 1), (2, 3, 2), (2, 3, 3), (2, 3, 4), (2, 3, 5), (2, 3, 6), (2, 3, 7),
        (2, 4, 1), (2, 4, 2), (2, 4, 3), (2, 4, 4), (2, 4, 5), (2, 4, 6), (2, 4, 7),
        (2, 5, 1), (2, 5, 2), (2, 5, 3), (2, 5, 4), (2, 5, 5), (2, 5, 6), (2, 5, 7),
        (2, 6, 1), (2, 6, 2), (2, 6, 3), (2, 6, 4), (2, 6, 5), (2, 6, 6), (2, 6, 7),
        (2, 7, 1), (2, 7, 2), (2, 7, 3), (2, 7, 4), (2, 7, 5), (2, 7, 6), (2, 7, 7),
        (3, 1, 1), (3, 1, 2), (3, 1, 3), (3, 1, 4), (3, 1, 5), (3, 1, 6), (3, 1, 7),
        (3, 2, 1), (3, 2, 2), (3, 2, 3), (3, 2, 4), (3, 2, 5), (3, 2, 6), (3, 2, 7),
        (3, 3, 1), (3, 3, 2), (3, 3, 3), (3, 3, 4), (3, 3, 5), (3, 3, 6), (3, 3, 7),
        (3, 4, 1), (3, 4, 2), (3, 4, 3), (3, 4, 4), (3, 4, 5), (3, 4, 6), (3, 4, 7),
        (3, 5, 1), (3, 5, 2), (3, 5, 3), (3, 5, 4), (3, 5, 5), (3, 5, 6), (3, 5, 7),
        (3, 6, 1), (3, 6, 2), (3, 6, 3), (3, 6, 4), (3, 6, 5), (3, 6, 6), (3, 6, 7),
        (3, 7, 1), (3, 7, 2), (3, 7, 3), (3, 7, 4), (3, 7, 5), (3, 7, 6), (3, 7, 7),
        (4, 1, 1), (4, 1, 2), (4, 1, 3), (4, 1, 4), (4, 1, 5), (4, 1, 6), (4, 1, 7),
        (4, 2, 1), (4, 2, 2), (4, 2, 3), (4, 2, 4), (4, 2, 5), (4, 2, 6), (4, 2, 7),
        (4, 3, 1), (4, 3, 2), (4, 3, 3), (4, 3, 4), (4, 3, 5), (4, 3, 6), (4, 3, 7),
        (4, 4, 1), (4, 4, 2), (4, 4, 3), (4, 4, 4), (4, 4, 5), (4, 4, 6), (4, 4, 7),
        (4, 5, 1), (4, 5, 2), (4, 5, 3), (4, 5, 4), (4, 5, 5), (4, 5, 6), (4, 5, 7),
        (4, 6, 1), (4, 6, 2), (4, 6, 3), (4, 6, 4), (4, 6, 5), (4, 6, 6), (4, 6, 7),
        (4, 7, 1), (4, 7, 2), (4, 7, 3), (4, 7, 4), (4, 7, 5), (4, 7, 6), (4, 7, 7),
        (5, 1, 1), (5, 1, 2), (5, 1, 3), (5, 1, 4), (5, 1, 5), (5, 1, 6), (5, 1, 7),
        (5, 2, 1), (5, 2, 2), (5, 2, 3), (5, 2, 4), (5, 2, 5), (5, 2, 6), (5, 2, 7),
        (5, 3, 1), (5, 3, 2), (5, 3, 3), (5, 3, 4), (5, 3, 5), (5, 3, 6), (5, 3, 7),
        (5, 4, 1), (5, 4, 2), (5, 4, 3), (5, 4, 4), (5, 4, 5), (5, 4, 6), (5, 4, 7),
        (5, 5, 1), (5, 5, 2), (5, 5, 3), (5, 5, 4), (5, 5, 5), (5, 5, 6), (5, 5, 7),
        (5, 6, 1), (5, 6, 2), (5, 6, 3), (5, 6, 4), (5, 6, 5), (5, 6, 6), (5, 6, 7),
        (5, 7, 1), (5, 7, 2), (5, 7, 3), (5, 7, 4), (5, 7, 5), (5, 7, 6), (5, 7, 7),
        (6, 1, 1), (6, 1, 2), (6, 1, 3), (6, 1, 4), (6, 1, 5), (6, 1, 6), (6, 1, 7),
        (6, 2, 1), (6, 2, 2), (6, 2, 3), (6, 2, 4), (6, 2, 5), (6, 2, 6), (6, 2, 7),
        (6, 3, 1), (6, 3, 2), (6, 3, 3), (6, 3, 4), (6, 3, 5), (6, 3, 6), (6, 3, 7),
        (6, 4, 1), (6, 4, 2), (6, 4, 3), (6, 4, 4), (6, 4, 5), (6, 4, 6), (6, 4, 7),
        (6, 5, 1), (6, 5, 2), (6, 5, 3), (6, 5, 4), (6, 5, 5), (6, 5, 6), (6, 5, 7),
        (6, 6, 1), (6, 6, 2), (6, 6, 3), (6, 6, 4), (6, 6, 5), (6, 6, 6), (6, 6, 7),
        (6, 7, 1), (6, 7, 2), (6, 7, 3), (6, 7, 4), (6, 7, 5), (6, 7, 6), (6, 7, 7),
        (7, 1, 1), (7, 1, 2), (7, 1, 3), (7, 1, 4), (7, 1, 5), (7, 1, 6), (7, 1, 7),
        (7, 2, 1), (7, 2, 2), (7, 2, 3), (7, 2, 4), (7, 2, 5), (7, 2, 6), (7, 2, 7),
        (7, 3, 1), (7, 3, 2), (7, 3, 3), (7, 3, 4), (7, 3, 5), (7, 3, 6), (7, 3, 7),
        (7, 4, 1), (7, 4, 2), (7, 4, 3), (7, 4, 4), (7, 4, 5), (7, 4, 6), (7, 4, 7),
        (7, 5, 1), (7, 5, 2), (7, 5, 3), (7, 5, 4), (7, 5, 5), (7, 5, 6), (7, 5, 7),
        (7, 6, 1), (7, 6, 2), (7, 6, 3), (7, 6, 4), (7, 6, 5), (7, 6, 6), (7, 6, 7),
        (7, 7, 1), (7, 7, 2), (7, 7, 3), (7, 7, 4), (7, 7, 5), (7, 7, 6), (7, 7, 7))
{
    const auto& a = mx1.partition<m, n>();
    const auto& b = mx2.partition<n, p>();
    auto& c = mx3.partition<m, p>();
    auto& d = mx4.partition<m, p>();

    mul(a, b, c);

    mul_fast(a, b, d);

    CHECK(c == d);
}


static constexpr uint16_t S = 4096;
static i16mx_t<S, S> i16mx_1, i16mx_2, i16mx_3;
static i32mx_t<S, S> i32mx_1, i32mx_2, i32mx_3;
static i64mx_t<S, S> i64mx_1, i64mx_2, i64mx_3;
static f32mx_t<S, S> f32mx_1, f32mx_2, f32mx_3;
static f64mx_t<S, S> f64mx_1, f64mx_2, f64mx_3;

struct static_init
{
    static_init()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(-200, +200);

        for (uint16_t i = 0; i < S; ++i)
        {
            for (uint16_t j = 0; j < S; ++j)
            {
                i16mx_1(i, j) = distrib(gen);
                i16mx_2(i, j) = distrib(gen);
                i32mx_1(i, j) = distrib(gen);
                i32mx_2(i, j) = distrib(gen);
                i64mx_1(i, j) = distrib(gen);
                i64mx_2(i, j) = distrib(gen);
                f32mx_1(i, j) = distrib(gen);
                f32mx_2(i, j) = distrib(gen);
                f64mx_1(i, j) = distrib(gen);
                f64mx_2(i, j) = distrib(gen);
            }
        }
    }
} static_init;


TEMPLATE_TEST_CASE_SIG("fmmi mul_fast benchmark", "[fmmi][mul_fast][benchmark]",
                       ((uint16_t m, uint16_t n, uint16_t p), m, n, p),
                       (1, 1, 1),
                       (2, 2, 2),
                       (3, 3, 3),
                       (4, 4, 4),
                       (7, 7, 7),
                       (8, 8, 8),
                       (15, 15, 15),
                       (16, 16, 16),
                       (31, 31, 31),
                       (32, 32, 32),
                       (63, 63, 63),
                       (64, 64, 64)
//                       (127, 127, 127),
//                       (128, 128, 128)
//                       (255, 255, 255)
//                       (256, 256, 256)
//                       (512, 512, 512),
//                       (1024, 1024, 1024)
                       )
{
    const auto& i16mx_a = i16mx_1.partition<m, n>();
    const auto& i16mx_b = i16mx_2.partition<n, p>();
    auto& i16mx_c = i16mx_3.partition<m, p>();

    BENCHMARK("i16 mul")
    {
        mul(i16mx_a, i16mx_b, i16mx_c);
    };

    BENCHMARK("i16 mul_fast")
    {
        mul_fast(i16mx_a, i16mx_b, i16mx_c);
    };

    const auto& i32mx_a = i32mx_1.partition<m, n>();
    const auto& i32mx_b = i32mx_2.partition<n, p>();
    auto& i32mx_c = i32mx_3.partition<m, p>();

    BENCHMARK("i32 mul")
    {
        mul(i32mx_a, i32mx_b, i32mx_c);
    };

    BENCHMARK("i32 mul_fast")
    {
        mul_fast(i32mx_a, i32mx_b, i32mx_c);
    };

    const auto& i64mx_a = i64mx_1.partition<m, n>();
    const auto& i64mx_b = i64mx_2.partition<n, p>();
    auto& i64mx_c = i64mx_3.partition<m, p>();

    BENCHMARK("i64 mul")
    {
        mul(i64mx_a, i64mx_b, i64mx_c);
    };

    BENCHMARK("i64 mul_fast")
    {
        mul_fast(i64mx_a, i64mx_b, i64mx_c);
    };

    const auto& f32mx_a = f32mx_1.partition<m, n>();
    const auto& f32mx_b = f32mx_2.partition<n, p>();
    auto& f32mx_c = f32mx_3.partition<m, p>();

    BENCHMARK("f32 mul")
    {
        mul(f32mx_a, f32mx_b, f32mx_c);
    };

    BENCHMARK("f32 mul_fast")
    {
        mul_fast(f32mx_a, f32mx_b, f32mx_c);
    };

    const auto& f64mx_a = f64mx_1.partition<m, n>();
    const auto& f64mx_b = f64mx_2.partition<n, p>();
    auto& f64mx_c = f64mx_3.partition<m, p>();

    BENCHMARK("f64 mul")
    {
        mul(f64mx_a, f64mx_b, f64mx_c);
    };

    BENCHMARK("f64 mul_fast")
    {
        mul_fast(f64mx_a, f64mx_b, f64mx_c);
    };
}


TEST_CASE("inv_fast", "[inverse]")
{
    f32mx_t<1, 1> a{
        5,
    };

    f32mx_t<1, 1> expected_a_inv{
        0.2,
    };

    f32mx_t<1, 1> identity1x1 = f32mx_t<1, 1>::identity();
    f32mx_t<1, 1> a_expected_a_inv;
    mul(a, expected_a_inv, a_expected_a_inv);
    CHECK(a_expected_a_inv.equals(identity1x1, 0.001));

    f32mx_t<1, 1> a_inv;
    inv_fast(a, a_inv);
    CHECK(a_inv.equals(expected_a_inv));

    f32mx_t<1, 1> a_a_inv;
    mul(a, a_inv, a_a_inv);
    CHECK(a_a_inv.equals(identity1x1));

    f32mx_t<2, 2> b{
        5, 6,
        2, 2,
    };

    f32mx_t<2, 2> expected_b_inv{
        -1, 3,
        1, -2.5,
    };

    f32mx_t<2, 2> identity2x2 = f32mx_t<2, 2>::identity();
    f32mx_t<2, 2> b_expected_b_inv;
    mul(b, expected_b_inv, b_expected_b_inv);
    CHECK(b_expected_b_inv.equals(identity2x2, 1e-7));

    f32mx_t<2, 2> b_inv;
    inv_fast(b, b_inv);
    CHECK(b_inv.equals(expected_b_inv, 1e-7));

    f32mx_t<2, 2> b_b_inv;
    mul(b, b_inv, b_b_inv);
    CHECK(b_b_inv.equals(identity2x2, 1e-7));

    f32mx_t<4, 4> c{
        5, 6, 6, 8,
        2, 2, 2, 8,
        6, 6, 2, 8,
        2, 3, 6, 7,
    };

    f32mx_t<4, 4> expected_c_inv{
        -17, -9, 12, 16,
        17, 8.75, -11.75, -16,
        -4, -2.25, 2.75, 4,
        1, 0.75, -0.75, -1,
    };

    f32mx_t<4, 4> identity4x4 = f32mx_t<4, 4>::identity();
    f32mx_t<4, 4> c_expected_c_inv;
    mul(c, expected_c_inv, c_expected_c_inv);
    CHECK(c_expected_c_inv.equals(identity4x4, 1e-7));

    f32mx_t<4, 4> c_inv;
    inv_fast(c, c_inv);
    CHECK(c_inv.equals(expected_c_inv, 1e-7));

    f32mx_t<4, 4> c_c_inv;
    mul(c, c_inv, c_c_inv);
    CHECK(c_c_inv.equals(identity4x4, 1e-7));

    f32mx_t<3, 3> d{
        0, -3, -2,
        1, -4, -2,
        -3, 4, 1,
    };

    f32mx_t<3, 3> expected_d_inv{
        4, -5, -2,
        5, -6, -2,
        -8, 9, 3,
    };

    f32mx_t<3, 3> identity3x3 = f32mx_t<3, 3>::identity();
    f32mx_t<3, 3> d_expected_d_inv;
    mul(d, expected_d_inv, d_expected_d_inv);
    CHECK(d_expected_d_inv.equals(identity3x3, 1e-7));

    f32mx_t<3, 3> d_inv;
    inv_fast(d, d_inv);
    CHECK(d_inv.equals(expected_d_inv, 1e-7));

    f32mx_t<3, 3> d_d_inv;
    mul(d, d_inv, d_d_inv);
    CHECK(d_d_inv.equals(identity3x3, 1e-7));
}
