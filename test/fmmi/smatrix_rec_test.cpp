#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch.hpp>

#include <random>

#include <fmmi/smatrix_rec.hpp>

#ifndef SMATRIX_REC_TEST_MAX_MATRIX_SIZE
#define SMATRIX_REC_TEST_MAX_MATRIX_SIZE 64
#endif

using namespace fmmi;

static const smatrix<double, 7, 7> smx1{
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
    8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0,
    15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0,
    22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0,
    29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0,
    36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0,
    43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0,
};

static const smatrix<double, 7, 7> smx2{
    -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0,
    -8.0, -9.0, -10.0, -11.0, -12.0, -13.0, -14.0,
    -15.0, -16.0, -17.0, -18.0, -19.0, -20.0, -21.0,
    -22.0, -23.0, -24.0, -25.0, -26.0, -27.0, -28.0,
    -29.0, -30.0, -31.0, -32.0, -33.0, -34.0, -35.0,
    -36.0, -37.0, -38.0, -39.0, -40.0, -41.0, -42.0,
    -43.0, -44.0, -45.0, -46.0, -47.0, -48.0, -49.0,
};

static smatrix<double, 7, 7> smx3;
static smatrix<double, 7, 7> smx4;

TEMPLATE_TEST_CASE_SIG("mul_rec smatrix equals up to 7x7", "[mul_rec][smatrix][equals]",
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
    const auto& sa = smx1.partition<m, n>();
    const auto& sb = smx2.partition<n, p>();
    auto& sc = smx3.partition<m, p>();
    auto& sd = smx4.partition<m, p>();

    mul(sa, sb, sc);

    mul_rec(sa, sb, sd);

    CHECK(sc == sd);
}


static constexpr uint16_t S = 4096;
static i16smx_t<S, S> i16smx_1, i16smx_2, i16smx_3;
static i32smx_t<S, S> i32smx_1, i32smx_2, i32smx_3;
static i64smx_t<S, S> i64smx_1, i64smx_2, i64smx_3;
static f32smx_t<S, S> f32smx_1, f32smx_2, f32smx_3;
static f64smx_t<S, S> f64smx_1, f64smx_2, f64smx_3;
static f64smx_t<S, S> f128smx_1, f128smx_2, f128smx_3;

static
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
                i16smx_1(i, j) = distrib(gen);
                i16smx_2(i, j) = distrib(gen);
                i32smx_1(i, j) = distrib(gen);
                i32smx_2(i, j) = distrib(gen);
                i64smx_1(i, j) = distrib(gen);
                i64smx_2(i, j) = distrib(gen);
                f32smx_1(i, j) = distrib(gen);
                f32smx_2(i, j) = distrib(gen);
                f64smx_1(i, j) = distrib(gen);
                f64smx_2(i, j) = distrib(gen);
                f128smx_1(i, j) = distrib(gen);
                f128smx_2(i, j) = distrib(gen);
            }
        }
    }
} static_init;


TEMPLATE_TEST_CASE_SIG("mul_rec smatrix benchmark", "[mul_rec][smatrix][benchmark]",
                       ((uint16_t m, uint16_t n, uint16_t p), m, n, p)
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 2
                      , (1, 1, 1)
                      , (2, 2, 2)
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 4
                      , (3, 3, 3)
                      , (4, 4, 4)
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 8
                      , (7, 7, 7)
                      , (8, 8, 8)
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 16
                      , (15, 15, 15)
                      , (16, 16, 16)
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 32
                      , (31, 31, 31)
                      , (32, 32, 32)
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 64
                      , (63, 63, 63)
                      , (64, 64, 64)
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 128
                      , (127, 127, 127)
                      , (128, 128, 128)
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 256
                      , (255, 255, 255)
                      , (256, 256, 256)
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 512
                      , (511, 511, 511)
                      , (512, 512, 512)
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 1024
                      , (1023, 1023, 1023)
                      , (1024, 1024, 1024)
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 2048
                      , (2047, 2047, 2047)
                      , (2048, 2048, 2048)
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 4096
                      , (4095, 4095, 4095)
                      , (4096, 4096, 4096)
#endif
                       )
{
    const auto& i16smx_a = i16smx_1.partition<m, n>();
    const auto& i16smx_b = i16smx_2.partition<n, p>();
    auto& i16smx_c = i16smx_3.partition<m, p>();

    BENCHMARK("mul i16smx_t")
    {
        mul(i16smx_a, i16smx_b, i16smx_c);
    };

    BENCHMARK("mul_rec i16smx_t")
    {
        mul_rec(i16smx_a, i16smx_b, i16smx_c);
    };

    const auto& i32smx_a = i32smx_1.partition<m, n>();
    const auto& i32smx_b = i32smx_2.partition<n, p>();
    auto& i32smx_c = i32smx_3.partition<m, p>();

    BENCHMARK("mul i32smx_t")
    {
        mul(i32smx_a, i32smx_b, i32smx_c);
    };

    BENCHMARK("mul_rec i32smx_t")
    {
        mul_rec(i32smx_a, i32smx_b, i32smx_c);
    };

    const auto& i64smx_a = i64smx_1.partition<m, n>();
    const auto& i64smx_b = i64smx_2.partition<n, p>();
    auto& i64smx_c = i64smx_3.partition<m, p>();

    BENCHMARK("mul i64smx_t")
    {
        mul(i64smx_a, i64smx_b, i64smx_c);
    };

    BENCHMARK("mul_rec i64smx_t")
    {
        mul_rec(i64smx_a, i64smx_b, i64smx_c);
    };

    const auto& f32smx_a = f32smx_1.partition<m, n>();
    const auto& f32smx_b = f32smx_2.partition<n, p>();
    auto& f32smx_c = f32smx_3.partition<m, p>();

    BENCHMARK("mul f32smx_t")
    {
        mul(f32smx_a, f32smx_b, f32smx_c);
    };

    BENCHMARK("mul_rec f32smx_t")
    {
        mul_rec(f32smx_a, f32smx_b, f32smx_c);
    };

    const auto& f64smx_a = f64smx_1.partition<m, n>();
    const auto& f64smx_b = f64smx_2.partition<n, p>();
    auto& f64smx_c = f64smx_3.partition<m, p>();

    BENCHMARK("mul f64smx_t")
    {
        mul(f64smx_a, f64smx_b, f64smx_c);
    };

    BENCHMARK("mul_rec f64smx_t")
    {
        mul_rec(f64smx_a, f64smx_b, f64smx_c);
    };

    const auto& f128smx_a = f128smx_1.partition<m, n>();
    const auto& f128smx_b = f128smx_2.partition<n, p>();
    auto& f128smx_c = f128smx_3.partition<m, p>();

    BENCHMARK("mul f128smx_t")
    {
        mul(f128smx_a, f128smx_b, f128smx_c);
    };

    BENCHMARK("mul_rec f128smx_t")
    {
        mul_rec(f128smx_a, f128smx_b, f128smx_c);
    };
}


TEST_CASE("inv_rec smatrix equals", "[inv_rec][smatrix][equals]")
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
    inv_rec(a, a_inv);
    CHECK(a_inv.equals(expected_a_inv));

    f32smx_t<1, 1> a_a_inv;
    mul(a, a_inv, a_a_inv);
    CHECK(a_a_inv.equals(identity1x1));

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
    inv_rec(b, b_inv);
    CHECK(b_inv.equals(expected_b_inv, 1e-7));

    f32smx_t<2, 2> b_b_inv;
    mul(b, b_inv, b_b_inv);
    CHECK(b_b_inv.equals(identity2x2, 1e-7));

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
    CHECK(c_expected_c_inv.equals(identity4x4, 1e-7));

    f32smx_t<4, 4> c_inv;
    inv_rec(c, c_inv);
    CHECK(c_inv.equals(expected_c_inv, 1e-7));

    f32smx_t<4, 4> c_c_inv;
    mul(c, c_inv, c_c_inv);
    CHECK(c_c_inv.equals(identity4x4, 1e-7));

    f32smx_t<4, 4> d{
        0, -3, -2, 0,
        1, -4, -2, 0,
        -3, 4, 1, 0,
        0, 0, 0, 1,
    };

    f32smx_t<4, 4> expected_d_inv{
        4, -5, -2, 0,
        5, -6, -2, 0,
        -8, 9, 3, 0,
        0, 0, 0, 1,
    };

    f32smx_t<4, 4> d_expected_d_inv;
    mul(d, expected_d_inv, d_expected_d_inv);
    CHECK(d_expected_d_inv.equals(identity4x4, 1e-6));

    f32smx_t<4, 4> d_inv;
    inv_rec(d, d_inv);
    CHECK(d_inv.equals(expected_d_inv, 1e-5));

    f32smx_t<4, 4> d_d_inv;
    mul(d, d_inv, d_d_inv);
    CHECK(d_d_inv.equals(identity4x4, 1e-5));

    f32smx_t<3, 3> e{
        0, -3, -2,
        1, -4, -2,
        -3, 4, 1,
    };

    f32smx_t<3, 3> expected_e_inv{
        4, -5, -2,
        5, -6, -2,
        -8, 9, 3,
    };

    f32smx_t<3, 3> identity3x3 = f32smx_t<3, 3>::identity();
    f32smx_t<3, 3> e_expected_e_inv;
    mul(e, expected_e_inv, e_expected_e_inv);
    CHECK(e_expected_e_inv.equals(identity3x3, 1e-6));

    f32smx_t<3, 3> e_inv;
    inv_rec(e, e_inv);
    CHECK(e_inv.equals(expected_e_inv, 1e-5));

    f32smx_t<3, 3> e_e_inv;
    mul(e, e_inv, e_e_inv);
    CHECK(e_e_inv.equals(identity3x3, 1e-5));
}


TEMPLATE_TEST_CASE_SIG("inv_rec smatrix benchmark", "[inv_rec][smatrix][benchmark]",
                       ((uint16_t m), m)
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 2
                       , 1
                       , 2
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 4
                       , 3
                       , 4
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 8
                       , 7
                       , 8
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 16
                       , 15
                       , 16
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 32
                       , 31
                       , 32
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 64
                       , 63
                       , 64
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 128
                       , 127
                       , 128
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 256
                       , 255
                       , 256
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 512
                       , 511
                       , 512
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 1024
                       , 1023
                       , 1024
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 2048
                       , 2047
                       , 2048
#endif
#if SMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 4096
                       , 4095
                       , 4096
#endif
                       )
{
    const auto& i16smx_a = i16smx_1.partition<m, m>();
    auto& i16smx_b = i16smx_2.partition<m, m>();

    BENCHMARK("inv i16smx_t")
    {
        inv(i16smx_a, i16smx_b);
    };

    BENCHMARK("inv_rec i16smx_t")
    {
        inv_rec(i16smx_a, i16smx_b);
    };

    const auto& i32smx_a = i32smx_1.partition<m, m>();
    auto& i32smx_b = i32smx_2.partition<m, m>();

    BENCHMARK("inv i32smx_t")
    {
        inv(i32smx_a, i32smx_b);
    };

    BENCHMARK("inv_rec i32smx_t")
    {
        inv_rec(i32smx_a, i32smx_b);
    };

    const auto& i64smx_a = i64smx_1.partition<m, m>();
    auto& i64smx_b = i64smx_2.partition<m, m>();

    BENCHMARK("inv i64smx_t")
    {
        inv(i64smx_a, i64smx_b);
    };

    BENCHMARK("inv_rec i64smx_t")
    {
        inv_rec(i64smx_a, i64smx_b);
    };

    const auto& f32smx_a = f32smx_1.partition<m, m>();
    auto& f32smx_b = f32smx_2.partition<m, m>();

    BENCHMARK("inv f32smx_t")
    {
        inv(f32smx_a, f32smx_b);
    };

    BENCHMARK("inv_rec f32smx_t")
    {
        inv_rec(f32smx_a, f32smx_b);
    };

    const auto& f64smx_a = f64smx_1.partition<m, m>();
    auto& f64smx_b = f64smx_2.partition<m, m>();

    BENCHMARK("inv f64smx_t")
    {
        inv(f64smx_a, f64smx_b);
    };

    BENCHMARK("inv_rec f64smx_t")
    {
        inv_rec(f64smx_a, f64smx_b);
    };

    const auto& f128smx_a = f128smx_1.partition<m, m>();
    auto& f128smx_b = f128smx_2.partition<m, m>();

    BENCHMARK("inv f128smx_t")
    {
        inv(f128smx_a, f128smx_b);
    };

    BENCHMARK("inv_rec f128smx_t")
    {
        inv_rec(f128smx_a, f128smx_b);
    };
}
