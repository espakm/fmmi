#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch.hpp>

#include <random>

#include <fmmi/dmatrix_rec.hpp>

#ifndef DMATRIX_REC_TEST_MAX_MATRIX_SIZE
#define DMATRIX_REC_TEST_MAX_MATRIX_SIZE 64
#endif

using namespace fmmi;

static const dmatrix<double> dmx1(7, 7, {
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
    8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0,
    15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0,
    22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0,
    29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0,
    36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0,
    43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0,
});

static const dmatrix<double> dmx2(7, 7, {
    -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0,
    -8.0, -9.0, -10.0, -11.0, -12.0, -13.0, -14.0,
    -15.0, -16.0, -17.0, -18.0, -19.0, -20.0, -21.0,
    -22.0, -23.0, -24.0, -25.0, -26.0, -27.0, -28.0,
    -29.0, -30.0, -31.0, -32.0, -33.0, -34.0, -35.0,
    -36.0, -37.0, -38.0, -39.0, -40.0, -41.0, -42.0,
    -43.0, -44.0, -45.0, -46.0, -47.0, -48.0, -49.0,
});

static dmatrix<double> dmx3(7, 7);
static dmatrix<double> dmx4(7, 7);

TEMPLATE_TEST_CASE_SIG("mul_rec smatrix/dmatrix equals up to 7x7", "[mul_rec][dmatrix][equals]",
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
    const auto& da = dmx1.partition(m, n);
    const auto& db = dmx2.partition(n, p);
    auto dc = dmx3.partition(m, p);
    auto dd = dmx4.partition(m, p);

    mul(da, db, dc);

    mul_rec(da, db, dd);

    CHECK(dc == dd);
}


static constexpr uint16_t S = 4096;
static i16dmx_t i16dmx_1(S, S), i16dmx_2(S, S), i16dmx_3(S, S);
static i32dmx_t i32dmx_1(S, S), i32dmx_2(S, S), i32dmx_3(S, S);
static i64dmx_t i64dmx_1(S, S), i64dmx_2(S, S), i64dmx_3(S, S);
static f32dmx_t f32dmx_1(S, S), f32dmx_2(S, S), f32dmx_3(S, S);
static f64dmx_t f64dmx_1(S, S), f64dmx_2(S, S), f64dmx_3(S, S);
static f64dmx_t f128dmx_1(S, S), f128dmx_2(S, S), f128dmx_3(S, S);

static
struct dmatric_rec_test_static_init
{
    dmatric_rec_test_static_init()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(-200, +200);

        for (uint16_t i = 0; i < S; ++i)
        {
            for (uint16_t j = 0; j < S; ++j)
            {
                i16dmx_1(i, j) = distrib(gen);
                i16dmx_2(i, j) = distrib(gen);
                i32dmx_1(i, j) = distrib(gen);
                i32dmx_2(i, j) = distrib(gen);
                i64dmx_1(i, j) = distrib(gen);
                i64dmx_2(i, j) = distrib(gen);
                f32dmx_1(i, j) = distrib(gen);
                f32dmx_2(i, j) = distrib(gen);
                f64dmx_1(i, j) = distrib(gen);
                f64dmx_2(i, j) = distrib(gen);
                f128dmx_1(i, j) = distrib(gen);
                f128dmx_2(i, j) = distrib(gen);
            }
        }
    }
} dmatric_rec_test_static_init;


TEMPLATE_TEST_CASE_SIG("mul_rec dmatrix benchmark", "[mul_rec][dmatrix][benchmark]",
                       ((uint16_t m, uint16_t n, uint16_t p), m, n, p)
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 2
                       , (1, 1, 1)
                       , (2, 2, 2)
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 4
                       , (3, 3, 3)
                       , (4, 4, 4)
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 8
                       , (7, 7, 7)
                       , (8, 8, 8)
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 16
                       , (15, 15, 15)
                       , (16, 16, 16)
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 32
                       , (31, 31, 31)
                       , (32, 32, 32)
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 64
                       , (63, 63, 63)
                       , (64, 64, 64)
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 128
                       , (127, 127, 127)
                       , (128, 128, 128)
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 256
                       , (255, 255, 255)
                       , (256, 256, 256)
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 512
                       , (511, 511, 511)
                       , (512, 512, 512)
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 1024
                       , (1023, 1023, 1023)
                       , (1024, 1024, 1024)
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 2048
                       , (2047, 2047, 2047)
                       , (2048, 2048, 2048)
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 4096
                       , (4095, 4095, 4095)
                       , (4096, 4096, 4096)
#endif
                       )
{
    const auto& i16dmx_a = i16dmx_1.partition(m, n);
    const auto& i16dmx_b = i16dmx_2.partition(n, p);
    auto i16dmx_c = i16dmx_3.partition(m, p);

    BENCHMARK("mul i16dmx_t")
    {
        mul(i16dmx_a, i16dmx_b, i16dmx_c);
    };

    BENCHMARK("mul_rec i16dmx_t")
    {
        mul_rec(i16dmx_a, i16dmx_b, i16dmx_c);
    };

    const auto& i32dmx_a = i32dmx_1.partition(m, n);
    const auto& i32dmx_b = i32dmx_2.partition(n, p);
    auto i32dmx_c = i32dmx_3.partition(m, p);

    BENCHMARK("mul i32dmx_t")
    {
        mul(i32dmx_a, i32dmx_b, i32dmx_c);
    };

    BENCHMARK("mul_rec i32dmx_t")
    {
        mul_rec(i32dmx_a, i32dmx_b, i32dmx_c);
    };

    const auto& i64dmx_a = i64dmx_1.partition(m, n);
    const auto& i64dmx_b = i64dmx_2.partition(n, p);
    auto i64dmx_c = i64dmx_3.partition(m, p);

    BENCHMARK("mul i64dmx_t")
    {
        mul(i64dmx_a, i64dmx_b, i64dmx_c);
    };

    BENCHMARK("mul_rec i64dmx_t")
    {
        mul_rec(i64dmx_a, i64dmx_b, i64dmx_c);
    };

    const auto& f32dmx_a = f32dmx_1.partition(m, n);
    const auto& f32dmx_b = f32dmx_2.partition(n, p);
    auto f32dmx_c = f32dmx_3.partition(m, p);

    BENCHMARK("mul f32dmx_t")
    {
        mul(f32dmx_a, f32dmx_b, f32dmx_c);
    };

    BENCHMARK("mul_rec f32dmx_t")
    {
        mul_rec(f32dmx_a, f32dmx_b, f32dmx_c);
    };

    const auto& f64dmx_a = f64dmx_1.partition(m, n);
    const auto& f64dmx_b = f64dmx_2.partition(n, p);
    auto f64dmx_c = f64dmx_3.partition(m, p);

    BENCHMARK("mul f64dmx_t")
    {
        mul(f64dmx_a, f64dmx_b, f64dmx_c);
    };

    BENCHMARK("mul_rec f64dmx_t")
    {
        mul_rec(f64dmx_a, f64dmx_b, f64dmx_c);
    };

    const auto& f128dmx_a = f128dmx_1.partition(m, n);
    const auto& f128dmx_b = f128dmx_2.partition(n, p);
    auto f128dmx_c = f128dmx_3.partition(m, p);

    BENCHMARK("mul f128dmx_t")
    {
        mul(f128dmx_a, f128dmx_b, f128dmx_c);
    };

    BENCHMARK("mul_rec f128dmx_t")
    {
        mul_rec(f128dmx_a, f128dmx_b, f128dmx_c);
    };
}


TEST_CASE("inv_rec dmatrix equals", "[inv_rec][dmatrix][equals]")
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
    inv_rec(a, a_inv);
    CHECK(a_inv.equals(expected_a_inv));

    f32dmx_t a_a_inv(1, 1);
    mul(a, a_inv, a_a_inv);
    CHECK(a_a_inv.equals(identity1x1));

    f32dmx_t a0(1, 1, {
        0,
    });
    f32dmx_t a0_inv(1, 1);
    CHECK_THROWS_WITH(inv_rec(a0, a0_inv), "Matrix is not invertible.");

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
    inv_rec(b, b_inv);
    CHECK(b_inv.equals(expected_b_inv, 1e-7));

    f32dmx_t b_b_inv(2, 2);
    mul(b, b_inv, b_b_inv);
    CHECK(b_b_inv.equals(identity2x2, 1e-7));

    f32dmx_t b0(2, 2, {
        2, 4,
        4, 8,
    });
    f32dmx_t b0_inv(2, 2);
    CHECK_THROWS_WITH(inv_rec(b0, b0_inv), "Matrix is not invertible.");

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
    inv_rec(c, c_inv);
    CHECK(c_inv.equals(expected_c_inv, 1e-7));

    f32dmx_t c_c_inv(4, 4);
    mul(c, c_inv, c_c_inv);
    CHECK(c_c_inv.equals(identity4x4, 1e-7));

    f32dmx_t c0(4, 4, {
        5, 6, 6, 8,
        2, 2, 2, 8,
        6, 6, 2, 8,
        -13, -14, -10, -24,
    });
    f32dmx_t c0_inv(4, 4);
    CHECK_THROWS_WITH(inv_rec(c0, c0_inv), "Matrix is not invertible.");

    f32dmx_t d(4, 4, {
        0, -3, -2, 0,
        1, -4, -2, 0,
        -3, 4, 1, 0,
        0, 0, 0, 1,
    });

    f32dmx_t expected_d_inv(4, 4, {
        4, -5, -2, 0,
        5, -6, -2, 0,
        -8, 9, 3, 0,
        0, 0, 0, 1,
    });

    f32dmx_t d_expected_d_inv(4, 4);
    mul(d, expected_d_inv, d_expected_d_inv);
    CHECK(d_expected_d_inv.equals(identity4x4, 1e-6));

    f32dmx_t d_inv(4, 4);
    inv_rec(d, d_inv);
    CHECK(d_inv.equals(expected_d_inv, 1e-5));

    f32dmx_t d_d_inv(4, 4);
    mul(d, d_inv, d_d_inv);
    CHECK(d_d_inv.equals(identity4x4, 1e-5));

    f32dmx_t e(3, 3, {
        0, -3, -2,
        1, -4, -2,
        -3, 4, 1,
    });

    f32dmx_t expected_e_inv(3, 3, {
        4, -5, -2,
        5, -6, -2,
        -8, 9, 3,
    });

    f32dmx_t identity3x3 = f32dmx_t::identity(3, 3);
    f32dmx_t e_expected_e_inv(3, 3);
    mul(e, expected_e_inv, e_expected_e_inv);
    CHECK(e_expected_e_inv.equals(identity3x3, 1e-6));

    f32dmx_t e_inv(3, 3);
    inv_rec(e, e_inv);
    CHECK(e_inv.equals(expected_e_inv, 1e-5));

    f32dmx_t e_e_inv(3, 3);
    mul(e, e_inv, e_e_inv);
    CHECK(e_e_inv.equals(identity3x3, 1e-5));
}


TEMPLATE_TEST_CASE_SIG("inv_rec dmatrix benchmark", "[inv_rec][dmatrix][benchmark]",
                       ((uint16_t m), m)
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 2
                       , 1
                       , 2
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 4
                       , 3
                       , 4
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 8
                       , 7
                       , 8
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 16
                       , 15
                       , 16
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 32
                       , 31
                       , 32
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 64
                       , 63
                       , 64
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 128
                       , 127
                       , 128
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 256
                       , 255
                       , 256
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 512
                       , 511
                       , 512
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 1024
                       , 1023
                       , 1024
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 2048
                       , 2047
                       , 2048
#endif
#if DMATRIX_REC_TEST_MAX_MATRIX_SIZE >= 4096
                       , 4095
                       , 4096
#endif
                       )
{
    const auto& i16dmx_a = i16dmx_1.partition(m, m);
    auto i16dmx_b = i16dmx_2.partition(m, m);

    BENCHMARK("inv i16dmx_t")
    {
        inv(i16dmx_a, i16dmx_b);
    };

    BENCHMARK("inv_rec i16dmx_t")
    {
        inv_rec(i16dmx_a, i16dmx_b);
    };

    const auto& i32dmx_a = i32dmx_1.partition(m, m);
    auto i32dmx_b = i32dmx_2.partition(m, m);

    BENCHMARK("inv i32dmx_t")
    {
        inv(i32dmx_a, i32dmx_b);
    };

    BENCHMARK("inv_rec i32dmx_t")
    {
        inv_rec(i32dmx_a, i32dmx_b);
    };

    const auto& i64dmx_a = i64dmx_1.partition(m, m);
    auto i64dmx_b = i64dmx_2.partition(m, m);

    BENCHMARK("inv i64dmx_t")
    {
        inv(i64dmx_a, i64dmx_b);
    };

    BENCHMARK("inv_rec i64dmx_t")
    {
        inv_rec(i64dmx_a, i64dmx_b);
    };

    const auto& f32dmx_a = f32dmx_1.partition(m, m);
    auto f32dmx_b = f32dmx_2.partition(m, m);

    BENCHMARK("inv f32dmx_t")
    {
        inv(f32dmx_a, f32dmx_b);
    };

    BENCHMARK("inv_rec f32dmx_t")
    {
        inv_rec(f32dmx_a, f32dmx_b);
    };

    const auto& f64dmx_a = f64dmx_1.partition(m, m);
    auto f64dmx_b = f64dmx_2.partition(m, m);

    BENCHMARK("inv f64dmx_t")
    {
        inv(f64dmx_a, f64dmx_b);
    };

    BENCHMARK("inv_rec f64dmx_t")
    {
        inv_rec(f64dmx_a, f64dmx_b);
    };

    const auto& f128dmx_a = f128dmx_1.partition(m, m);
    auto f128dmx_b = f128dmx_2.partition(m, m);

    BENCHMARK("inv f128dmx_t")
    {
        inv(f128dmx_a, f128dmx_b);
    };

    BENCHMARK("inv_rec f128dmx_t")
    {
        inv_rec(f128dmx_a, f128dmx_b);
    };
}
