#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch.hpp>

#include <random>

#include "fmmi/fmmi.hpp"

using namespace fmmi;

TEST_CASE("fmmi mul_fast")
{
    matrix<double, 6, 6> mx1{
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0,
        25.0, 26.0, 27.0, 28.0, 29.0, 30.0,
        31.0, 32.0, 33.0, 34.0, 35.0, 36.0,
    };
    matrix<double, 6, 6> mx2{
        -1.0, -2.0, -3.0, -4.0, -5.0, -6.0,
        -7.0, -8.0, -9.0, -10.0, -11.0, -12.0,
        -13.0, -14.0, -15.0, -16.0, -17.0, -18.0,
        -19.0, -20.0, -21.0, -22.0, -23.0, -24.0,
        -25.0, -26.0, -27.0, -28.0, -29.0, -30.0,
        -31.0, -32.0, -33.0, -34.0, -35.0, -36.0,
    };

    SECTION("2x2 times 2x2")
    {
        const auto& a = mx1.partition<2, 2>();
        const auto& b = mx2.partition<2, 2>();

        matrix<double, 2, 2> c;
        mul(a, b, c);

        matrix<double, 2, 2> d;
        mul_fast(a, b, d);

        CHECK(c == d);
    }

    SECTION("4x4 times 4x4")
    {
        const auto& a = mx1.partition<4, 4>();
        const auto& b = mx2.partition<4, 4>();

        matrix<double, 4, 4> c;
        mul(a, b, c);

        matrix<double, 4, 4> d;
        mul_fast(a, b, d);

        CHECK(c == d);
    }

    SECTION("3x3 times 3x3")
    {
        const auto& a = mx1.partition<3, 3>();
        const auto& b = mx2.partition<3, 3>();

        matrix<double, 3, 3> c;
        mul(a, b, c);

        matrix<double, 3, 3> d;
        mul_fast(a, b, d);

        CHECK(c == d);
    }
}


TEMPLATE_TEST_CASE_SIG("fmmi mul_fast up to 6x6", "[template][product][nttp]",
                       ((uint16_t m, uint16_t n, uint16_t p), m, n, p),
                       (1, 1, 1),
                       (1, 1, 2),
                       (1, 1, 3),
                       (1, 1, 4),
                       (1, 1, 5),
                       (1, 1, 6),
                       (1, 2, 1),
                       (1, 2, 2),
                       (1, 2, 3),
                       (1, 2, 4),
                       (1, 2, 5),
                       (1, 2, 6),
                       (1, 3, 1),
                       (1, 3, 2),
                       (1, 3, 3),
                       (1, 3, 4),
                       (1, 3, 5),
                       (1, 3, 6),
                       (1, 4, 1),
                       (1, 4, 2),
                       (1, 4, 3),
                       (1, 4, 4),
                       (1, 4, 5),
                       (1, 4, 6),
                       (1, 5, 1),
                       (1, 5, 2),
                       (1, 5, 3),
                       (1, 5, 4),
                       (1, 5, 5),
                       (1, 5, 6),
                       (1, 6, 1),
                       (1, 6, 2),
                       (1, 6, 3),
                       (1, 6, 4),
                       (1, 6, 5),
                       (1, 6, 6),
                       (2, 1, 1),
                       (2, 1, 2),
                       (2, 1, 3),
                       (2, 1, 4),
                       (2, 1, 5),
                       (2, 1, 6),
                       (2, 2, 1),
                       (2, 2, 2),
                       (2, 2, 3),
                       (2, 2, 4),
                       (2, 2, 5),
                       (2, 2, 6),
                       (2, 3, 1),
                       (2, 3, 2),
                       (2, 3, 3),
                       (2, 3, 4),
                       (2, 3, 5),
                       (2, 3, 6),
                       (2, 4, 1),
                       (2, 4, 2),
                       (2, 4, 3),
                       (2, 4, 4),
                       (2, 4, 5),
                       (2, 4, 6),
                       (2, 5, 1),
                       (2, 5, 2),
                       (2, 5, 3),
                       (2, 5, 4),
                       (2, 5, 5),
                       (2, 5, 6),
                       (2, 6, 1),
                       (2, 6, 2),
                       (2, 6, 3),
                       (2, 6, 4),
                       (2, 6, 5),
                       (2, 6, 6),
                       (3, 1, 1),
                       (3, 1, 2),
                       (3, 1, 3),
                       (3, 1, 4),
                       (3, 1, 5),
                       (3, 1, 6),
                       (3, 2, 1),
                       (3, 2, 2),
                       (3, 2, 3),
                       (3, 2, 4),
                       (3, 2, 5),
                       (3, 2, 6),
                       (3, 3, 1),
                       (3, 3, 2),
                       (3, 3, 3),
                       (3, 3, 4),
                       (3, 3, 5),
                       (3, 3, 6),
                       (3, 4, 1),
                       (3, 4, 2),
                       (3, 4, 3),
                       (3, 4, 4),
                       (3, 4, 5),
                       (3, 4, 6),
                       (3, 5, 1),
                       (3, 5, 2),
                       (3, 5, 3),
                       (3, 5, 4),
                       (3, 5, 5),
                       (3, 5, 6),
                       (3, 6, 1),
                       (3, 6, 2),
                       (3, 6, 3),
                       (3, 6, 4),
                       (3, 6, 5),
                       (3, 6, 6),
                       (4, 1, 1),
                       (4, 1, 2),
                       (4, 1, 3),
                       (4, 1, 4),
                       (4, 1, 5),
                       (4, 1, 6),
                       (4, 2, 1),
                       (4, 2, 2),
                       (4, 2, 3),
                       (4, 2, 4),
                       (4, 2, 5),
                       (4, 2, 6),
                       (4, 3, 1),
                       (4, 3, 2),
                       (4, 3, 3),
                       (4, 3, 4),
                       (4, 3, 5),
                       (4, 3, 6),
                       (4, 4, 1),
                       (4, 4, 2),
                       (4, 4, 3),
                       (4, 4, 4),
                       (4, 4, 5),
                       (4, 4, 6),
                       (4, 5, 1),
                       (4, 5, 2),
                       (4, 5, 3),
                       (4, 5, 4),
                       (4, 5, 5),
                       (4, 5, 6),
                       (4, 6, 1),
                       (4, 6, 2),
                       (4, 6, 3),
                       (4, 6, 4),
                       (4, 6, 5),
                       (4, 6, 6),
                       (5, 1, 1),
                       (5, 1, 2),
                       (5, 1, 3),
                       (5, 1, 4),
                       (5, 1, 5),
                       (5, 1, 6),
                       (5, 2, 1),
                       (5, 2, 2),
                       (5, 2, 3),
                       (5, 2, 4),
                       (5, 2, 5),
                       (5, 2, 6),
                       (5, 3, 1),
                       (5, 3, 2),
                       (5, 3, 3),
                       (5, 3, 4),
                       (5, 3, 5),
                       (5, 3, 6),
                       (5, 4, 1),
                       (5, 4, 2),
                       (5, 4, 3),
                       (5, 4, 4),
                       (5, 4, 5),
                       (5, 4, 6),
                       (5, 5, 1),
                       (5, 5, 2),
                       (5, 5, 3),
                       (5, 5, 4),
                       (5, 5, 5),
                       (5, 5, 6),
                       (5, 6, 1),
                       (5, 6, 2),
                       (5, 6, 3),
                       (5, 6, 4),
                       (5, 6, 5),
                       (5, 6, 6),
                       (6, 1, 1),
                       (6, 1, 2),
                       (6, 1, 3),
                       (6, 1, 4),
                       (6, 1, 5),
                       (6, 1, 6),
                       (6, 2, 1),
                       (6, 2, 2),
                       (6, 2, 3),
                       (6, 2, 4),
                       (6, 2, 5),
                       (6, 2, 6),
                       (6, 3, 1),
                       (6, 3, 2),
                       (6, 3, 3),
                       (6, 3, 4),
                       (6, 3, 5),
                       (6, 3, 6),
                       (6, 4, 1),
                       (6, 4, 2),
                       (6, 4, 3),
                       (6, 4, 4),
                       (6, 4, 5),
                       (6, 4, 6),
                       (6, 5, 1),
                       (6, 5, 2),
                       (6, 5, 3),
                       (6, 5, 4),
                       (6, 5, 5),
                       (6, 5, 6),
                       (6, 6, 1),
                       (6, 6, 2),
                       (6, 6, 3),
                       (6, 6, 4),
                       (6, 6, 5),
                       (6, 6, 6))
{
    matrix<double, 6, 6> mx1{
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0,
        25.0, 26.0, 27.0, 28.0, 29.0, 30.0,
        31.0, 32.0, 33.0, 34.0, 35.0, 36.0,
    };
    matrix<double, 6, 6> mx2{
        -1.0, -2.0, -3.0, -4.0, -5.0, -6.0,
        -7.0, -8.0, -9.0, -10.0, -11.0, -12.0,
        -13.0, -14.0, -15.0, -16.0, -17.0, -18.0,
        -19.0, -20.0, -21.0, -22.0, -23.0, -24.0,
        -25.0, -26.0, -27.0, -28.0, -29.0, -30.0,
        -31.0, -32.0, -33.0, -34.0, -35.0, -36.0,
    };

    const auto& a = mx1.partition<m, n>();
    const auto& b = mx2.partition<n, p>();

    matrix<double, m, p> c;
    mul(a, b, c);

    matrix<double, m, p> d;
    mul_fast(a, b, d);

    CHECK(c == d);
}


static constexpr uint16_t S = 256;
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


TEMPLATE_TEST_CASE_SIG("fmmi mul_fast benchmark int16_t", "",
                       ((uint16_t m, uint16_t n, uint16_t p), m, n, p),
                       (1, 1, 1),
                       (2, 2, 2),
                       (4, 4, 4),
                       (8, 8, 8),
                       (16, 16, 16),
                       (32, 32, 32),
                       (64, 64, 64),
                       (128, 128, 128),
                       (255, 255, 255),
                       (256, 256, 256)
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
