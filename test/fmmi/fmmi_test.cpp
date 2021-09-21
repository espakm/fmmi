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
        const auto a = mx1.partition<0, 0, 2, 2>();
        const auto b = mx2.partition<0, 0, 2, 2>();

        matrix<double, 2, 2> c;
        mul(a, b, c);

        matrix<double, 2, 2> d;
        mul_fast(a, b, d);

        CHECK(c == d);
    }

    SECTION("4x4 times 4x4")
    {
        const auto a = mx1.partition<0, 0, 4, 4>();
        const auto b = mx2.partition<0, 0, 4, 4>();

        matrix<double, 4, 4> c;
        mul(a, b, c);

        matrix<double, 4, 4> d;
        mul_fast(a, b, d);

        CHECK(c == d);
    }

    SECTION("3x3 times 3x3")
    {
        const auto a = mx1.partition<0, 0, 3, 3>();
        const auto b = mx2.partition<0, 0, 3, 3>();

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

    const auto a = mx1.partition<0, 0, m, n>();
    const auto b = mx2.partition<0, 0, n, p>();

    matrix<double, m, p> c;
    mul(a, b, c);

    matrix<double, m, p> d;
    mul_fast(a, b, d);

    CHECK(c == d);
}


TEST_CASE("fmmi mul_fast benchmark")
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(-200, +200);

    constexpr uint16_t s = 256;
    matrix<double, s, s> mx1, mx2, mx3;

    for (uint16_t i = 0; i < s; ++i)
    {
        for (uint16_t j = 0; j < s; ++j)
        {
            mx1(i, j) = distrib(gen);
            mx2(i, j) = distrib(gen);
        }
    }

    SECTION("2x2")
    {
        constexpr uint16_t n = 2;
        auto a = mx1.template partition<0, 0, n, n>();
        auto b = mx2.template partition<0, 0, n, n>();
        auto c = mx3.template partition<0, 0, n, n>();

        BENCHMARK("mul")
        {
            mul(a, b, c);
        };

        BENCHMARK("mul_fast")
        {
            mul_fast(a, b, c);
        };
    }

    SECTION("4x4")
    {
        constexpr uint16_t n = 4;
        auto a = mx1.template partition<0, 0, n, n>();
        auto b = mx2.template partition<0, 0, n, n>();
        auto c = mx3.template partition<0, 0, n, n>();

        BENCHMARK("mul")
        {
            mul(a, b, c);
        };

        BENCHMARK("mul_fast")
        {
            mul_fast(a, b, c);
        };
    }

    SECTION("8x8")
    {
        constexpr uint16_t n = 8;
        auto a = mx1.template partition<0, 0, n, n>();
        auto b = mx2.template partition<0, 0, n, n>();
        auto c = mx3.template partition<0, 0, n, n>();

        BENCHMARK("mul")
        {
            mul(a, b, c);
        };

        BENCHMARK("mul_fast")
        {
            mul_fast(a, b, c);
        };
    }

    SECTION("16x16")
    {
        constexpr uint16_t n = 16;
        auto a = mx1.template partition<0, 0, n, n>();
        auto b = mx2.template partition<0, 0, n, n>();
        auto c = mx3.template partition<0, 0, n, n>();

        BENCHMARK("mul")
        {
            mul(a, b, c);
        };

        BENCHMARK("mul_fast")
        {
            mul_fast(a, b, c);
        };
    }
}


TEST_CASE("fmmi mul_fast benchmark int")
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(-200, +200);

    constexpr uint16_t s = 256;
    matrix<int, s, s> mx1, mx2, mx3;

    for (uint16_t i = 0; i < s; ++i)
    {
        for (uint16_t j = 0; j < s; ++j)
        {
            mx1(i, j) = distrib(gen);
            mx2(i, j) = distrib(gen);
        }
    }

    SECTION("2x2")
    {
        constexpr uint16_t n = 2;
        auto a = mx1.template partition<0, 0, n, n>();
        auto b = mx2.template partition<0, 0, n, n>();
        auto c = mx3.template partition<0, 0, n, n>();

        BENCHMARK("mul")
        {
            mul(a, b, c);
        };

        BENCHMARK("mul_fast")
        {
            mul_fast(a, b, c);
        };
    }

    SECTION("4x4")
    {
        constexpr uint16_t n = 4;
        auto a = mx1.template partition<0, 0, n, n>();
        auto b = mx2.template partition<0, 0, n, n>();
        auto c = mx3.template partition<0, 0, n, n>();

        BENCHMARK("mul")
        {
            mul(a, b, c);
        };

        BENCHMARK("mul_fast")
        {
            mul_fast(a, b, c);
        };
    }

    SECTION("8x8")
    {
        constexpr uint16_t n = 8;
        auto a = mx1.template partition<0, 0, n, n>();
        auto b = mx2.template partition<0, 0, n, n>();
        auto c = mx3.template partition<0, 0, n, n>();

        BENCHMARK("mul")
        {
            mul(a, b, c);
        };

        BENCHMARK("mul_fast")
        {
            mul_fast(a, b, c);
        };
    }

    SECTION("16x16")
    {
        constexpr uint16_t n = 16;
        auto a = mx1.template partition<0, 0, n, n>();
        auto b = mx2.template partition<0, 0, n, n>();
        auto c = mx3.template partition<0, 0, n, n>();

        BENCHMARK("mul")
        {
            mul(a, b, c);
        };

        BENCHMARK("mul_fast")
        {
            mul_fast(a, b, c);
        };
    }
}
