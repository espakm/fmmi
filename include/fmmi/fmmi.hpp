#ifndef FMMI_FMMI_HPP
#define FMMI_FMMI_HPP

#include "fmmi/matrix.hpp"

namespace fmmi
{

template <typename T, uint16_t m, uint16_t n, uint16_t p,
          uint16_t y0_a, uint16_t x0_a, uint16_t stride_a,
          uint16_t y0_b, uint16_t x0_b, uint16_t stride_b,
          uint16_t y0_c, uint16_t x0_c, uint16_t stride_c>
//inline __attribute__((always_inline))
void mul_fast(const matrix<T, m, n, y0_a, x0_a, stride_a>& a,
              const matrix<T, n, p, y0_b, x0_b, stride_b>& b,
              matrix<T, m, p, y0_c, x0_c, stride_c>& c)
{
    if constexpr (m == 1 && n == 1 && p == 1)
    {
        c(0, 0) = a(0, 0) * b(0, 0);
    }
    else if constexpr (m == 1 || n == 1 || p == 1)
    {
        mul(a, b, c);
    }
    else if constexpr (m == 2 && n == 2 && p == 2)
    {
        const auto p1 = (a(0, 0) + a(1, 1)) * (b(0, 0) + b(1, 1));
        const auto p2 = (a(1, 0) + a(1, 1)) * b(0, 0);
        const auto p3 = a(0, 0) * (b(0, 1) - b(1, 1));
        const auto p4 = a(1, 1) * (b(1, 0) - b(0, 0));
        const auto p5 = (a(0, 0) + a(0, 1)) * b(1, 1);
        const auto p6 = (a(1, 0) - a(0, 0)) * (b(0, 0) + b(0, 1));
        const auto p7 = (a(0, 1) - a(1, 1)) * (b(1, 0) + b(1, 1));

        c(0, 0) = p1 + p4 - p5 + p7;
        c(0, 1) = p3 + p5;
        c(1, 0) = p2 + p4;
        c(1, 1) = p1 + p3 - p2 + p6;
    }
    else if constexpr ((m % 2) == 0 && (n % 2) == 0 && (p % 2) == 0)
    {
        constexpr uint16_t mh = m / 2;
        constexpr uint16_t nh = n / 2;
        constexpr uint16_t ph = p / 2;

        const auto& a00 = a.template partition<mh, nh, 0, 0>();
        const auto& a01 = a.template partition<mh, nh, 0, nh>();
        const auto& a10 = a.template partition<mh, nh, mh, 0>();
        const auto& a11 = a.template partition<mh, nh, mh, nh>();

        const auto& b00 = b.template partition<nh, ph, 0, 0>();
        const auto& b01 = b.template partition<nh, ph, 0, ph>();
        const auto& b10 = b.template partition<nh, ph, nh, 0>();
        const auto& b11 = b.template partition<nh, ph, nh, ph>();

        auto& c00 = c.template partition<mh, ph, 0, 0>();
        auto& c01 = c.template partition<mh, ph, 0, ph>();
        auto& c10 = c.template partition<mh, ph, mh, 0>();
        auto& c11 = c.template partition<mh, ph, mh, ph>();

        static matrix<T, mh, nh> tmp1, tmp2, tmp3, tmp4, tmp5;
        static matrix<T, nh, ph> tmp6, tmp7, tmp8, tmp9, tmp10;

        for (uint16_t i = 0; i < nh; ++i)
        {
            for (uint16_t j = 0; j < mh; ++j)
            {
                tmp1(j, i) = a00(j, i) + a11(j, i);
                tmp2(j, i) = a10(j, i) + a11(j, i);
                tmp3(j, i) = a00(j, i) + a01(j, i);
                tmp4(j, i) = a10(j, i) - a00(j, i);
                tmp5(j, i) = a01(j, i) - a11(j, i);
            }
            for (uint16_t j = 0; j < ph; ++j)
            {
                tmp6(i, j) = b00(i, j) + b11(i, j);
                tmp7(i, j) = b01(i, j) - b11(i, j);
                tmp8(i, j) = b10(i, j) - b00(i, j);
                tmp9(i, j) = b00(i, j) + b01(i, j);
                tmp10(i, j) = b10(i, j) + b11(i, j);
            }
        }

        static matrix<T, mh, ph> p1, p2, p3, p4, p5, p6, p7;
        mul_fast(tmp1, tmp6, p1);
        mul_fast(tmp2, b00, p2);
        mul_fast(a00, tmp7, p3);
        mul_fast(a11, tmp8, p4);
        mul_fast(tmp3, b11, p5);
        mul_fast(tmp4, tmp9, p6);
        mul_fast(tmp5, tmp10, p7);

        for (uint16_t i = 0; i < mh; ++i)
        {
            for (uint16_t j = 0; j < ph; ++j)
            {
                c00(i, j) = p1(i, j) + p4(i, j) - p5(i, j) + p7(i, j);
                c01(i, j) = p3(i, j) + p5(i, j);
                c10(i, j) = p2(i, j) + p4(i, j);
                c11(i, j) = p1(i, j) + p3(i, j) - p2(i, j) + p6(i, j);
            }
        }
    }
    else if constexpr ((m % 2) == 1)
    {
        const auto& a0 = a.template partition<1, n, 0, 0>();
        const auto& a1 = a.template partition<m - 1, n, 1, 0>();
        auto& c0 = c.template partition<1, p, 0, 0>();
        auto& c1 = c.template partition<m - 1, p, 1, 0>();

        mul_fast(a0, b, c0);
        mul_fast(a1, b, c1);
    }
    else if constexpr ((p % 2) == 1)
    {
        const auto& b0 = b.template partition<n, 1, 0, 0>();
        const auto& b1 = b.template partition<n, p - 1, 0, 1>();
        auto& c0 = c.template partition<m, 1, 0, 0>();
        auto& c1 = c.template partition<m, p - 1, 0, 1>();

        mul_fast(a, b0, c0);
        mul_fast(a, b1, c1);
    }
    else // if constexpr ((n % 2) == 1)
    {
        const auto& a0 = a.template partition<m, 1, 0, 0>();
        const auto& a1 = a.template partition<m, n - 1, 0, 1>();
        const auto& b0 = b.template partition<1, p, 0, 0>();
        const auto& b1 = b.template partition<n - 1, p, 1, 0>();

        matrix<T, m, p> tmp1;
        mul(a0, b0, tmp1);

        mul_fast(a1, b1, c);
        add(c, tmp1, c);
    }
}

}

#endif
