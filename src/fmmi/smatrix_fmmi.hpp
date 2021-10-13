#ifndef FMMI_SMATRIX_FMMI_HPP
#define FMMI_SMATRIX_FMMI_HPP

#include "fmmi/smatrix.hpp"

namespace fmmi
{

template <typename T, uint16_t m, uint16_t n, uint16_t p,
          uint16_t y0_a, uint16_t x0_a, uint16_t stride_a,
          uint16_t y0_b, uint16_t x0_b, uint16_t stride_b,
          uint16_t y0_c, uint16_t x0_c, uint16_t stride_c>
inline
void mul_rec(const smatrix<T, m, n, y0_a, x0_a, stride_a>& a,
             const smatrix<T, n, p, y0_b, x0_b, stride_b>& b,
             smatrix<T, m, p, y0_c, x0_c, stride_c>& c)
{
    if constexpr (m == 2 && n == 2 && p == 2)
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

        smatrix<T, mh, nh> tmp1, tmp2, tmp3, tmp4, tmp5;
        smatrix<T, nh, ph> tmp6, tmp7, tmp8, tmp9, tmp10;

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

        smatrix<T, mh, ph> p1, p2, p3, p4, p5, p6, p7;

        mul_rec(tmp1, tmp6, p1);
        mul_rec(tmp2, b00, p2);
        mul_rec(a00, tmp7, p3);
        mul_rec(a11, tmp8, p4);
        mul_rec(tmp3, b11, p5);
        mul_rec(tmp4, tmp9, p6);
        mul_rec(tmp5, tmp10, p7);

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
        if constexpr (m != 1)
        {
            const auto& a1 = a.template partition<m - 1, n, 1, 0>();
            auto& c1 = c.template partition<m - 1, p, 1, 0>();
            mul_rec(a1, b, c1);
        }

        for (uint16_t j = 0; j < p; ++j)
        {
            T sum{};
            for (uint16_t k = 0; k < n; ++k)
            {
                sum += a(0, k) * b(k, j);
            }
            c(0, j) = sum;
        }
    }
    else if constexpr ((p % 2) == 1)
    {
        if constexpr (p != 1)
        {
            const auto& b1 = b.template partition<n, p - 1, 0, 1>();
            auto& c1 = c.template partition<m, p - 1, 0, 1>();
            mul_rec(a, b1, c1);
        }

        for (uint16_t i = 0; i < m; ++i)
        {
            T sum{};
            for (uint16_t k = 0; k < n; ++k)
            {
                sum += a(i, k) * b(k, 0);
            }
            c(i, 0) = sum;
        }
    }
    else /// (n % 2) == 1
    {
        if constexpr (n != 1)
        {
            const auto& a1 = a.template partition<m, n - 1, 0, 1>();
            const auto& b1 = b.template partition<n - 1, p, 1, 0>();
            mul_rec(a1, b1, c);

            for (uint16_t i = 0; i < m; ++i)
            {
                for (uint16_t j = 0; j < p; ++j)
                {
                    c(i, j) += a(i, 0) * b(0, j);
                }
            }
        }
        else
        {
            for (uint16_t i = 0; i < m; ++i)
            {
                for (uint16_t j = 0; j < p; ++j)
                {
                    c(i, j) = a(i, 0) * b(0, j);
                }
            }
        }
    }
}


template <typename T, uint16_t m,
          uint16_t y0_a, uint16_t x0_a, uint16_t stride_a,
          uint16_t y0_ainv, uint16_t x0_ainv, uint16_t stride_ainv>
inline
void inv_rec(const smatrix<T, m, m, y0_a, x0_a, stride_a>& a,
              smatrix<T, m, m, y0_ainv, x0_ainv, stride_ainv>& ainv)
{
    if constexpr (m == 1)
    {
        ainv(0, 0) = 1 / a(0, 0);
    }
    else if constexpr (m == 2)
    {
        const T coeff = 1.0 / (a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0));
        ainv(0, 0) = a(1, 1) * coeff;
        ainv(0, 1) = -a(0, 1) * coeff;
        ainv(1, 0) = -a(1, 0) * coeff;
        ainv(1, 1) = a(0, 0) * coeff;
    }
    else
    {
        constexpr uint16_t n = (m % 2) == 0 ? m / 2 : 1;
        constexpr uint16_t p = m - n;

        const auto& a00 = a.template partition<n, n, 0, 0>();
        const auto& a01 = a.template partition<n, p, 0, n>();
        const auto& a10 = a.template partition<p, n, n, 0>();
        const auto& a11 = a.template partition<p, p, n, n>();

        auto& ainv00 = ainv.template partition<n, n, 0, 0>();
        auto& ainv01 = ainv.template partition<n, p, 0, n>();
        auto& ainv10 = ainv.template partition<p, n, n, 0>();
        auto& ainv11 = ainv.template partition<p, p, n, n>();

        smatrix<T, m, m> tmp;
        auto& tmp00 = tmp.template partition<n, n, 0, 0>();
        auto& tmp01 = tmp.template partition<n, p, 0, n>();
        auto& tmp10 = tmp.template partition<p, n, n, 0>();
        auto& tmp11 = tmp.template partition<p, p, n, n>();

        inv_rec(a00, tmp00);
        mul_rec(a10, tmp00, tmp10);
        mul_rec(tmp10, a01, tmp11);
        sub(a11, tmp11);
        inv_rec(tmp11, ainv11);

        mul_rec(tmp00, a01, tmp01);
        mul_rec(tmp01, ainv11, ainv01);

        mul_rec(ainv01, tmp10, ainv00);
        add(tmp00, ainv00);

        mul_rec(ainv11, tmp10, ainv10);

        for (uint16_t i = 0; i < n; ++i)
        {
            for (uint16_t j = 0; j < p; ++j)
            {
                ainv01(i, j) *= -1;
                ainv10(j, i) *= -1;
            }
        }
    }
}

}

#endif
