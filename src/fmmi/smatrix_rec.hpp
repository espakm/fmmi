#ifndef FMMI_SMATRIX_REC_HPP
#define FMMI_SMATRIX_REC_HPP

#include "fmmi/smatrix.hpp"

namespace fmmi
{

template <typename T, uint16_t m, uint16_t n, uint16_t p,
          uint16_t a_y0, uint16_t a_x0, uint16_t a_stride,
          uint16_t b_y0, uint16_t b_x0, uint16_t b_stride,
          uint16_t c_y0, uint16_t c_x0, uint16_t c_stride>
inline
void mul_rec(const smatrix<T, m, n, a_y0, a_x0, a_stride>& a,
             const smatrix<T, n, p, b_y0, b_x0, b_stride>& b,
             smatrix<T, m, p, c_y0, c_x0, c_stride>& c)
{
    constexpr uint16_t max_size = m < n ? (n < p ? p : n) : (m < p ? p : m);
    constexpr uint16_t pm = pow_2_upper_bound(max_size);

    if constexpr (m == pm && n == pm && p == pm)
    {
        mul_rec_pow_2(a, b, c);
    }
    else if constexpr (m == pm && n == pm)
    {
        smatrix<T, pm, pm> b2 = b;
        for (uint16_t i = 0; i < pm; ++i)
        {
            std::fill_n(&b2(i, p), pm - p, 0);
        }

        smatrix<T, pm, pm> c2;
        mul_rec_pow_2(a, b2, c2);
        c = c2.template partition<m, p>();
    }
    else if constexpr (n == pm && p == pm)
    {
        smatrix<T, pm, pm> a2 = a;
        for (uint16_t i = m; i < pm; ++i)
        {
            std::fill_n(&a2(i, 0), pm, 0);
        }

        smatrix<T, pm, pm> c2;
        mul_rec_pow_2(a2, b, c2);
        c = c2.template partition<m, p>();
    }
    else
    {
        smatrix<T, pm, pm> a2 = a;
        smatrix<T, pm, pm> b2 = b;
        uint16_t i;
        for (i = 0; i < m; ++i)
        {
            std::fill_n(&a2(i, n), pm - n, 0);
        }
        for (; i < pm; ++i)
        {
            std::fill_n(&a2(i, 0), pm, 0);
        }
        for (i = 0; i < n; ++i)
        {
            std::fill_n(&b2(i, p), pm - p, 0);
        }
        for (; i < pm; ++i)
        {
            std::fill_n(&b2(i, 0), pm, 0);
        }

        if constexpr (m == pm && p == pm)
        {
            mul_rec_pow_2(a2, b2, c);
        }
        else
        {
            smatrix<T, pm, pm> c2;
            mul_rec_pow_2(a2, b2, c2);
            c = c2.template partition<m, p>();
        }
    }
}


template <typename T, uint16_t m,
          uint16_t a_y0, uint16_t a_x0, uint16_t a_stride,
          uint16_t b_y0, uint16_t b_x0, uint16_t b_stride,
          uint16_t c_y0, uint16_t c_x0, uint16_t c_stride>
inline
void mul_rec_pow_2(const smatrix<T, m, m, a_y0, a_x0, a_stride>& a,
                   const smatrix<T, m, m, b_y0, b_x0, b_stride>& b,
                   smatrix<T, m, m, c_y0, c_x0, c_stride>& c)
{
    static_assert(m == pow_2_upper_bound(m),
            "mul_rec_pow_2 can only be applied on quadratic matrices of size power of 2.");

    if constexpr (m == 1)
    {
        c(0, 0) = a(0, 0) * b(0, 0);
    }
    else if constexpr (m == 2)
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
    else
    {
        constexpr uint16_t mh = m / 2;

        const auto& a00 = a.template partition<mh, mh, 0, 0>();
        const auto& a01 = a.template partition<mh, mh, 0, mh>();
        const auto& a10 = a.template partition<mh, mh, mh, 0>();
        const auto& a11 = a.template partition<mh, mh, mh, mh>();

        const auto& b00 = b.template partition<mh, mh, 0, 0>();
        const auto& b01 = b.template partition<mh, mh, 0, mh>();
        const auto& b10 = b.template partition<mh, mh, mh, 0>();
        const auto& b11 = b.template partition<mh, mh, mh, mh>();

        auto& c00 = c.template partition<mh, mh, 0, 0>();
        auto& c01 = c.template partition<mh, mh, 0, mh>();
        auto& c10 = c.template partition<mh, mh, mh, 0>();
        auto& c11 = c.template partition<mh, mh, mh, mh>();

        smatrix<T, mh, mh> tmp1, tmp2, tmp3, tmp4, tmp5;
        smatrix<T, mh, mh> tmp6, tmp7, tmp8, tmp9, tmp10;

        for (uint16_t i = 0; i < mh; ++i)
        {
            for (uint16_t j = 0; j < mh; ++j)
            {
                tmp1(j, i) = a00(j, i) + a11(j, i);
                tmp2(j, i) = a10(j, i) + a11(j, i);
                tmp3(j, i) = a00(j, i) + a01(j, i);
                tmp4(j, i) = a10(j, i) - a00(j, i);
                tmp5(j, i) = a01(j, i) - a11(j, i);

                tmp6(i, j) = b00(i, j) + b11(i, j);
                tmp7(i, j) = b01(i, j) - b11(i, j);
                tmp8(i, j) = b10(i, j) - b00(i, j);
                tmp9(i, j) = b00(i, j) + b01(i, j);
                tmp10(i, j) = b10(i, j) + b11(i, j);
            }
        }

        smatrix<T, mh, mh> p1, p2, p3, p4, p5, p6, p7;

        mul_rec_pow_2(tmp1, tmp6, p1);
        mul_rec_pow_2(tmp2, b00, p2);
        mul_rec_pow_2(a00, tmp7, p3);
        mul_rec_pow_2(a11, tmp8, p4);
        mul_rec_pow_2(tmp3, b11, p5);
        mul_rec_pow_2(tmp4, tmp9, p6);
        mul_rec_pow_2(tmp5, tmp10, p7);

        for (uint16_t i = 0; i < mh; ++i)
        {
            for (uint16_t j = 0; j < mh; ++j)
            {
                c00(i, j) = p1(i, j) + p4(i, j) - p5(i, j) + p7(i, j);
                c01(i, j) = p3(i, j) + p5(i, j);
                c10(i, j) = p2(i, j) + p4(i, j);
                c11(i, j) = p1(i, j) + p3(i, j) - p2(i, j) + p6(i, j);
            }
        }
    }
}


template <typename T, uint16_t m, uint16_t n, uint16_t p,
          uint16_t a_y0, uint16_t a_x0, uint16_t a_stride,
          uint16_t b_y0, uint16_t b_x0, uint16_t b_stride,
          uint16_t c_y0, uint16_t c_x0, uint16_t c_stride>
inline
void mul_rec_part(const smatrix<T, m, n, a_y0, a_x0, a_stride>& a,
                  const smatrix<T, n, p, b_y0, b_x0, b_stride>& b,
                  smatrix<T, m, p, c_y0, c_x0, c_stride>& c)
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

        mul_rec_part(tmp1, tmp6, p1);
        mul_rec_part(tmp2, b00, p2);
        mul_rec_part(a00, tmp7, p3);
        mul_rec_part(a11, tmp8, p4);
        mul_rec_part(tmp3, b11, p5);
        mul_rec_part(tmp4, tmp9, p6);
        mul_rec_part(tmp5, tmp10, p7);

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
            mul_rec_part(a1, b, c1);
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
            mul_rec_part(a, b1, c1);
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
            mul_rec_part(a1, b1, c);

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
          uint16_t a_y0, uint16_t a_x0, uint16_t a_stride,
          uint16_t ainv_y0, uint16_t ainv_x0, uint16_t ainv_stride>
inline
void inv_rec(const smatrix<T, m, m, a_y0, a_x0, a_stride>& a,
              smatrix<T, m, m, ainv_y0, ainv_x0, ainv_stride>& ainv)
{
    constexpr uint16_t p = pow_2_upper_bound(m);

    if constexpr (p == m)
    {
        inv_rec_pow_2(a, ainv);
    }
    else
    {
        /// Add 0 and 1 elements from the right and the bottom, like in the
        /// unity matrix. Matrix is copied that is inefficient.
        smatrix<T, p, p> a2 = a;
        for (uint16_t i = m; i < p; ++i)
        {
            uint16_t j = 0;
            for (; j < m; ++j)
            {
                a2(i, j) = 0;
                a2(j, i) = 0;
            }
            for (; j < p; ++j)
            {
                a2(i, j) = i == j;
            }
        }

        smatrix<T, p, p> ainv2;
        inv_rec_pow_2(a2, ainv2);

        /// The result is in the mxm top-left partition of the new matrix.
        /// Sadly, it needs to be copied back to the output argument that
        /// is inefficient.
        ainv = ainv2.template partition<m, m>();
    }
}


template <typename T, uint16_t m,
          uint16_t a_y0, uint16_t a_x0, uint16_t a_stride,
          uint16_t ainv_y0, uint16_t ainv_x0, uint16_t ainv_stride>
inline
void inv_rec_pow_2(const smatrix<T, m, m, a_y0, a_x0, a_stride>& a,
                   smatrix<T, m, m, ainv_y0, ainv_x0, ainv_stride>& ainv)
{
    if constexpr (m == 1)
    {
        const T f = a(0, 0);
        if (f == 0)
        {
            throw std::logic_error("Matrix is not invertible.");
        }
        ainv(0, 0) = 1 / f;
    }
    else if constexpr (m == 2)
    {
        const T f = a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);
        if (f == 0)
        {
            throw std::logic_error("Matrix is not invertible.");
        }
        ainv(0, 0) = a(1, 1) / f;
        ainv(0, 1) = -a(0, 1) / f;
        ainv(1, 0) = -a(1, 0) / f;
        ainv(1, 1) = a(0, 0) / f;
    }
    else
    {
        constexpr uint16_t mh = m / 2;
        static_assert (2 * mh == m, "Matrix size is not quadratic and power of 2.");

        const auto& a00 = a.template partition<mh, mh, 0, 0>();
        const auto& a01 = a.template partition<mh, mh, 0, mh>();
        const auto& a10 = a.template partition<mh, mh, mh, 0>();
        const auto& a11 = a.template partition<mh, mh, mh, mh>();

        auto& ainv00 = ainv.template partition<mh, mh, 0, 0>();
        auto& ainv01 = ainv.template partition<mh, mh, 0, mh>();
        auto& ainv10 = ainv.template partition<mh, mh, mh, 0>();
        auto& ainv11 = ainv.template partition<mh, mh, mh, mh>();

        smatrix<T, m, m> tmp;
        auto& tmp00 = tmp.template partition<mh, mh, 0, 0>();
        auto& tmp01 = tmp.template partition<mh, mh, 0, mh>();
        auto& tmp10 = tmp.template partition<mh, mh, mh, 0>();
        auto& tmp11 = tmp.template partition<mh, mh, mh, mh>();

        inv_rec_pow_2(a00, tmp00);
        mul_rec_pow_2(a10, tmp00, tmp10);
        mul_rec_pow_2(tmp10, a01, tmp11);
        sub(a11, tmp11);
        inv_rec_pow_2(tmp11, ainv11);

        mul_rec_pow_2(tmp00, a01, tmp01);
        mul_rec_pow_2(tmp01, ainv11, ainv01);

        mul_rec_pow_2(ainv01, tmp10, ainv00);
        add(tmp00, ainv00);

        mul_rec_pow_2(ainv11, tmp10, ainv10);

        for (uint16_t i = 0; i < mh; ++i)
        {
            for (uint16_t j = 0; j < mh; ++j)
            {
                ainv01(i, j) *= -1;
                ainv10(j, i) *= -1;
            }
        }
    }
}

}

#endif
