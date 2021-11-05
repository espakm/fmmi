#ifndef FMMI_DMATRIX_REC_HPP
#define FMMI_DMATRIX_REC_HPP

#include "fmmi/common.hpp"
#include "fmmi/dmatrix.hpp"

namespace fmmi
{

template <typename T, bool a_managed, bool b_managed, bool c_managed>
inline
void mul_rec(const dmatrix<T, a_managed>& a,
             const dmatrix<T, b_managed>& b,
             dmatrix<T, c_managed>& c)
{
    assert(a.width() == b.height()
           && a.height() == c.height()
           && b.width() == c.width());

    const uint16_t m = a.height();
    const uint16_t n = a.width();
    const uint16_t p = b.width();

    const uint16_t max_size = m < n ? (n < p ? p : n) : (m < p ? p : m);
    const uint16_t pm = pow_2_upper_bound(max_size);

    if (m == pm && n == pm && p == pm)
    {
        mul_rec_pow_2(a, b, c);
    }
    else if (m == pm && n == pm)
    {
        dmatrix<T> b2(pm, pm);
        b2.assign(b);
        for (uint16_t i = 0; i < pm; ++i)
        {
            std::fill_n(&b2(i, p), pm - p, 0);
        }

        dmatrix<T> c2(pm, pm);
        mul_rec_pow_2(a, b2, c2);
        c.assign(c2.partition(m, p));
    }
    else if (n == pm && p == pm)
    {
        dmatrix<T> a2(pm, pm);
        a2.assign(a);
        for (uint16_t i = m; i < pm; ++i)
        {
            std::fill_n(&a2(i, 0), pm, 0);
        }

        dmatrix<T> c2(pm, pm);
        mul_rec_pow_2(a2, b, c2);
        c.assign(c2.partition(m, p));
    }
    else
    {
        dmatrix<T> a2(pm, pm);
        dmatrix<T> b2(pm, pm);
        a2.assign(a);
        b2.assign(b);
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

        if (m == pm && p == pm)
        {
            mul_rec_pow_2(a2, b2, c);
        }
        else
        {
            dmatrix<T> c2(pm, pm);
            mul_rec_pow_2(a2, b2, c2);
            c.assign(c2.partition(m, p));
        }
    }
}




template <typename T, bool a_managed, bool b_managed, bool c_managed>
inline
void mul_rec_pow_2(const dmatrix<T, a_managed>& a,
                   const dmatrix<T, b_managed>& b,
                   dmatrix<T, c_managed>& c)
{
    assert(a.height() == a.width()
           && b.height() == b.width()
           && c.height() == c.width()
           && a.height() == b.height()
           && a.height() == c.height());

    const uint16_t m = a.height();

    if (m == 1)
    {
        c(0, 0) = a(0, 0) * b(0, 0);
    }
    else if (m == 2)
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
        const uint16_t mh = m / 2;
        assert(2 * mh == m);

        const auto& a00 = a.partition(mh, mh, 0, 0);
        const auto& a01 = a.partition(mh, mh, 0, mh);
        const auto& a10 = a.partition(mh, mh, mh, 0);
        const auto& a11 = a.partition(mh, mh, mh, mh);

        const auto& b00 = b.partition(mh, mh, 0, 0);
        const auto& b01 = b.partition(mh, mh, 0, mh);
        const auto& b10 = b.partition(mh, mh, mh, 0);
        const auto& b11 = b.partition(mh, mh, mh, mh);

        auto c00 = c.partition(mh, mh, 0, 0);
        auto c01 = c.partition(mh, mh, 0, mh);
        auto c10 = c.partition(mh, mh, mh, 0);
        auto c11 = c.partition(mh, mh, mh, mh);

        dmatrix<T> tmp1(mh, mh), tmp2(mh, mh), tmp3(mh, mh), tmp4(mh, mh), tmp5(mh, mh);
        dmatrix<T> tmp6(mh, mh), tmp7(mh, mh), tmp8(mh, mh), tmp9(mh, mh), tmp10(mh, mh);

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

        dmatrix<T> p1(mh, mh), p2(mh, mh), p3(mh, mh), p4(mh, mh), p5(mh, mh), p6(mh, mh), p7(mh, mh);
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


template <typename T, bool a_managed, bool b_managed, bool c_managed>
inline
void mul_rec_part(const dmatrix<T, a_managed>& a,
             const dmatrix<T, b_managed>& b,
             dmatrix<T, c_managed>& c)
{
    assert(a.width() == b.height()
           && a.height() == c.height()
           && b.width() == c.width());

    const uint16_t m = a.height();
    const uint16_t n = a.width();
    const uint16_t p = b.width();

    if (m == 2 && n == 2 && p == 2)
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
    else if ((m % 2) == 0 && (n % 2) == 0 && (p % 2) == 0)
    {
        uint16_t mh = m / 2;
        uint16_t nh = n / 2;
        uint16_t ph = p / 2;

        const auto& a00 = a.partition(mh, nh, 0, 0);
        const auto& a01 = a.partition(mh, nh, 0, nh);
        const auto& a10 = a.partition(mh, nh, mh, 0);
        const auto& a11 = a.partition(mh, nh, mh, nh);

        const auto& b00 = b.partition(nh, ph, 0, 0);
        const auto& b01 = b.partition(nh, ph, 0, ph);
        const auto& b10 = b.partition(nh, ph, nh, 0);
        const auto& b11 = b.partition(nh, ph, nh, ph);

        auto c00 = c.partition(mh, ph, 0, 0);
        auto c01 = c.partition(mh, ph, 0, ph);
        auto c10 = c.partition(mh, ph, mh, 0);
        auto c11 = c.partition(mh, ph, mh, ph);

        dmatrix<T> tmp1(mh, nh), tmp2(mh, nh), tmp3(mh, nh), tmp4(mh, nh), tmp5(mh, nh);
        dmatrix<T> tmp6(nh, ph), tmp7(nh, ph), tmp8(nh, ph), tmp9(nh, ph), tmp10(nh, ph);

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

        dmatrix<T> p1(mh, ph), p2(mh, ph), p3(mh, ph), p4(mh, ph), p5(mh, ph), p6(mh, ph), p7(mh, ph);
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
    else if ((m % 2) == 1)
    {
        if (m != 1)
        {
            const auto& a1 = a.partition(m - 1, n, 1, 0);
            auto c1 = c.partition(m - 1, p, 1, 0);
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
    else if ((p % 2) == 1)
    {
        if (p != 1)
        {
            const auto& b1 = b.partition(n, p - 1, 0, 1);
            auto c1 = c.partition(m, p - 1, 0, 1);
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
        if (n != 1)
        {
            const auto& a1 = a.partition(m, n - 1, 0, 1);
            const auto& b1 = b.partition(n - 1, p, 1, 0);
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


template <typename T, bool a_managed, bool ainv_managed>
inline
void inv_rec(const dmatrix<T, a_managed>& a, dmatrix<T, ainv_managed>& ainv)
{
    assert(a.height() == a.width()
           && a.height() == ainv.height()
           && ainv.height() == ainv.width());

    const uint16_t m = a.height();
    const uint16_t p = pow_2_upper_bound(m);

    if (p == m)
    {
        inv_rec_pow_2(a, ainv);
    }
    else
    {
        /// Add 0 and 1 elements from the right and the bottom, like in the
        /// unity matrix. Matrix is copied that is inefficient.
        dmatrix<T> a2(p, p);
        a2.assign(a);
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

        dmatrix<T> ainv2(p, p);
        inv_rec_pow_2(a2, ainv2);

        /// The result is in the mxm top-left partition of the new matrix.
        /// Sadly, it needs to be copied back to the output argument that
        /// is inefficient.
        ainv.assign(ainv2.partition(m, m));
    }
}


template <typename T, bool a_managed, bool ainv_managed>
inline
void inv_rec_pow_2(const dmatrix<T, a_managed>& a, dmatrix<T, ainv_managed>& ainv)
{
    assert(a.height() == a.width()
           && a.height() == ainv.height()
           && ainv.height() == ainv.width());

    const uint16_t m = a.height();

    if (m == 1)
    {
        const T f = a(0, 0);
        if (f == 0)
        {
            throw std::logic_error("Matrix is not invertible.");
        }
        ainv(0, 0) = 1 / f;
    }
    else if (m == 2)
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
        const uint16_t mh = m / 2;
        assert(2 * mh == m);

        const auto& a00 = a.partition(mh, mh, 0, 0);
        const auto& a01 = a.partition(mh, mh, 0, mh);
        const auto& a10 = a.partition(mh, mh, mh, 0);
        const auto& a11 = a.partition(mh, mh, mh, mh);

        auto ainv00 = ainv.partition(mh, mh, 0, 0);
        auto ainv01 = ainv.partition(mh, mh, 0, mh);
        auto ainv10 = ainv.partition(mh, mh, mh, 0);
        auto ainv11 = ainv.partition(mh, mh, mh, mh);

        dmatrix<T> tmp(m, m);
        auto tmp00 = tmp.partition(mh, mh, 0, 0);
        auto tmp01 = tmp.partition(mh, mh, 0, mh);
        auto tmp10 = tmp.partition(mh, mh, mh, 0);
        auto tmp11 = tmp.partition(mh, mh, mh, mh);

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
