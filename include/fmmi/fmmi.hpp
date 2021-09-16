#ifndef FMMI_FMMI_HPP
#define FMMI_FMMI_HPP

#include "fmmi/matrix.hpp"

namespace fmmi
{

template <uint16_t m, uint16_t n, uint16_t p, uint16_t stride1, typename D1, uint16_t stride2, typename D2, uint16_t stride3, typename D3>
void mul_fast(const matrix<m, n, stride1, D1>& a, const matrix<n, p, stride2, D2>& b, matrix<m, p, stride3, D3>& c)
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

        return;
    }

    else if constexpr ((m % 2) == 0 && (n % 2) == 0 && (p % 2) == 0)
    {
        constexpr uint16_t mh = m / 2;
        constexpr uint16_t nh = n / 2;
        constexpr uint16_t ph = p / 2;

        const auto a00 = a.template partition<0, 0, mh, nh>();
        const auto a01 = a.template partition<0, nh, mh, nh>();
        const auto a10 = a.template partition<mh, 0, mh, nh>();
        const auto a11 = a.template partition<mh, nh, mh, nh>();


        const auto b00 = b.template partition<0, 0, nh, ph>();
        const auto b01 = b.template partition<0, ph, nh, ph>();
        const auto b10 = b.template partition<nh, 0, nh, ph>();
        const auto b11 = b.template partition<nh, ph, nh, ph>();

        auto c00 = c.template partition<0, 0, mh, ph>();
        auto c01 = c.template partition<0, ph, mh, ph>();
        auto c10 = c.template partition<mh, 0, mh, ph>();
        auto c11 = c.template partition<mh, ph, mh, ph>();

        matrix<mh, nh> tmp1;
        matrix<nh, ph> tmp2;

        matrix<mh, ph> p1;
        add(a00, a11, tmp1);
        add(b00, b11, tmp2);
        mul_fast(tmp1, tmp2, p1);

        matrix<mh, ph> p2;
        add(a10, a11, tmp1);
        mul_fast(tmp1, b00, p2);

        matrix<mh, ph> p3;
        sub(b01, b11, tmp2);
        mul_fast(a00, tmp2, p3);

        matrix<mh, ph> p4;
        sub(b10, b00, tmp2);
        mul_fast(a11, tmp2, p4);

        matrix<mh, ph> p5;
        add(a00, a01, tmp1);
        mul_fast(tmp1, b11, p5);

        matrix<mh, ph> p6;
        sub(a10, a00, tmp1);
        add(b00, b01, tmp2);
        mul_fast(tmp1, tmp2, p6);

        matrix<mh, ph> p7;
        sub(a01, a11, tmp1);
        add(b10, b11, tmp2);
        mul_fast(tmp1, tmp2, p7);

        add(p1, p4, c00);
        sub(c00, p5, c00);
        add(c00, p7, c00);

        add(p3, p5, c01);

        add(p2, p4, c10);

        add(p1, p3, c11);
        sub(c11, p2, c11);
        add(c11, p6, c11);

        return;
    }

    throw std::runtime_error("baj van");
}

}

#endif
