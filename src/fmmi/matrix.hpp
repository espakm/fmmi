#ifndef FMMI_MATRIX_HPP
#define FMMI_MATRIX_HPP

#include <algorithm>
#include <cstdint>

namespace fmmi
{

constexpr uint16_t log_2(uint16_t n)
{
    return n < 2 ? 0 : log_2(n / 2) + 1;
}


constexpr std::size_t exp_2(uint16_t n)
{
    return n == 0 ? 1 : 2 * exp_2(n - 1);
}


constexpr std::size_t padded_size(uint16_t n)
{
    return n < 2 ? 1 : exp_2(log_2(n - 1) + 1);
}


constexpr std::size_t padded_size(uint16_t height, uint16_t width)
{
    auto ps = padded_size(height > width ? height : width);
    return ps * ps;
}


template <typename T, uint16_t height, uint16_t width,
          uint16_t y0, uint16_t x0, uint16_t stride>
class matrix;


template <uint16_t height, uint16_t width,
          uint16_t y0 = 0, uint16_t x0 = 0, uint16_t stride = width>
using i16mx_t = matrix<int16_t, height, width, y0, x0, stride>;


template <uint16_t height, uint16_t width,
          uint16_t y0 = 0, uint16_t x0 = 0, uint16_t stride = width>
using i32mx_t = matrix<int32_t, height, width, y0, x0, stride>;


template <uint16_t height, uint16_t width,
          uint16_t y0 = 0, uint16_t x0 = 0, uint16_t stride = width>
using i64mx_t = matrix<int64_t, height, width, y0, x0, stride>;


template <uint16_t height, uint16_t width,
          uint16_t y0 = 0, uint16_t x0 = 0, uint16_t stride = width>
using f32mx_t = matrix<float, height, width, y0, x0, stride>;


template <uint16_t height, uint16_t width,
          uint16_t y0 = 0, uint16_t x0 = 0, uint16_t stride = width>
using f64mx_t = matrix<double, height, width, y0, x0, stride>;


template <typename T, uint16_t height, uint16_t width,
          uint16_t y0 = 0, uint16_t x0 = 0, uint16_t stride = width>
class matrix
{
public:
    matrix();

    matrix(std::initializer_list<T> init_list);

    inline
    const T& operator()(uint16_t y, uint16_t x) const;

    inline
    T& operator()(uint16_t y, uint16_t x);

    template <uint16_t p_height, uint16_t p_width, uint16_t p_y0 = 0, uint16_t p_x0 = 0>
    inline
    const matrix<T, p_height, p_width, y0 + p_y0, x0 + p_x0, stride>& partition() const;

    template <uint16_t p_height, uint16_t p_width, uint16_t p_y0 = 0, uint16_t p_x0 = 0>
    inline
    matrix<T, p_height, p_width, y0 + p_y0, x0 + p_x0, stride>& partition();

    /// Matrix equality.
    /// Note that this implementation does not allow any 'epsilon' difference
    /// between elements.
    template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
    bool operator==(const matrix<T, height, width, other_y0, other_x0, other_stride>& other) const;

    template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
    bool operator!=(const matrix<T, height, width, other_y0, other_x0, other_stride>& other) const;

    template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
    bool equals(const matrix<T, height, width, other_y0, other_x0, other_stride>& other, double margin = 0.0) const;

    static
    matrix<T, height, width> identity();

private:
    static constexpr
    std::size_t offset_ = y0 * stride + x0;

    T data_[height * width];
};


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
matrix<T, height, width, y0, x0, stride>::matrix()
//    : data_{}
{
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
matrix<T, height, width, y0, x0, stride>::matrix(std::initializer_list<T> init_list)
    : data_{}
{
    auto it_data = &data_[offset_];
    if (init_list.size() == 0)
    {
        for (uint16_t y = 0; y < height; ++y)
        {
            std::fill_n(it_data, width, 0);
            it_data += stride;
        }
    }
    else
    {
        auto it_init = std::begin(init_list);
        for (uint16_t y = 0; y < height; ++y)
        {
            std::copy_n(it_init, width, it_data);
            it_init += width;
            it_data += stride;
        }
    }
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
inline
const T& matrix<T, height, width, y0, x0, stride>::operator()(uint16_t y, uint16_t x) const
{
    return data_[offset_ + y * stride + x];
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
inline
T& matrix<T, height, width, y0, x0, stride>::operator()(uint16_t y, uint16_t x)
{
    return data_[offset_ + y * stride + x];
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
template <uint16_t p_height, uint16_t p_width, uint16_t p_y0, uint16_t p_x0>
inline
const matrix<T, p_height, p_width, y0 + p_y0, x0 + p_x0, stride>& matrix<T, height, width, y0, x0, stride>::partition() const
{
    return *reinterpret_cast<const matrix<T, p_height, p_width, y0 + p_y0, x0 + p_x0, stride>*>(this);
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
template <uint16_t p_height, uint16_t p_width, uint16_t p_y0, uint16_t p_x0>
inline
matrix<T, p_height, p_width, y0 + p_y0, x0 + p_x0, stride>& matrix<T, height, width, y0, x0, stride>::partition()
{
    return *reinterpret_cast<matrix<T, p_height, p_width, y0 + p_y0, x0 + p_x0, stride>*>(this);
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
bool matrix<T, height, width, y0, x0, stride>::operator==(const matrix<T, height, width, other_y0, other_x0, other_stride>& other) const
{
    if constexpr (width == stride && width == other_stride)
    {
        return std::equal(&(*this)(0, 0), &(*this)(height, 0), &other(0, 0));
    }

    for (uint16_t y = 0; y < height; ++y)
    {
        if (!std::equal(&(*this)(y, 0), &(*this)(y, width), &other(y, 0)))
        {
            return false;
        }
    }

    return true;
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
bool matrix<T, height, width, y0, x0, stride>::operator!=(const matrix<T, height, width, other_y0, other_x0, other_stride>& other) const
{
    return !(*this == other);
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
bool matrix<T, height, width, y0, x0, stride>::equals(
        const matrix<T, height, width, other_y0, other_x0, other_stride>& other,
        double margin) const
{
    for (uint16_t y = 0; y < height; ++y)
    {
        for (uint16_t x = 0; x < width; ++x)
        {
            if (std::abs((*this)(y, x) - other(y, x)) > margin)
            {
                return false;
            }
        }
    }

    return true;
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
matrix<T, height, width> matrix<T, height, width, y0, x0, stride>::identity()
{
    matrix<T, height, width> ident;
    for (uint16_t y = 0; y < height; ++y)
    {
        for (uint16_t x = 0; x < width; ++x)
        {
            ident(y, x) = y == x ? 1 : 0;
        }
    }
    return ident;
}


template <typename T, uint16_t m, uint16_t n,
          uint16_t y0_a, uint16_t x0_a, uint16_t stride_a,
          uint16_t y0_b, uint16_t x0_b, uint16_t stride_b,
          uint16_t y0_c, uint16_t x0_c, uint16_t stride_c>
void add(const matrix<T, m, n, y0_a, x0_a, stride_a>& a,
         const matrix<T, m, n, y0_b, x0_b, stride_b>& b,
         matrix<T, m, n, y0_c, x0_c, stride_c>& c)
{
    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < n; ++j)
        {
            c(i, j) = a(i, j) + b(i, j);
        }
    }
}


template <typename T, uint16_t m, uint16_t n,
          uint16_t y0_a, uint16_t x0_a, uint16_t stride_a,
          uint16_t y0_b, uint16_t x0_b, uint16_t stride_b,
          uint16_t y0_c, uint16_t x0_c, uint16_t stride_c>
void sub(const matrix<T, m, n, y0_a, x0_a, stride_a>& a,
         const matrix<T, m, n, y0_b, x0_b, stride_b>& b,
         matrix<T, m, n, y0_c, x0_c, stride_c>& c)
{
    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < n; ++j)
        {
            c(i, j) = a(i, j) - b(i, j);
        }
    }
}


template <typename T, uint16_t m, uint16_t n,
          uint16_t y0_a, uint16_t x0_a, uint16_t stride_a,
          uint16_t y0_b, uint16_t x0_b, uint16_t stride_b>
void add(const matrix<T, m, n, y0_a, x0_a, stride_a>& a,
         matrix<T, m, n, y0_b, x0_b, stride_b>& b)
{
    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < n; ++j)
        {
            b(i, j) = a(i, j) + b(i, j);
        }
    }
}


template <typename T, uint16_t m, uint16_t n,
          uint16_t y0_a, uint16_t x0_a, uint16_t stride_a,
          uint16_t y0_b, uint16_t x0_b, uint16_t stride_b>
void sub(const matrix<T, m, n, y0_a, x0_a, stride_a>& a,
         matrix<T, m, n, y0_b, x0_b, stride_b>& b)
{
    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < n; ++j)
        {
            b(i, j) = a(i, j) - b(i, j);
        }
    }
}


template <typename T, uint16_t m, uint16_t n, uint16_t p,
          uint16_t y0_a, uint16_t x0_a, uint16_t stride_a,
          uint16_t y0_b, uint16_t x0_b, uint16_t stride_b,
          uint16_t y0_c, uint16_t x0_c, uint16_t stride_c>
void mul(const matrix<T, m, n, y0_a, x0_a, stride_a>& a,
         const matrix<T, n, p, y0_b, x0_b, stride_b>& b,
         matrix<T, m, p, y0_c, x0_c, stride_c>& c)
{
    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < p; ++j)
        {
            T sum{};
            for (uint16_t k = 0; k < n; ++k)
            {
                sum += a(i, k) * b(k, j);
            }
            c(i, j) = sum;
        }
    }
}


template <typename T, uint16_t m, uint16_t n,
          uint16_t y0_a, uint16_t x0_a, uint16_t stride_a,
          uint16_t y0_b, uint16_t x0_b, uint16_t stride_b>
void transpose(const matrix<T, m, n, y0_a, x0_a, stride_a>& a,
         matrix<T, n, m, y0_b, x0_b, stride_b>& b)
{
    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < n; ++j)
        {
            b(j, i) = a(i, j);
        }
    }
}


template <typename T, uint16_t m,
          uint16_t y0_a, uint16_t x0_a, uint16_t stride_a,
          uint16_t y0_ainv, uint16_t x0_ainv, uint16_t stride_ainv>
inline
void inv(const matrix<T, m, m, y0_a, x0_a, stride_a>& a,
         matrix<T, m, m, y0_ainv, x0_ainv, stride_ainv>& ainv)
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

        matrix<T, m, m> tmp;
        auto& tmp00 = tmp.template partition<n, n, 0, 0>();
        auto& tmp01 = tmp.template partition<n, p, 0, n>();
        auto& tmp10 = tmp.template partition<p, n, n, 0>();
        auto& tmp11 = tmp.template partition<p, p, n, n>();

        inv(a00, tmp00);
        mul(a10, tmp00, tmp10);
        mul(tmp10, a01, tmp11);
        sub(a11, tmp11);
        inv(tmp11, ainv11);

        mul(tmp00, a01, tmp01);
        mul(tmp01, ainv11, ainv01);

        mul(ainv01, tmp10, ainv00);
        add(tmp00, ainv00);

        mul(ainv11, tmp10, ainv10);

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
