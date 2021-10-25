#ifndef FMMI_SMATRIX_HPP
#define FMMI_SMATRIX_HPP

#include <algorithm>
#include <cstdint>
#include <iostream>

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
class smatrix;


template <uint16_t height, uint16_t width,
          uint16_t y0 = 0, uint16_t x0 = 0, uint16_t stride = width>
using i16smx_t = smatrix<int16_t, height, width, y0, x0, stride>;


template <uint16_t height, uint16_t width,
          uint16_t y0 = 0, uint16_t x0 = 0, uint16_t stride = width>
using i32smx_t = smatrix<int32_t, height, width, y0, x0, stride>;


template <uint16_t height, uint16_t width,
          uint16_t y0 = 0, uint16_t x0 = 0, uint16_t stride = width>
using i64smx_t = smatrix<int64_t, height, width, y0, x0, stride>;


template <uint16_t height, uint16_t width,
          uint16_t y0 = 0, uint16_t x0 = 0, uint16_t stride = width>
using f32smx_t = smatrix<float, height, width, y0, x0, stride>;


template <uint16_t height, uint16_t width,
          uint16_t y0 = 0, uint16_t x0 = 0, uint16_t stride = width>
using f64smx_t = smatrix<double, height, width, y0, x0, stride>;


template <uint16_t height, uint16_t width,
          uint16_t y0 = 0, uint16_t x0 = 0, uint16_t stride = width>
using f128smx_t = smatrix<long double, height, width, y0, x0, stride>;


template <typename T, uint16_t height, uint16_t width,
          uint16_t y0 = 0, uint16_t x0 = 0, uint16_t stride = width>
class smatrix
{
public:
    smatrix();

    smatrix(std::initializer_list<T> init_list);

    template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
    smatrix(const smatrix<T, height, width, other_y0, other_x0, other_stride>& other);

    inline
    const T& operator()(uint16_t y, uint16_t x) const;

    inline
    T& operator()(uint16_t y, uint16_t x);

    template <uint16_t p_height, uint16_t p_width, uint16_t p_y0 = 0, uint16_t p_x0 = 0>
    inline
    const smatrix<T, p_height, p_width, y0 + p_y0, x0 + p_x0, stride>& partition() const;

    template <uint16_t p_height, uint16_t p_width, uint16_t p_y0 = 0, uint16_t p_x0 = 0>
    inline
    smatrix<T, p_height, p_width, y0 + p_y0, x0 + p_x0, stride>& partition();

    /// smatrix equality.
    /// Note that this implementation does not allow any 'epsilon' difference
    /// between elements.
    template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
    bool operator==(const smatrix<T, height, width, other_y0, other_x0, other_stride>& other) const;

    template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
    bool operator!=(const smatrix<T, height, width, other_y0, other_x0, other_stride>& other) const;

    template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
    bool equals(const smatrix<T, height, width, other_y0, other_x0, other_stride>& other, double margin = 0.0) const;

    template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
    smatrix<T, height, width, y0, x0, stride>& operator=(const smatrix<T, height, width, other_y0, other_x0, other_stride>& other);

    void print() const;

    static
    smatrix<T, height, width> identity();

private:
    static constexpr
    std::size_t offset_ = y0 * stride + x0;

    T data_[height * width];
};


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
smatrix<T, height, width, y0, x0, stride>::smatrix()
//    : data_{}
{
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
smatrix<T, height, width, y0, x0, stride>::smatrix(std::initializer_list<T> init_list)
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
template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
smatrix<T, height, width, y0, x0, stride>::smatrix(const smatrix<T, height, width, other_y0, other_x0, other_stride>& other)
{
    if constexpr (width == stride && width == other_stride)
    {
        std::copy(&other(0, 0), &other(height, 0), &(*this)(0, 0));
    }

    for (uint16_t y = 0; y < height; ++y)
    {
        std::copy(&other(y, 0), &other(y, width), &(*this)(y, 0));
    }
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
inline
const T& smatrix<T, height, width, y0, x0, stride>::operator()(uint16_t y, uint16_t x) const
{
    return data_[offset_ + y * stride + x];
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
inline
T& smatrix<T, height, width, y0, x0, stride>::operator()(uint16_t y, uint16_t x)
{
    return data_[offset_ + y * stride + x];
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
template <uint16_t p_height, uint16_t p_width, uint16_t p_y0, uint16_t p_x0>
inline
const smatrix<T, p_height, p_width, y0 + p_y0, x0 + p_x0, stride>& smatrix<T, height, width, y0, x0, stride>::partition() const
{
    return *reinterpret_cast<const smatrix<T, p_height, p_width, y0 + p_y0, x0 + p_x0, stride>*>(this);
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
template <uint16_t p_height, uint16_t p_width, uint16_t p_y0, uint16_t p_x0>
inline
smatrix<T, p_height, p_width, y0 + p_y0, x0 + p_x0, stride>& smatrix<T, height, width, y0, x0, stride>::partition()
{
    return *reinterpret_cast<smatrix<T, p_height, p_width, y0 + p_y0, x0 + p_x0, stride>*>(this);
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
bool smatrix<T, height, width, y0, x0, stride>::operator==(const smatrix<T, height, width, other_y0, other_x0, other_stride>& other) const
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
bool smatrix<T, height, width, y0, x0, stride>::operator!=(const smatrix<T, height, width, other_y0, other_x0, other_stride>& other) const
{
    return !(*this == other);
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
bool smatrix<T, height, width, y0, x0, stride>::equals(
        const smatrix<T, height, width, other_y0, other_x0, other_stride>& other,
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
template <uint16_t other_y0, uint16_t other_x0, uint16_t other_stride>
smatrix<T, height, width, y0, x0, stride>& smatrix<T, height, width, y0, x0, stride>::operator=(const smatrix<T, height, width, other_y0, other_x0, other_stride>& other)
{
    if constexpr (width == stride && width == other_stride)
    {
        std::copy(&other(0, 0), &other(height, 0), &(*this)(0, 0));
    }

    for (uint16_t y = 0; y < height; ++y)
    {
        std::copy(&other(y, 0), &other(y, width), &(*this)(y, 0));
    }

    return *this;
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
void smatrix<T, height, width, y0, x0, stride>::print() const
{
    for (uint16_t y = 0; y < height; ++y)
    {
        uint16_t x = 0;
        for (; x < width - 1; ++x)
        {
            std::cout << (*this)(y, x) << ", ";
        }
        std::cout << (*this)(y, x) << std::endl;
    }
}


template <typename T, uint16_t height, uint16_t width, uint16_t y0, uint16_t x0, uint16_t stride>
smatrix<T, height, width> smatrix<T, height, width, y0, x0, stride>::identity()
{
    smatrix<T, height, width> ident;
    for (uint16_t y = 0; y < height; ++y)
    {
        for (uint16_t x = 0; x < width; ++x)
        {
            ident(y, x) = y == x;
        }
    }
    return ident;
}


template <typename T, uint16_t m, uint16_t n,
          uint16_t a_y0, uint16_t a_x0, uint16_t a_stride,
          uint16_t b_y0, uint16_t b_x0, uint16_t b_stride,
          uint16_t c_y0, uint16_t c_x0, uint16_t c_stride>
void add(const smatrix<T, m, n, a_y0, a_x0, a_stride>& a,
         const smatrix<T, m, n, b_y0, b_x0, b_stride>& b,
         smatrix<T, m, n, c_y0, c_x0, c_stride>& c)
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
          uint16_t a_y0, uint16_t a_x0, uint16_t a_stride,
          uint16_t b_y0, uint16_t b_x0, uint16_t b_stride,
          uint16_t c_y0, uint16_t c_x0, uint16_t c_stride>
void sub(const smatrix<T, m, n, a_y0, a_x0, a_stride>& a,
         const smatrix<T, m, n, b_y0, b_x0, b_stride>& b,
         smatrix<T, m, n, c_y0, c_x0, c_stride>& c)
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
          uint16_t a_y0, uint16_t a_x0, uint16_t a_stride,
          uint16_t b_y0, uint16_t b_x0, uint16_t b_stride>
void add(const smatrix<T, m, n, a_y0, a_x0, a_stride>& a,
         smatrix<T, m, n, b_y0, b_x0, b_stride>& b)
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
          uint16_t a_y0, uint16_t a_x0, uint16_t a_stride,
          uint16_t b_y0, uint16_t b_x0, uint16_t b_stride>
void sub(const smatrix<T, m, n, a_y0, a_x0, a_stride>& a,
         smatrix<T, m, n, b_y0, b_x0, b_stride>& b)
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
          uint16_t a_y0, uint16_t a_x0, uint16_t a_stride,
          uint16_t b_y0, uint16_t b_x0, uint16_t b_stride,
          uint16_t c_y0, uint16_t c_x0, uint16_t c_stride>
void mul(const smatrix<T, m, n, a_y0, a_x0, a_stride>& a,
         const smatrix<T, n, p, b_y0, b_x0, b_stride>& b,
         smatrix<T, m, p, c_y0, c_x0, c_stride>& c)
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
          uint16_t a_y0, uint16_t a_x0, uint16_t a_stride,
          uint16_t b_y0, uint16_t b_x0, uint16_t b_stride>
void transpose(const smatrix<T, m, n, a_y0, a_x0, a_stride>& a,
         smatrix<T, n, m, b_y0, b_x0, b_stride>& b)
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
          uint16_t a_y0, uint16_t a_x0, uint16_t a_stride,
          uint16_t ainv_y0, uint16_t ainv_x0, uint16_t ainv_stride>
void inv(const smatrix<T, m, m, a_y0, a_x0, a_stride>& a,
         smatrix<T, m, m, ainv_y0, ainv_x0, ainv_stride>& ainv)
{
    smatrix<T, m, m> atmp = a;
    ainv = smatrix<T, m, m>::identity();

    for (uint16_t row = 0, lead = 0; row < m && lead < 2 * m; ++row, ++lead) {
        uint16_t i = row;
        while ((lead < m ? atmp(i, lead) : ainv(i, lead - m)) == 0)
        {
            if (++i == m)
            {
                i = row;
                if (++lead == 2 * m)
                {
                    return;
                }
            }
        }
        for (uint16_t column = 0; column < m; ++column)
        {
            std::swap(atmp(i, column), atmp(row, column));
            std::swap(ainv(i, column), ainv(row, column));
        }
        if ((lead < m ? atmp(row, lead) : ainv(row, lead - m)) != 0)
        {
            auto f = lead < m ? atmp(row, lead) : ainv(row, lead - m);
            for (uint16_t column = 0; column < m; ++column)
            {
                atmp(row, column) /= f;
                ainv(row, column) /= f;
            }
        }
        for (uint16_t j = 0; j < m; ++j)
        {
            if (j == row)
            {
                continue;
            }
            auto f = lead < m ? atmp(j, lead) : ainv(j, lead - m);
            for (uint16_t column = 0; column < m; ++column)
            {
                atmp(j, column) -= f * atmp(row, column);
                ainv(j, column) -= f * ainv(row, column);
            }
        }
    }
}

}

#endif
