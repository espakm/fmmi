#ifndef FMMI_DMATRIX_HPP
#define FMMI_DMATRIX_HPP

#include <algorithm>

namespace fmmi
{

template <typename T, bool managed>
class dmatrix;


using i16dmx_t = dmatrix<int16_t, true>;
using i32dmx_t = dmatrix<int32_t, true>;
using i64dmx_t = dmatrix<int64_t, true>;
using f32dmx_t = dmatrix<float, true>;
using f64dmx_t = dmatrix<double, true>;
using f128dmx_t = dmatrix<long double, true>;


template <typename T, bool managed = true>
class dmatrix
{
public:
    dmatrix(uint16_t height, uint16_t width);

    dmatrix(uint16_t height, uint16_t width, std::initializer_list<T> init_list);

    ~dmatrix();

    inline
    uint16_t height() const;

    inline
    uint16_t width() const;

    inline
    const T& operator()(uint16_t y, uint16_t x) const;

    inline
    T& operator()(uint16_t y, uint16_t x);

    inline
    const dmatrix<T, false> partition(uint16_t p_height, uint16_t p_width, uint16_t p_y0 = 0, uint16_t p_x0 = 0) const;

    inline
    dmatrix<T, false> partition(uint16_t p_height, uint16_t p_width, uint16_t p_y0 = 0, uint16_t p_x0 = 0);

    /// Matrix equality.
    /// Note that this implementation does not allow any 'epsilon' difference
    /// between elements.
    template <bool managed_other>
    bool operator==(const dmatrix<T, managed_other>& other) const;

    template <bool managed_other>
    bool operator!=(const dmatrix<T, managed_other>& other) const;

    template <bool managed_other>
    bool equals(const dmatrix<T, managed_other>& other, double margin = 0.0) const;

    static
    dmatrix<T> identity(uint16_t height, uint16_t width);

private:
    dmatrix(const T* data, uint16_t height, uint16_t width, uint16_t stride);
    dmatrix(T* data, uint16_t height, uint16_t width, uint16_t stride);

    T* const data_;

    const uint16_t height_;
    const uint16_t width_;
    const uint16_t stride_;

    template <typename, bool>
    friend class dmatrix;
};


template <typename T, bool managed>
dmatrix<T, managed>::dmatrix(uint16_t height, uint16_t width)
    : data_(new T[height * width]),
      height_(height),
      width_(width),
      stride_(width)
{
}


template <typename T, bool managed>
dmatrix<T, managed>::dmatrix(uint16_t height, uint16_t width, std::initializer_list<T> init_list)
    : dmatrix<T>(height, width)
{
    auto it_init = std::begin(init_list);
    auto it_data = data_;
    std::copy_n(it_init, height * width, it_data);
}


template <typename T, bool managed>
dmatrix<T, managed>::~dmatrix()
{
    if (managed)
    {
        delete[] data_;
    }
}


template <typename T, bool managed>
dmatrix<T, managed>::dmatrix(const T* data, uint16_t height, uint16_t width, uint16_t stride)
    : dmatrix<T, managed>(const_cast<T*>(data), height, width, stride)
{
}


template <typename T, bool managed>
dmatrix<T, managed>::dmatrix(T* data, uint16_t height, uint16_t width, uint16_t stride)
    : data_(data),
      height_(height),
      width_(width),
      stride_(stride)
{
}


template <typename T, bool managed>
inline
uint16_t dmatrix<T, managed>::height() const
{
    return height_;
}


template <typename T, bool managed>
inline
uint16_t dmatrix<T, managed>::width() const
{
    return width_;
}


template <typename T, bool managed>
inline
const T& dmatrix<T, managed>::operator()(uint16_t y, uint16_t x) const
{
    return data_[y * stride_ + x];
}


template <typename T, bool managed>
inline
T& dmatrix<T, managed>::operator()(uint16_t y, uint16_t x)
{
    return data_[y * stride_ + x];
}


template <typename T, bool managed>
inline
const dmatrix<T, false> dmatrix<T, managed>::partition(uint16_t height, uint16_t width, uint16_t y0, uint16_t x0) const
{
    return dmatrix<T, false>(&(*this)(y0, x0), height, width, stride_);
}


template <typename T, bool managed>
inline
dmatrix<T, false> dmatrix<T, managed>::partition(uint16_t height, uint16_t width, uint16_t y0, uint16_t x0)
{
    return dmatrix<T, false>(&(*this)(y0, x0), height, width, stride_);
}


template <typename T, bool managed>
template <bool managed_other>
bool dmatrix<T, managed>::operator==(const dmatrix<T, managed_other>& other) const
{
    if (width_ == stride_ && width_ == other.width_ && width_ == other.stride_)
    {
        return std::equal(&(*this)(0, 0), &(*this)(height_, 0), &other(0, 0));
    }

    for (uint16_t y = 0; y < height_; ++y)
    {
        if (!std::equal(&(*this)(y, 0), &(*this)(y, width_), &other(y, 0)))
        {
            return false;
        }
    }

    return true;
}


template <typename T, bool managed>
template <bool managed_other>
bool dmatrix<T, managed>::operator!=(const dmatrix<T, managed_other>& other) const
{
    return !(*this == other);
}


template <typename T, bool managed>
template <bool managed_other>
bool dmatrix<T, managed>::equals(const dmatrix<T, managed_other>& other, double margin) const
{
    for (uint16_t y = 0; y < height_; ++y)
    {
        for (uint16_t x = 0; x < width_; ++x)
        {
            if (std::abs((*this)(y, x) - other(y, x)) > margin)
            {
                return false;
            }
        }
    }

    return true;
}


template <typename T, bool managed>
dmatrix<T> dmatrix<T, managed>::identity(uint16_t height, uint16_t width)
{
    dmatrix<T> ident(height, width);
    for (uint16_t y = 0; y < height; ++y)
    {
        for (uint16_t x = 0; x < width; ++x)
        {
            ident(y, x) = y == x ? 1 : 0;
        }
    }
    return ident;
}


template <typename T, bool managed_a, bool managed_b, bool managed_c>
void add(const dmatrix<T, managed_a>& a,
         const dmatrix<T, managed_b>& b,
         dmatrix<T, managed_c>& c)
{
    assert(a.height() == b.height() && a.width() == b.width());
    assert(a.height() == c.height() && a.width() == c.width());

    const uint16_t m = a.height();
    const uint16_t n = a.width();

    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < n; ++j)
        {
            c(i, j) = a(i, j) + b(i, j);
        }
    }
}


template <typename T, bool managed_a, bool managed_b, bool managed_c>
void sub(const dmatrix<T, managed_a>& a,
         const dmatrix<T, managed_b>& b,
         dmatrix<T, managed_c>& c)
{
    assert(a.height() == b.height() && a.width() == b.width());
    assert(a.height() == c.height() && a.width() == c.width());

    const uint16_t m = a.height();
    const uint16_t n = a.width();

    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < n; ++j)
        {
            c(i, j) = a(i, j) - b(i, j);
        }
    }
}


template <typename T, bool managed_a, bool managed_b>
void add(const dmatrix<T, managed_a>& a, dmatrix<T, managed_b>& b)
{
    assert(a.height() == b.height() && a.width() == b.width());

    const uint16_t m = a.height();
    const uint16_t n = a.width();

    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < n; ++j)
        {
            b(i, j) = a(i, j) + b(i, j);
        }
    }
}


template <typename T, bool managed_a, bool managed_b>
void sub(const dmatrix<T, managed_a>& a, dmatrix<T, managed_b>& b)
{
    assert(a.height() == b.height() && a.width() == b.width());

    const uint16_t m = a.height();
    const uint16_t n = a.width();

    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < n; ++j)
        {
            b(i, j) = a(i, j) - b(i, j);
        }
    }
}


template <typename T, bool managed_a, bool managed_b, bool managed_c>
void mul(const dmatrix<T, managed_a>& a, const dmatrix<T, managed_b>& b, dmatrix<T, managed_c>& c)
{
    assert(a.width() == b.height()
           && a.height() == c.height()
           && b.width() == c.width());

    const uint16_t m = a.height();
    const uint16_t n = a.width();
    const uint16_t p = b.width();

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


template <typename T, bool managed_a, bool managed_b>
void transpose(const dmatrix<T, managed_a>& a, dmatrix<T, managed_b>& b)
{
    assert(a.height() == b.width()
           && a.width() == b.height());

    const uint16_t m = a.height();
    const uint16_t n = a.width();

    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < n; ++j)
        {
            b(j, i) = a(i, j);
        }
    }
}


template <typename T, bool managed_a, bool managed_ainv>
void inv(const dmatrix<T, managed_a>& a, dmatrix<T, managed_ainv>& ainv)
{
    assert(a.width() == a.height()
           && ainv.width() == ainv.height()
           && a.width() == ainv.width());

    const uint16_t m = a.width();

    if (m == 1)
    {
        ainv(0, 0) = 1 / a(0, 0);
    }
    else if (m == 2)
    {
        const T coeff = 1.0 / (a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0));
        ainv(0, 0) = a(1, 1) * coeff;
        ainv(0, 1) = -a(0, 1) * coeff;
        ainv(1, 0) = -a(1, 0) * coeff;
        ainv(1, 1) = a(0, 0) * coeff;
    }
    else
    {
        const uint16_t n = (m % 2) == 0 ? m / 2 : 1;
        const uint16_t p = m - n;

        const auto& a00 = a.partition(n, n, 0, 0);
        const auto& a01 = a.partition(n, p, 0, n);
        const auto& a10 = a.partition(p, n, n, 0);
        const auto& a11 = a.partition(p, p, n, n);

        auto ainv00 = ainv.partition(n, n, 0, 0);
        auto ainv01 = ainv.partition(n, p, 0, n);
        auto ainv10 = ainv.partition(p, n, n, 0);
        auto ainv11 = ainv.partition(p, p, n, n);

        dmatrix<T> tmp(m, m);
        auto tmp00 = tmp.partition(n, n, 0, 0);
        auto tmp01 = tmp.partition(n, p, 0, n);
        auto tmp10 = tmp.partition(p, n, n, 0);
        auto tmp11 = tmp.partition(p, p, n, n);

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
