#ifndef FMMI_DMATRIX_HPP
#define FMMI_DMATRIX_HPP

#include <algorithm>
#include <iostream>

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
    template <bool other_managed>
    bool operator==(const dmatrix<T, other_managed>& other) const;

    template <bool other_managed>
    bool operator!=(const dmatrix<T, other_managed>& other) const;

    template <bool other_managed>
    bool equals(const dmatrix<T, other_managed>& other, double margin = 0.0) const;

    /// Copy over values from another matrix of the same size.
    /// Using the assignment operator (=) would be more elegant but it is
    /// implicitly removed by the compiler.
    template <bool other_managed>
    void assign(const dmatrix<T, other_managed>& other);

    void print() const;

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
template <bool other_managed>
bool dmatrix<T, managed>::operator==(const dmatrix<T, other_managed>& other) const
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
template <bool other_managed>
bool dmatrix<T, managed>::operator!=(const dmatrix<T, other_managed>& other) const
{
    return !(*this == other);
}


template <typename T, bool managed>
template <bool other_managed>
bool dmatrix<T, managed>::equals(const dmatrix<T, other_managed>& other, double margin) const
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
template <bool other_managed>
void dmatrix<T, managed>::assign(const dmatrix<T, other_managed>& other)
{
    if (width_ == stride_ && width_ == other.width_ && width_ == other.stride_)
    {
        std::copy(&other(0, 0), &other(height_, 0), &(*this)(0, 0));
    }
    else
    {
        for (uint16_t y = 0; y < height_; ++y)
        {
            std::copy(&other(y, 0), &other(y, width_), &(*this)(y, 0));
        }
    }
}


template <typename T, bool managed>
void dmatrix<T, managed>::print() const
{
    for (uint16_t y = 0; y < height_; ++y)
    {
        uint16_t x = 0;
        for (; x < width_ - 1; ++x)
        {
            std::cout << (*this)(y, x) << ", ";
        }
        std::cout << (*this)(y, x) << std::endl;
    }
}


template <typename T, bool managed>
dmatrix<T> dmatrix<T, managed>::identity(uint16_t height, uint16_t width)
{
    dmatrix<T> ident(height, width);
    for (uint16_t y = 0; y < height; ++y)
    {
        for (uint16_t x = 0; x < width; ++x)
        {
            ident(y, x) = y == x;
        }
    }
    return ident;
}


template <typename T, bool a_managed, bool b_managed, bool c_managed>
void add(const dmatrix<T, a_managed>& a,
         const dmatrix<T, b_managed>& b,
         dmatrix<T, c_managed>& c)
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


template <typename T, bool a_managed, bool b_managed, bool c_managed>
void sub(const dmatrix<T, a_managed>& a,
         const dmatrix<T, b_managed>& b,
         dmatrix<T, c_managed>& c)
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


template <typename T, bool a_managed, bool b_managed>
void add(const dmatrix<T, a_managed>& a, dmatrix<T, b_managed>& b)
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


template <typename T, bool a_managed, bool b_managed>
void sub(const dmatrix<T, a_managed>& a, dmatrix<T, b_managed>& b)
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


template <typename T, bool a_managed, bool b_managed, bool c_managed>
void mul(const dmatrix<T, a_managed>& a, const dmatrix<T, b_managed>& b, dmatrix<T, c_managed>& c)
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


template <typename T, bool a_managed, bool b_managed>
void transpose(const dmatrix<T, a_managed>& a, dmatrix<T, b_managed>& b)
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


template <typename T, bool a_managed, bool ainv_managed>
void inv(const dmatrix<T, a_managed>& a, dmatrix<T, ainv_managed>& ainv)
{
    assert(a.width() == a.height()
           && ainv.width() == ainv.height()
           && a.width() == ainv.width());

    const uint16_t m = a.width();

    dmatrix<T> atmp(a.height(), a.width());
    atmp.assign(a);
    ainv.assign(dmatrix<T>::identity(m, m));

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
