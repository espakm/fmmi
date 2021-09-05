#ifndef FMMI_MATRIX_HPP
#define FMMI_MATRIX_HPP

#include <algorithm>
#include <cstdint>

namespace fmmi
{

template <uint16_t height, uint16_t width, uint16_t stride = width, typename D = double[height * width]>
struct matrix
{
    template <typename ...T>
    matrix(T... init_list)
        : data_{init_list...}
    {
    }

    inline
    double operator()(uint16_t y, uint16_t x) const;

    inline
    double& operator()(uint16_t y, uint16_t x);

    /// Matrix equality.
    /// Note that this implementation does not allow any 'epsilon' difference
    /// between elements.
    template <uint16_t stride2, typename D2>
    inline
    bool operator==(const matrix<height, width, stride2, D2>& mx) const;

    template <uint16_t stride2, typename D2>
    inline
    bool operator!=(const matrix<height, width, stride2, D2>& mx) const;

    template <uint16_t y, uint16_t x, uint16_t p_height, uint16_t p_width>
    inline
    matrix<p_height, p_width, stride, double*> partition();

    template <uint16_t y, uint16_t x, uint16_t p_height, uint16_t p_width>
    inline
    const matrix<p_height, p_width, stride, double*> partition() const;

    matrix(double* data)
        : data_(data)
    {
    }

    D data_;
};


template <uint16_t height, uint16_t width, uint16_t stride, typename D>
inline
double matrix<height, width, stride, D>::operator()(uint16_t y, uint16_t x) const
{
    return data_[y * stride + x];
}


template <uint16_t height, uint16_t width, uint16_t stride, typename D>
inline
double& matrix<height, width, stride, D>::operator()(uint16_t y, uint16_t x)
{
    return data_[y * stride + x];
}


template <uint16_t height, uint16_t width, uint16_t stride, typename D>
template <uint16_t stride2, typename D2>
inline
bool matrix<height, width, stride, D>::operator==(const matrix<height, width, stride2, D2>& mx) const
{
    if constexpr (width == stride && width == stride2)
    {
        return std::equal(&data_[0], &data_[height * width], &mx.data_[0]);
    }

    for (uint16_t y = 0; y < height; ++y)
    {
        if (!std::equal(&data_[y * stride], &data_[y * stride + width], &mx.data_[y * stride2]))
        {
            return false;
        }
    }
    return true;
}


template <uint16_t height, uint16_t width, uint16_t stride, typename D>
template <uint16_t stride2, typename D2>
inline
bool matrix<height, width, stride, D>::operator!=(const matrix<height, width, stride2, D2>& mx) const
{
    return !(*this == mx);
}


template <uint16_t height, uint16_t width, uint16_t stride, typename D>
template <uint16_t y, uint16_t x, uint16_t p_height, uint16_t p_width>
inline
matrix<p_height, p_width, stride, double*> matrix<height, width, stride, D>::partition()
{
    return matrix<p_height, p_width, stride, double*>(&data_[y * stride + x]);
}


template <uint16_t height, uint16_t width, uint16_t stride, typename D>
template <uint16_t y, uint16_t x, uint16_t p_height, uint16_t p_width>
inline
const matrix<p_height, p_width, stride, double*> matrix<height, width, stride, D>::partition() const
{
    return matrix<p_height, p_width, stride, double*>(&data_[y * stride + x]);
}


template <uint16_t m, uint16_t n, uint16_t stride1, typename D1, uint16_t stride2, typename D2, uint16_t stride3, typename D3>
void add(const matrix<m, n, stride1, D1>& mx1, const matrix<m, n, stride2, D2>& mx2, matrix<m, n, stride3, D3>& mx3)
{
    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < n; ++j)
        {
            mx3(i, j) = mx1(i, j) + mx2(i, j);
        }
    }
}


template <uint16_t m, uint16_t n, uint16_t p, uint16_t stride1, typename D1, uint16_t stride2, typename D2, uint16_t stride3, typename D3>
void mul(const matrix<m, n, stride1, D1>& mx1, const matrix<n, p, stride2, D2>& mx2, matrix<m, p, stride3, D3>& mx3)
{
    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < p; ++j)
        {
            double sum = 0.0;
            for (uint16_t k = 0; k < n; ++k)
            {
                sum += mx1(i, k) * mx2(k, j);
            }
            mx3(i, j) = sum;
        }
    }
}

}

#endif
