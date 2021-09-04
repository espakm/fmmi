#ifndef FMMI_MATRIX_HPP
#define FMMI_MATRIX_HPP

#include <cstdint>

namespace fmmi
{

template <uint16_t height, uint16_t width>
class matrix
{
public:
    template <typename ...T>
    matrix(T... init_list)
        : data_{init_list...}
    {
    }

    inline
    double operator()(uint16_t rowIdx, uint16_t columnIdx) const;

    inline
    double& operator()(uint16_t rowIdx, uint16_t columnIdx);

private:
    double data_[width * height];
};


template <uint16_t height, uint16_t width>
inline
double matrix<height, width>::operator()(uint16_t rowIdx, uint16_t columnIdx) const
{
    return data_[rowIdx * width + columnIdx];
}


template <uint16_t height, uint16_t width>
inline
double& matrix<height, width>::operator()(uint16_t rowIdx, uint16_t columnIdx)
{
    return data_[rowIdx * width + columnIdx];
}


template <uint16_t m, uint16_t n>
void add(const matrix<m, n>& mx1, const matrix<m, n>& mx2, matrix<m, n>& mx3)
{
    for (uint16_t i = 0; i < m; ++i)
    {
        for (uint16_t j = 0; j < n; ++j)
        {
            mx3(i, j) = mx1(i, j) + mx2(i, j);
        }
    }
}


template <uint16_t m, uint16_t n, uint16_t p>
void mul(const matrix<m, n>& mx1, const matrix<n, p>& mx2, matrix<m, p>& mx3)
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
