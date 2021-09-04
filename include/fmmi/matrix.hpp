#ifndef FMMI_MATRIX_HPP
#define FMMI_MATRIX_HPP

#include <cstdint>


template <uint16_t width, uint16_t height>
class matrix
{
public:
    template <typename ...T>
    matrix(T... init_list)
        : data_{init_list...}
    {
    }

    inline
    double operator()(uint16_t columnIdx, uint16_t rowIdx) const;

    inline
    double& operator()(uint16_t columnIdx, uint16_t rowIdx);

private:
    double data_[width * height];
};


template <uint16_t width, uint16_t height>
inline
double matrix<width, height>::operator()(uint16_t columnIdx, uint16_t rowIdx) const
{
    return data_[rowIdx * width + columnIdx];
}


template <uint16_t width, uint16_t height>
inline
double& matrix<width, height>::operator()(uint16_t columnIdx, uint16_t rowIdx)
{
    return data_[rowIdx * width + columnIdx];
}

#endif
