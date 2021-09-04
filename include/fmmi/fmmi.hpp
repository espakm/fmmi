#ifndef FMMI_FMMI_HPP
#define FMMI_FMMI_HPP

#include "fmmi/matrix.hpp"

namespace fmmi
{

template <uint16_t m, uint16_t n, uint16_t p>
void mul_fast(const matrix<m, n>& mx1, const matrix<n, p>& mx2, matrix<m, p>& mx3);

template <>
void mul_fast(const matrix<2, 2>& mx1, const matrix<2, 2>& mx2, matrix<2, 2>& mx3);

}

#endif
