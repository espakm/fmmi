#include "fmmi/fmmi.hpp"

namespace fmmi
{

int dummy;

//template <>
//void mul_fast(const matrix<2, 2>& a, const matrix<2, 2>& b, matrix<2, 2>& c)
//{
//    const auto p1 = (a(0, 0) + a(1, 1)) * (b(0, 0) + b(1, 1));
//    const auto p2 = (a(1, 0) + a(1, 1)) * b(0, 0);
//    const auto p3 = a(0, 0) * (b(0, 1) - b(1, 1));
//    const auto p4 = a(1, 1) * (b(1, 0) - b(0, 0));
//    const auto p5 = (a(0, 0) + a(0, 1)) * b(1, 1);
//    const auto p6 = (a(1, 0) - a(0, 0)) * (b(0, 0) + b(0, 1));
//    const auto p7 = (a(0, 1) - a(1, 1)) * (b(1, 0) + b(1, 1));

//    c(0, 0) = p1 + p4 - p5 + p7;
//    c(0, 1) = p3 + p5;
//    c(1, 0) = p2 + p4;
//    c(1, 1) = p1 + p3 - p2 + p6;
//}

}
