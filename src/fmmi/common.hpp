#ifndef FMMI_COMMON_HPP
#define FMMI_COMMON_HPP

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


constexpr std::size_t pow_2_upper_bound(uint16_t n)
{
    return n < 2 ? 1 : exp_2(log_2(n - 1) + 1);
}


constexpr std::size_t pow_2_upper_bound(uint16_t height, uint16_t width)
{
    auto pow_2 = pow_2_upper_bound(height > width ? height : width);
    return pow_2 * pow_2;
}

}

#endif
