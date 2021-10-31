#include <catch.hpp>
#include <fmmi/common.hpp>

using namespace fmmi;

TEST_CASE("log_2")
{
    CHECK(log_2(1) == 0);
    CHECK(log_2(2) == 1);
    CHECK(log_2(3) == 1);
    CHECK(log_2(4) == 2);
    CHECK(log_2(5) == 2);
    CHECK(log_2(6) == 2);
    CHECK(log_2(7) == 2);
    CHECK(log_2(8) == 3);
    CHECK(log_2(9) == 3);
}


TEST_CASE("exp_2")
{
    CHECK(exp_2(0) == 1);
    CHECK(exp_2(1) == 2);
    CHECK(exp_2(2) == 4);
    CHECK(exp_2(3) == 8);
    CHECK(exp_2(4) == 16);
    CHECK(exp_2(5) == 32);
}


TEST_CASE("pow_2_upper_bound")
{
    CHECK(pow_2_upper_bound(1) == 1);
    CHECK(pow_2_upper_bound(2) == 2);
    CHECK(pow_2_upper_bound(3) == 4);
    CHECK(pow_2_upper_bound(4) == 4);
    CHECK(pow_2_upper_bound(5) == 8);
    CHECK(pow_2_upper_bound(6) == 8);
    CHECK(pow_2_upper_bound(7) == 8);
    CHECK(pow_2_upper_bound(8) == 8);
    CHECK(pow_2_upper_bound(9) == 16);

    CHECK(pow_2_upper_bound(3, 4) == 4 * 4);
    CHECK(pow_2_upper_bound(3, 5) == 8 * 8);
    CHECK(pow_2_upper_bound(5, 3) == 8 * 8);
}
