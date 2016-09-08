// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "../sanity-checks.hpp"

#include <Eigen/Dense>
#include <gtest/gtest.h>

using Matrix = Eigen::Matrix2cd;
using Complex = std::complex<double>;

TEST(sanityChecks, isZero) {
    Matrix m;
    m << 0, 0, 0, 0;
    ASSERT_TRUE(is_zero(m));

    m << 1, 0, 1, 0;
    ASSERT_FALSE(is_zero(m));
}

TEST(sanityChecks, isEqual) {
    Matrix m1;
    m1 << 0, 1, 0, 2;
    Matrix m2;
    m2 << 0, 1, 3, 2;

    ASSERT_TRUE(is_equal(m1, m1));
    ASSERT_FALSE(is_equal(m1, m2));
}

TEST(sanityChecks, isHermitian) {
    Matrix m;
    m << Complex{0, 0}, Complex{0, -1}, Complex{0, 1}, Complex{0, 0};
    ASSERT_TRUE(is_hermitian(m));

    m << 1, 0, 1, 0;
    ASSERT_FALSE(is_hermitian(m));
}
