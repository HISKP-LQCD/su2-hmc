// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "../hybrid-monte-carlo.hpp"
#include "../sanity-checks.hpp"

#include <gtest/gtest.h>

#include <random>

TEST(hybridMonteCarlo, generate_from_gaussian) {
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    for (int i = 0; i < 10; ++i) {
        auto const mat = generate_from_gaussian(engine, dist);
        ASSERT_TRUE(is_hermitian(mat));
        auto const exponent = std::complex<double>{0, 1} * mat;
        ASSERT_TRUE(is_unitary(exponent.exp()));
    }
}
