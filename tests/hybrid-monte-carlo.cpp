// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "../hybrid-monte-carlo.hpp"
#include "../sanity-checks.hpp"

#include <gtest/gtest.h>

#include <random>

TEST(hybridMonteCarlo, generateFromAlgebra) {
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    for (int i = 0; i < 10; ++i) {
        auto const mat = random_from_algebra(engine, dist);
        ASSERT_TRUE(is_hermitian(mat)) << "Happened at i = " << i << "\n" << mat;
        ASSERT_TRUE(is_traceless(mat)) << "Happened at i = " << i << "\n" << mat;
    }
}

TEST(hybridMonteCarlo, generateFromGroup) {
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    for (int i = 0; i < 10; ++i) {
        auto const mat = random_from_group(engine, dist);
        ASSERT_TRUE(is_unitary(mat)) << "Happened at i = " << i << "\n" << mat;
    }
}

TEST(hybridMonteCarlo, randomizeAlgebra) {
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    Configuration config(10, 10);
    randomize_algebra(config, engine, dist);
    for (int i = 0; i < config.get_size(); ++i) {
        ASSERT_TRUE(is_hermitian(config[i])) << "Happened at i = " << i << "\n"
                                             << config[i];
        ASSERT_TRUE(is_traceless(config[i])) << "Happened at i = " << i << "\n"
                                             << config[i];
    }
}

TEST(hybridMonteCarlo, randomizeGroup) {
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    Configuration config(10, 10);
    randomize_algebra(config, engine, dist);
    for (int i = 0; i < config.get_size(); ++i) {
        ASSERT_TRUE(is_unitary(config[i])) << "Happened at i = " << i << "\n"
                                           << config[i];
    }
}
