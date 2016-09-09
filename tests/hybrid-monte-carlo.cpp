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

TEST(hybridMonteCarlo, groupFromAlgebra) {
    auto const pi = 4 * atan(1);
    ASSERT_LT(3.141, pi);
    ASSERT_LT(pi, 3.142);

    Matrix m;
    m << 0, 0, 0, 0;
    ASSERT_TRUE(is_unity(group_from_algebra(m)));

    m << 0, pi, pi, 0;
    auto const r2 = - group_from_algebra(m);
    ASSERT_TRUE(is_unity(r2)) << r2;
}

TEST(hybridMonteCarlo, randomizeGroup) {
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    Configuration config(10, 10);
    randomize_group(config, engine, dist);
    for (int i = 0; i < config.get_size(); ++i) {
        ASSERT_TRUE(is_unitary(config[i])) << "Happened at i = " << i << "\n"
                                           << config[i] << "\nU U^\\dagger:\n"
                                           << (config[i] * config[i].adjoint().eval());
    }
}

TEST(hybridMonteCarlo, globalGaugeInvariance) {
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    auto links = make_hot_start(10, 10, 1, 0);

    auto const old_plaquette = get_average_plaquette(links);
    auto const old_energy = get_link_energy(links);

    auto const transformation = random_from_group(engine, dist);
    global_gauge_transformation(transformation, links);

    auto const new_plaquette = get_average_plaquette(links);
    auto const new_energy = get_link_energy(links);

    ASSERT_DOUBLE_EQ(old_plaquette.real(), new_plaquette.real());
    ASSERT_DOUBLE_EQ(old_plaquette.imag(), new_plaquette.imag());
    ASSERT_DOUBLE_EQ(old_energy, new_energy);
}

TEST(hybridMonteCarlo, coldStartAveragePlaquette) {
    Configuration links(10, 10);
    for (int i = 0; i < links.get_size(); ++i) {
        links[i] = Matrix::Identity();
    }
    auto const average_plaquette = get_average_plaquette(links);
    ASSERT_DOUBLE_EQ(1.0, average_plaquette.real());
    ASSERT_DOUBLE_EQ(0.0, average_plaquette.imag());
}
