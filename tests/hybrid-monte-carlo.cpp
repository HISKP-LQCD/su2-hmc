// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "../hybrid-monte-carlo.hpp"
#include "../sanity-checks.hpp"

#include <gtest/gtest.h>

#include <random>

TEST(hybridMonteCarlo, generateFromAlgebra) {
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    for (int i = 0; i < 10; ++i) {
        Matrix const mat = random_from_algebra(engine, dist);
        ASSERT_TRUE(is_hermitian(mat)) << "Happened at i = " << i << "\n" << mat;
        ASSERT_TRUE(is_traceless(mat)) << "Happened at i = " << i << "\n" << mat;
    }
}

TEST(hybridMonteCarlo, generateFromGroup) {
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    for (int i = 0; i < 10; ++i) {
        Matrix const mat = random_from_group(engine, dist);
        ASSERT_TRUE(is_unitary(mat)) << "Happened at i = " << i << "\n" << mat;
        ASSERT_TRUE(is_unit_determinant(mat)) << "Happened at i = " << i << "\n" << mat;
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
    Matrix const r2 = - group_from_algebra(m);
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
        ASSERT_TRUE(is_unit_determinant(config[i]))
            << "Happened at i = " << i << "\n"
            << config[i] << "\nU U^\\dagger:\n"
            << (config[i] * config[i].adjoint().eval());
    }
}

TEST(hybridMonteCarlo, globalGaugeInvariance) {
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    Configuration links = make_hot_start(10, 10, 1, 0);

    auto const old_plaquette = get_average_plaquette(links);
    auto const old_energy = get_link_energy(links);

    Matrix const transformation = random_from_group(engine, dist);
    global_gauge_transformation(transformation, links);

    auto const new_plaquette = get_average_plaquette(links);
    auto const new_energy = get_link_energy(links);

    global_gauge_transformation(transformation.adjoint(), links);

    auto const old2_plaquette = get_average_plaquette(links);
    auto const old2_energy = get_link_energy(links);

    auto const error = 1e-10;

    ASSERT_NEAR(old_plaquette.real(), old2_plaquette.real(), error);
    ASSERT_NEAR(old_plaquette.imag(), old2_plaquette.imag(), error);
    ASSERT_NEAR(old_energy, old2_energy, error);

    ASSERT_NEAR(old_plaquette.real(), new_plaquette.real(), error);
    ASSERT_NEAR(old_plaquette.imag(), new_plaquette.imag(), error);
    ASSERT_NEAR(old_energy, new_energy, error);
}

TEST(hybridMonteCarlo, globalGaugeTransformation) {
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    auto links = make_hot_start(10, 10, 1, 0);

    Matrix const transformation = random_from_group(engine, dist);

    auto const old_links = links;

    for (int i = 0; i < links.get_size(); ++i) {
        ASSERT_TRUE(is_unitary(links[i]));
        ASSERT_TRUE(is_unit_determinant(links[i]));
    }
    global_gauge_transformation(transformation, links);
    for (int i = 0; i < links.get_size(); ++i) {
        ASSERT_TRUE(is_unitary(links[i]));
        ASSERT_TRUE(is_unit_determinant(links[i]));
    }
    global_gauge_transformation(transformation.adjoint(), links);
    for (int i = 0; i < links.get_size(); ++i) {
        ASSERT_TRUE(is_unitary(links[i]));
        ASSERT_TRUE(is_unit_determinant(links[i]));
        ASSERT_TRUE(is_equal(old_links[i], links[i]));
    }

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

TEST(hybridMonteCarlo, singlePlaquette) {
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    // Generate an empty configuration of links.
    Configuration links(10, 10);

    Matrix const m1 = random_from_group(engine, dist);
    Matrix const m2 = random_from_group(engine, dist);
    Matrix const m3 = random_from_group(engine, dist);
    Matrix const m4 = random_from_group(engine, dist);


    links(0, 0, 0, 0, 0) = m1;
    links(1, 0, 0, 0, 1) = m2;
    links(0, 1, 0, 0, 0) = m3;
    links(0, 0, 0, 0, 1) = m4;

    Matrix const expected_plaquette = m1 * m2 * m3.adjoint() * m4.adjoint();
    Matrix const actual_plaquette = get_plaquette(0, 0, 0, 0, 0, 1, links, false);

    ASSERT_TRUE(is_equal(expected_plaquette, actual_plaquette))
        << "Expected plaquette:\n"
        << expected_plaquette << "\nActual plaquette:\n"
        << actual_plaquette << "\n"
        << "m1:\n"
        << m1 << "\n"
        << "m2:\n"
        << m2 << "\n"
        << "m3:\n"
        << m3 << "\n"
        << "m4:\n"
        << m4 << "\n";
}
