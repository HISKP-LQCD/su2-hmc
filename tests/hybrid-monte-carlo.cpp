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

/**
  Tests whether a plaquette built from four matrices is gauge invariant.
  */
TEST(plaquette, simpleInvariance) {
    // Random generator setup.
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    // Generate four random SU(2) matrices.
    Matrix const m1 = random_from_group(engine, dist);
    Matrix const m2 = random_from_group(engine, dist);
    Matrix const m3 = random_from_group(engine, dist);
    Matrix const m4 = random_from_group(engine, dist);

    // Generate a random SU(2) matrix to be used as a transformation.
    Matrix const transformation = random_from_group(engine, dist);
    ASSERT_TRUE(is_unitary(transformation));

    // Compute the plaquette with the original matrices.
    Matrix const old_plaquette = m1 * m2 * m3.adjoint() * m4.adjoint();

    // Transform the four matrices.
    Matrix const m1t = transformation * m1;
    Matrix const m2t = transformation * m2;
    Matrix const m3t = transformation * m3;
    Matrix const m4t = transformation * m4;

    // Compute the new plaquette. This should not change!
    Matrix const new_plaquette = m1t * m2t * m3t.adjoint() * m4t.adjoint();

    // Compare old and new plaquette.
    ASSERT_TRUE(is_equal(old_plaquette, new_plaquette))
        << "old plaquette:\n"
        << old_plaquette << "\n"
        << "new plaquette:\n"
        << new_plaquette << "\n"
        << "transformation:\n"
        << transformation << "\n"
        << "transformation * transformation^\\dagger:\n"
        << (transformation * transformation.adjoint().eval()) << "\n"
        << "m1:\n"
        << m1 << "\n"
        << "m2:\n"
        << m2 << "\n"
        << "m3:\n"
        << m3 << "\n"
        << "m4:\n"
        << m4 << "\n"
        << "m1t:\n"
        << m1t << "\n"
        << "m2t:\n"
        << m2t << "\n"
        << "m3t:\n"
        << m3t << "\n"
        << "m4t:\n"
        << m4t << "\n";
}

TEST(plaquette, globalGaugeInvarianceReversibility) {
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


TEST(plaquette, globalGaugeInvariancePlaquetteInvariance) {
    std::mt19937 engine(0);
    std::normal_distribution<double> dist(0, 1);

    Configuration links = make_hot_start(10, 10, 1, 0);
    Configuration const old_links = links;

    Matrix const transformation = random_from_group(engine, dist);
    global_gauge_transformation(transformation, links);

    for (int n1 = 0; n1 < links.length_time; ++n1) {
        for (int n2 = 0; n2 < links.length_space; ++n2) {
            for (int n3 = 0; n3 < links.length_space; ++n3) {
                for (int n4 = 0; n4 < links.length_space; ++n4) {
                    for (int mu = 0; mu < 4; ++mu) {
                        for (int nu = 0; nu < 4; ++nu) {
                            Matrix const before =
                                get_plaquette(n1, n2, n3, n4, mu, nu, old_links);
                            Matrix const after =
                                get_plaquette(n1, n2, n3, n4, mu, nu, links);

                            ASSERT_TRUE(is_equal(before, after))
                                << "Before:\n"
                                << before << "\n"
                                << "After:\n"
                                << after << "\n"
                                << "At: t=" << n1 << ", x=" << n2 << ", y=" << n3 << ", z=" << n4
                                << "; mu=" << mu << ", nu=" << nu;
                        }
                    }
                }
            }
        }
    }
}

TEST(plaquette, globalGaugeInvarianceAverages) {
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

    // Compare values that are transformed forth and back with original values.
    ASSERT_NEAR(old_plaquette.real(), old2_plaquette.real(), error);
    ASSERT_NEAR(old_plaquette.imag(), old2_plaquette.imag(), error);
    ASSERT_NEAR(old_energy, old2_energy, error);

    // The plaquettes should be qauge invariant. Therefore those should be the same as
    // well.
    ASSERT_NEAR(old_plaquette.real(), new_plaquette.real(), error);
    ASSERT_NEAR(old_plaquette.imag(), new_plaquette.imag(), error);
    ASSERT_NEAR(old_energy, new_energy, error);
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
        << m4;
}
