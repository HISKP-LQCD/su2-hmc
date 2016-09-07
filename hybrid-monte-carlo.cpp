// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "hybrid-monte-carlo.hpp"

#include "pauli-matrices.hpp"

#include <unsupported/Eigen/MatrixFunctions>

#include <iostream>
#include <random>

Eigen::Matrix2cd generate_from_gaussian(std::mt19937 &engine,
                                 std::normal_distribution<double> &dist) {
    std::vector<double> coefficients = {dist(engine), dist(engine),
                                        dist(engine)};
    auto const pauli_matrices = PauliMatrices::get_instance();
    Eigen::Matrix2cd algebra_element;
    for (int i = 0; i != 3; ++i) {
        Eigen::Matrix2cd scaled_generator =
            coefficients[i] * pauli_matrices.get(i);
        algebra_element += scaled_generator;
    }

    return algebra_element;
}

Configuration make_hot_start(int const length_space,
                             int const length_time,
                             double const std,
                             int const seed) {
    std::mt19937 engine(seed);
    std::normal_distribution<double> dist(0, std);

    Configuration links(length_space, length_time);
    randomize(links, engine, dist);
    return links;
}

void randomize(Configuration &links,
               std::mt19937 &engine,
               std::normal_distribution<double> &dist) {
    for (int n1 = 0; n1 != links.length_time; ++n1) {
        std::cout << "n1 = " << n1 << std::endl;
        for (int n2 = 0; n2 != links.length_space; ++n2) {
            for (int n3 = 0; n3 != links.length_space; ++n3) {
                for (int n4 = 0; n4 != links.length_space; ++n4) {
                    for (int mu = 0; mu != 4; ++mu) {
                        links(n1, n2, n3, n4, mu) =
                            generate_from_gaussian(engine, dist);
                    }
                }
            }
        }
    }
}

void md_step(Configuration &links,
             Configuration &momenta,
             Configuration &momenta_half,
             std::mt19937 &engine,
             std::normal_distribution<double> &dist,
             double const time_step,
             double const beta) {
    // Update `momenta_half`.
    for (int n1 = 0; n1 != links.length_time; ++n1) {
        std::cout << "n1 = " << n1 << std::endl;
        for (int n2 = 0; n2 != links.length_space; ++n2) {
            for (int n3 = 0; n3 != links.length_space; ++n3) {
                for (int n4 = 0; n4 != links.length_space; ++n4) {
                    for (int mu = 0; mu != 4; ++mu) {
                        momenta_half(n1, n2, n3, n4, mu) =
                            compute_new_momentum(n1, n2, n3, n4, mu, links,
                                                 momenta, time_step, beta);
                    }
                }
            }
        }
    }
    // Update `links`.
    for (int n1 = 0; n1 != links.length_time; ++n1) {
        std::cout << "n1 = " << n1 << std::endl;
        for (int n2 = 0; n2 != links.length_space; ++n2) {
            for (int n3 = 0; n3 != links.length_space; ++n3) {
                for (int n4 = 0; n4 != links.length_space; ++n4) {
                    for (int mu = 0; mu != 4; ++mu) {
                        links(n1, n2, n3, n4, mu) =
                            compute_new_link(n1, n2, n3, n4, mu, links,
                                             momenta_half, time_step);
                    }
                }
            }
        }
    }
    // Update `momenta`.
    for (int n1 = 0; n1 != links.length_time; ++n1) {
        std::cout << "n1 = " << n1 << std::endl;
        for (int n2 = 0; n2 != links.length_space; ++n2) {
            for (int n3 = 0; n3 != links.length_space; ++n3) {
                for (int n4 = 0; n4 != links.length_space; ++n4) {
                    for (int mu = 0; mu != 4; ++mu) {
                        momenta(n1, n2, n3, n4, mu) =
                            compute_new_momentum(n1, n2, n3, n4, mu, links,
                                                 momenta_half, time_step, beta);
                    }
                }
            }
        }
    }
}

Eigen::Matrix2cd compute_new_momentum(int const n1,
                                           int const n2,
                                           int const n3,
                                           int const n4,
                                           int const mu,
                                           Configuration &links,
                                           Configuration &momenta,
                                           double const time_step,
                                           double const beta) {
    // Copy old momentum.
    Eigen::Matrix2cd result = momenta(n1, n2, n3, n4, mu);
    result += time_step / 2 *
              compute_momentum_derivative(n1, n2, n3, n4, mu, links, beta);
    return result;
}

Eigen::Matrix2cd compute_momentum_derivative(int const n1,
                                             int const n2,
                                             int const n3,
                                             int const n4,
                                             int const mu,
                                             Configuration &links,
                                             double const beta) {
    // Compute staples.
    Eigen::Matrix2cd staples;
    std::vector<int> const old_coords{n1, n2, n3, n4};
    for (int nu = 0; nu != 4; ++nu) {
        if (nu == mu) {
            continue;
        }
        auto coords = old_coords;

        auto &link3 = links(coords, nu);
        ++coords[mu];
        auto &link1 = links(coords, nu);
        --coords[mu];
        ++coords[nu];
        auto &link2 = links(coords, mu);

        staples += link1 * link2.adjoint() * link3.adjoint();

        coords = old_coords;

        --coords[nu];
        auto &link6 = links(coords, nu);
        auto &link5 = links(coords, mu);
        ++coords[mu];
        auto &link4 = links(coords, nu);

        staples += link4 * link5.adjoint() * link6.adjoint();
    }

    Eigen::Matrix2cd result = -beta / 6 * links(n1, n2, n3, n4, mu) * staples;
    result -= result.adjoint();

    return result;
}

Eigen::Matrix2cd compute_new_link(int const n1,
                                  int const n2,
                                  int const n3,
                                  int const n4,
                                  int const mu,
                                  Configuration &links,
                                  Configuration &momenta_half,
                                  double const time_step) {
    auto const exponent = std::complex<double>{0, 1} * time_step *
                          momenta_half(n1, n2, n3, n4, mu);
    auto const rotation = exponent.exp();
    return rotation * links(n1, n2, n3, n4, mu);
}
