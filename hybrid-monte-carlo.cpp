// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "hybrid-monte-carlo.hpp"

#include "pauli-matrices.hpp"

#include <unsupported/Eigen/MatrixFunctions>

#include <iostream>
#include <random>

SU2Matrix generate_from_gaussian(std::mt19937 &engine,
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

    //std::cout << algebra_element << std::endl;

    Eigen::Matrix2cd const exponentiated = algebra_element.exp();

    return make_su2matrix(exponentiated);
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
             std::normal_distribution<double> &dist) {
    // Update `momenta_half`.
    for (int n1 = 0; n1 != links.length_time; ++n1) {
        std::cout << "n1 = " << n1 << std::endl;
        for (int n2 = 0; n2 != links.length_space; ++n2) {
            for (int n3 = 0; n3 != links.length_space; ++n3) {
                for (int n4 = 0; n4 != links.length_space; ++n4) {
                    for (int mu = 0; mu != 4; ++mu) {
                        momenta_half(n1, n2, n3, n4, mu) =
                            compute_new_momentum_half(n1, n2, n3, n4, mu, links, momenta)
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
                            compute_new_link(n1, n2, n3, n4, mu, links, momenta)
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
                            compute_momentum_half(n1, n2, n3, n4, mu, links, momenta)
                    }
                }
            }
        }
    }
}
