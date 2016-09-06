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

    const SU2Matrix result(exponentiated(0, 0), exponentiated(1, 0));
    return result;
}

Configuration make_hot_start(int const length_space,
                             int const length_time,
                             double const std,
                             int const seed) {
    std::mt19937 engine(seed);
    std::normal_distribution<double> dist(0, std);

    Configuration links(length_space, length_time);

    for (int n1 = 0; n1 != length_time; ++n1) {
        std::cout << "n1 = " << n1 << std::endl;
        for (int n2 = 0; n2 != length_space; ++n2) {
            for (int n3 = 0; n3 != length_space; ++n3) {
                for (int n4 = 0; n4 != length_space; ++n4) {
                    for (int mu = 0; mu != 4; ++mu) {
                        links(n1, n2, n3, n4, mu) =
                            generate_from_gaussian(engine, dist);
                    }
                }
            }
        }
    }

    return links;
}
