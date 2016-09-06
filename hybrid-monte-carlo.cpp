// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "hybrid-monte-carlo.hpp"

#include "pauli-matrices.hpp"

#include <unsupported/Eigen/MatrixFunctions>

#include <iostream>

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
