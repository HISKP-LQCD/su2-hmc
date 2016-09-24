// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "pauli-matrices.hpp"

#include "sanity-checks.hpp"

PauliMatrices::PauliMatrices() : matrices(3) {
    matrices[0] << value_type{0, 0}, value_type{1, 0}, value_type{1, 0}, value_type{0, 0};
    matrices[1] << value_type{0, 0}, value_type{0, -1}, value_type{0, 1},
        value_type{0, 0};
    matrices[2] << value_type{1, 0}, value_type{0, 0}, value_type{0, 0},
        value_type{-1, 0};
}

Matrix random_from_algebra(std::mt19937 &engine, std::normal_distribution<double> &dist) {
    std::vector<double> coefficients = {dist(engine), dist(engine), dist(engine)};
    PauliMatrices const &pauli_matrices = PauliMatrices::get_instance();
    Matrix algebra_element;
    for (int i = 0; i < 3; ++i) {
        Matrix scaled_generator = coefficients[i] * pauli_matrices.get(i);
        algebra_element += scaled_generator;
    }

    assert(is_hermitian(algebra_element));
    assert(is_traceless(algebra_element));
    return algebra_element;
}

Matrix group_from_algebra(Matrix const &algebra_element) {
    Matrix const exponent = imag_unit * algebra_element;
    Matrix const group_element = exponent.exp();
    assert(is_unitary(group_element));
    return group_element;
}

Matrix random_from_group(std::mt19937 &engine, std::normal_distribution<double> &dist) {
    Matrix const algebra_element = random_from_algebra(engine, dist);
    Matrix const group_element = group_from_algebra(algebra_element);
    assert(is_unitary(group_element));
    return group_element;
}
