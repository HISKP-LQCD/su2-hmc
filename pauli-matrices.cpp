// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "pauli-matrices.hpp"

PauliMatrices::PauliMatrices() : matrices(3) {
    matrices[0] << value_type{0, 0}, value_type{1, 0}, value_type{1, 0}, value_type{1, 0};
    matrices[1] << value_type{0, 0}, value_type{0, -1}, value_type{0, 1},
        value_type{0, 0};
    matrices[2] << value_type{1, 0}, value_type{0, 0}, value_type{0, 0},
        value_type{-1, 0};
}
