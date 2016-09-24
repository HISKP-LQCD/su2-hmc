// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#pragma once

#include "matrix.hpp"

#include <complex>
#include <random>
#include <vector>

class PauliMatrices {
    using value_type = std::complex<double>;

  public:
    Matrix const get(int i) const { return matrices[i]; }

    static PauliMatrices &get_instance() {
        static PauliMatrices instance;
        return instance;
    }

  private:
    PauliMatrices();

    std::vector<Matrix> matrices;
};

Matrix random_from_algebra(std::mt19937 &engine, std::normal_distribution<double> &dist);
Matrix random_from_group(std::mt19937 &engine, std::normal_distribution<double> &dist);
Matrix group_from_algebra(Matrix const &algebra_element);
