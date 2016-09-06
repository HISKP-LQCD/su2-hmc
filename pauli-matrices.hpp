// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#pragma once

#include <Eigen/Dense>

#include <complex>
#include <vector>

class PauliMatrices {
    using value_type = std::complex<double>;

  public:
    Eigen::Matrix2cd const get(int i) const { return matrices[i]; }

    static PauliMatrices &get_instance() {
        static PauliMatrices instance;
        return instance;
    }

  private:
    PauliMatrices();

    std::vector<Eigen::Matrix2cd> matrices;
};
