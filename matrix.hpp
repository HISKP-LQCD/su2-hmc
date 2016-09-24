// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#pragma once

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <complex>

using Matrix = Eigen::Matrix2cd;
using Complex = std::complex<double>;

int constexpr number_of_colors = 2;
Complex constexpr imag_unit{0, 1};
