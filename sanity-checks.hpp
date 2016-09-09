// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#pragma once

#include <Eigen/Dense>

bool is_equal(Eigen::Matrix2cd const &mat1, Eigen::Matrix2cd const &mat2);
bool is_equal(double const d1, double const d2);
bool is_unit_determinant(Eigen::Matrix2cd const &mat);
bool is_hermitian(Eigen::Matrix2cd const &mat);
bool is_real(std::complex<double> const &c);
bool is_traceless(Eigen::Matrix2cd const &mat);
bool is_unitary(Eigen::Matrix2cd const &mat);
bool is_unity(Eigen::Matrix2cd const &mat);
bool is_zero(double const &d);
bool is_zero(Eigen::Matrix2cd const &mat);
bool is_zero(std::complex<double> const &c);
