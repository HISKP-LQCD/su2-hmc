// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#pragma once

#include "matrix.hpp"


bool is_equal(double const d1, double const d2);
bool is_equal(Matrix const &mat1, Matrix const &mat2);
bool is_hermitian(Matrix const &mat);
bool is_real(std::complex<double> const &c);
bool is_traceless(Matrix const &mat);
bool is_unit_determinant(Matrix const &mat);
bool is_unitary(Matrix const &mat);
bool is_unity(Matrix const &mat);
bool is_zero(double const &d);
bool is_zero(Matrix const &mat);
bool is_zero(std::complex<double> const &c);
