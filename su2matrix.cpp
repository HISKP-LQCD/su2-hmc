// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "su2matrix.hpp"

#include <unsupported/Eigen/MatrixFunctions>

#include <cassert>

SU2Matrix SU2Matrix::operator+(SU2Matrix const &right) const {
    SU2Matrix result;
    result.a1 = a1 + right.a1;
    result.a2 = a2 + right.a2;
    return result;
}

SU2Matrix SU2Matrix::operator-(SU2Matrix const &right) const {
    SU2Matrix result;
    result.a1 = a1 - right.a1;
    result.a2 = a2 - right.a2;
    return result;
}

SU2Matrix SU2Matrix::operator*(SU2Matrix const &right) const {
    SU2Matrix result;
    result.a1 = a1 * right.a1 - std::conj(a2) * right.a2;
    result.a2 = a2 * right.a1 + std::conj(a1) * right.a2;
    return result;
}

std::complex<double> SU2Matrix::operator()(int const row, int const col) const {
    assert(row == 0 || row == 1);
    assert(col == 0 || col == 1);

    if (row == 0 && col == 0)
        return a1;
    else if (row == 0 && col == 1)
        return -std::conj(a2);
    else if (row == 1 && col == 0)
        return a2;
    else
        return std::conj(a1);
}

SU2Matrix make_su2matrix(Eigen::Matrix2cd const &mat) {
    return {mat(0, 0), mat(1, 0)};
}

Eigen::Matrix2cd make_eigen(SU2Matrix const &mat) {
    Eigen::Matrix2cd out;
    out << mat(0, 0), mat(0, 1), mat(1, 0), mat(0, 1);
    return out;
}

SU2Matrix exp(SU2Matrix const &mat) {
    return make_su2matrix(make_eigen(mat).exp());
}
