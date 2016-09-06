// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "su2matrix.hpp"

#include <cassert>

SU2Matrix SU2Matrix::operator+(SU2Matrix const &right) {
    SU2Matrix result;
    result.a1 = a1 + right.a1;
    result.a2 = a2 + right.a2;
    return result;
}

SU2Matrix SU2Matrix::operator-(SU2Matrix const &right) {
    SU2Matrix result;
    result.a1 = a1 - right.a1;
    result.a2 = a2 - right.a2;
    return result;
}

SU2Matrix SU2Matrix::operator*(SU2Matrix const &right) {
    SU2Matrix result;
    result.a1 = a1 * right.a1 - std::conj(a2) * right.a2;
    result.a2 = a2 * right.a1 + std::conj(a1) * right.a2;
    return result;
}

std::complex<double> SU2Matrix::operator()(int const row, int const col) {
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
