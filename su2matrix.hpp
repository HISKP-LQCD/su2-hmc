// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#pragma once

#include <unsupported/Eigen/MatrixFunctions>

#include <complex>

/**
  SU(2) matrix class.

  Internally, it uses the following representation with \f$ a_1 \f$ and \f$ a_2
  \f$ complex numbers:
  \f[
      \begin{pmatrix}
          a_1 & - a_2^* \\
          a_2 & a_1^*
      \end{pmatrix} \,.
  \f]

  https://en.wikipedia.org/wiki/Special_unitary_group#n_.3D_2
  */
class SU2Matrix {
  public:
    using value_type = std::complex<double>;

    SU2Matrix() : a1(value_type{0, 0}), a2(value_type{0, 0}){};
    SU2Matrix(value_type const &a1, value_type const &a2) : a1(a1), a2(a2){};

    value_type operator()(int const row, int const col) const;

    SU2Matrix operator+(SU2Matrix const &right) const;
    SU2Matrix operator-(SU2Matrix const &right) const;
    SU2Matrix operator*(SU2Matrix const &right) const;

  private:
    value_type a1, a2;
};

SU2Matrix const zero(SU2Matrix::value_type{0, 0}, SU2Matrix::value_type{0, 0});
SU2Matrix const unity(SU2Matrix::value_type{1, 0}, SU2Matrix::value_type{0, 0});

SU2Matrix make_su2matrix(Eigen::Matrix2cd const &mat);

Eigen::Matrix2cd make_eigen(SU2Matrix const &mat);

SU2Matrix exp(SU2Matrix const &mat);
