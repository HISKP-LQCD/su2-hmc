// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#pragma once

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
    std::complex<double> operator()(int const row, int const col);

    SU2Matrix operator+(SU2Matrix const &right);
    SU2Matrix operator-(SU2Matrix const &right);
    SU2Matrix operator*(SU2Matrix const &right);

  private:
    std::complex<double> a1, a2;
};
