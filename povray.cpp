// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "hybrid-monte-carlo.hpp"
#include "pauli-matrices.hpp"

#include <cmath>
#include <iostream>

int get_color(Matrix const &link, int const index) {
    PauliMatrices const &pauli_matrices = PauliMatrices::get_instance();
    double const pi = std::acos(-1);

    Matrix const algebra = link.log();
    auto const trace = (algebra * pauli_matrices.get(index)).trace();
    auto const real = trace.imag();
    auto const abs = std::abs(real);

    auto const result = abs / (2 * pi) * 255;

    std::cerr << result << "\n";
    return result;
}

int main(int argc, char **argv) {
    Configuration links(8, 8);
    links.load(argv[1]);

    for (int n1 = 0; n1 < links.length_time; ++n1) {
        for (int n2 = 0; n2 < links.length_space; ++n2) {
            for (int n3 = 0; n3 < links.length_space; ++n3) {
                for (int n4 = 0; n4 < links.length_space; ++n4) {
                    double sum = 0.0;
                    for (int mu = 1; mu < 4; ++mu) {
                        for (int nu = 1; nu < 4; ++nu) {
                            sum += get_plaquette(n1, n2, n3, n4, mu, nu, links)
                                       .trace()
                                       .real();
                        }
                    }
                    std::cout << sum << "\t";
                }
            }
        }
    }

    std::cout << std::flush;
}
