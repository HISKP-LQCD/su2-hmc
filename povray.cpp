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

int main() {
    Configuration links(8, 8);
    links.load("gauge-links-0000.bin");

    int const max = 5;


    // for (int n1 = 0; n1 < links.length_time; ++n1) {
    for (int n2 = 0; n2 < max; ++n2) {
        for (int n3 = 0; n3 < max; ++n3) {
            for (int n4 = 0; n4 < max; ++n4) {
                for (int mu = 1; mu < 4; ++mu) {
                    int const x = n2;
                    int const y = n3;
                    int const z = n4;

                    auto const red = get_color(links(0, n2, n3, n4, mu), 0);
                    auto const green = get_color(links(0, n2, n3, n4, mu), 1);
                    auto const blue = get_color(links(0, n2, n3, n4, mu), 2);

                    std::cout << "cylinder {\n"
                                 "<0, 0, 0> <0, 0, 1> rad\n"
                                 "pigment { color rgb<"
                              << red << ", " << green << ", " << blue << "> }\n";
                    if (mu == 2) {
                        std::cout << "rotate <-90, 0, 0>\n";
                    } else if (mu == 3) {
                        std::cout << "rotate <0, 90, 0>\n";
                    }
                    std::cout << "translate <" << x << ", " << y << ", " << z << ">\n";

                    //std::cout << "finish {phong 1}\n";
                    std::cout << "no_shadow\n"
                                 "} \n";
                }
            }
        }
    }

    std::cout << std::flush;
}
