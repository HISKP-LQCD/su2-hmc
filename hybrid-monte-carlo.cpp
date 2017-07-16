// Copyright © 2016-2017 Martin Ueding <dev@martin-ueding.de>

#include "hybrid-monte-carlo.hpp"

#include "pauli-matrices.hpp"
#include "sanity-checks.hpp"

#include <unsupported/Eigen/MatrixFunctions>

#include <cassert>
#include <iostream>
#include <random>

double md_evolution(Configuration &links,
                    std::mt19937 &engine,
                    std::normal_distribution<double> &dist,
                    double const time_step,
                    int const md_steps,
                    double const beta) {
    Configuration momenta(links.length_space, links.length_time);
    randomize_algebra(momenta, engine, dist);

    double const old_energy = get_energy(links, momenta, beta);

    md_momentum_step(links, momenta, engine, dist, time_step / 2, beta);
    md_link_step(links, momenta, engine, dist, time_step, beta);
    for (int md_step_idx = 1; md_step_idx != md_steps; ++md_step_idx) {
        md_momentum_step(links, momenta,engine, dist, time_step, beta);
        md_link_step(links, momenta, engine, dist, time_step, beta);
    }
    md_momentum_step(links, momenta, engine, dist, time_step / 2, beta);

    double const new_energy = get_energy(links, momenta, beta);
    double const energy_difference = new_energy - old_energy;

    std::cout << "HMD Energy: " << old_energy << " → " << new_energy
              << "; ΔE = " << energy_difference << std::endl;

    return energy_difference;
}

void md_link_step(Configuration &links,
             Configuration &momenta,
             std::mt19937 &engine,
             std::normal_distribution<double> &dist,
             double const time_step,
             double const beta) {
    auto const links_old = links;
#pragma omp parallel for collapse(4)
    for (int n1 = 0; n1 < links.length_time; ++n1) {
        for (int n2 = 0; n2 < links.length_space; ++n2) {
            for (int n3 = 0; n3 < links.length_space; ++n3) {
                for (int n4 = 0; n4 < links.length_space; ++n4) {
                    for (int mu = 0; mu < 4; ++mu) {
                        links(n1, n2, n3, n4, mu) = compute_new_link(
                            n1, n2, n3, n4, mu, links_old, momenta, time_step);
                    }
                }
            }
        }
    }
}

void md_momentum_step(Configuration &links,
                      Configuration &momenta,
                      std::mt19937 &engine,
                      std::normal_distribution<double> &dist,
                      double const time_step,
                      double const beta) {
    auto const momenta_old = momenta;
#pragma omp parallel for collapse(4)
    for (int n1 = 0; n1 < links.length_time; ++n1) {
        for (int n2 = 0; n2 < links.length_space; ++n2) {
            for (int n3 = 0; n3 < links.length_space; ++n3) {
                for (int n4 = 0; n4 < links.length_space; ++n4) {
                    for (int mu = 0; mu < 4; ++mu) {
                        momenta(n1, n2, n3, n4, mu) = compute_new_momentum(
                            n1, n2, n3, n4, mu, links, momenta_old, time_step, beta);
                    }
                }
            }
        }
    }
}

Matrix compute_new_momentum(int const n1,
                                      int const n2,
                                      int const n3,
                                      int const n4,
                                      int const mu,
                                      Configuration const &links,
                                      Configuration const &momenta,
                                      double const time_step,
                                      double const beta) {
    // Copy old momentum.
    Matrix result = momenta(n1, n2, n3, n4, mu);
    assert(is_traceless(result));
    assert(is_hermitian(result));
    result += time_step * compute_momentum_derivative(n1, n2, n3, n4, mu, links, beta);
    assert(is_traceless(result));
    assert(is_hermitian(result));
    return result;
}

Matrix get_staples(int const n1,
                             int const n2,
                             int const n3,
                             int const n4,
                             int const mu,
                             Configuration const &links) {
    Matrix staples;
    // XXX Perhaps use std::array here.
    std::vector<int> const old_coords{n1, n2, n3, n4};
    for (int nu = 0; nu < 4; ++nu) {
        if (nu == mu) {
            continue;
        }
        auto coords = old_coords;

        Matrix const &link3 = links(coords, nu);
        ++coords[mu];
        Matrix const &link1 = links(coords, nu);
        --coords[mu];
        ++coords[nu];
        Matrix const &link2 = links(coords, mu);

        staples += link1 * link2.adjoint() * link3.adjoint();

        coords = old_coords;

        --coords[nu];
        Matrix const &link6 = links(coords, nu);
        Matrix const &link5 = links(coords, mu);
        ++coords[mu];
        Matrix const &link4 = links(coords, nu);

        staples += link4.adjoint() * link5.adjoint() * link6;
    }
    return staples;
}

Matrix get_stout_exponential(int const n1,
                             int const n2,
                             int const n3,
                             int const n4,
                             int const mu,
                             Configuration const &links) {
    Matrix staples = get_staples(n1, n2, n3, n4, mu, links);
    Matrix const links_staples = links(n1, n2, n3, n4, mu) * staples;
    Matrix const minus_adjoint = links_staples - links_staples.adjoint().eval();
    return minus_adjoint;
}

Matrix compute_momentum_derivative(int const n1,
                                   int const n2,
                                   int const n3,
                                   int const n4,
                                   int const mu,
                                   Configuration const &links,
                                   double const beta) {
    auto minus_adjoint = get_stout_exponential(n1, n2, n3, n4, mu, links);
    Matrix const derivative =
        imag_unit * beta / static_cast<double>(number_of_colors) * minus_adjoint;

    assert(is_traceless(derivative));
    assert(is_hermitian(derivative));

    return derivative;
}

Matrix compute_new_link(int const n1,
                        int const n2,
                        int const n3,
                        int const n4,
                        int const mu,
                        Configuration const &links,
                        Configuration const &momenta_half,
                        double const time_step) {
    Matrix const exponent = imag_unit * time_step * momenta_half(n1, n2, n3, n4, mu);
    Matrix const rotation = exponent.exp();
    assert(is_unitary(rotation));

    Matrix const new_link = rotation * links(n1, n2, n3, n4, mu);
    assert(is_unitary(new_link));

    return new_link;
}

Matrix get_plaquette(int const n1,
                     int const n2,
                     int const n3,
                     int const n4,
                     int const mu,
                     int const nu,
                     Configuration const &links,
                     bool const debug) {
    std::vector<int> coords{n1, n2, n3, n4};

    Matrix const &link1 = links(coords, mu);
    Matrix const &link4 = links(coords, nu);
    ++coords[mu];
    Matrix const &link2 = links(coords, nu);
    --coords[mu];
    ++coords[nu];
    Matrix const &link3 = links(coords, mu);

    if (debug) {
        std::cout << "link1:\n" << link1 << "\n";
        std::cout << "link2:\n" << link2 << "\n";
        std::cout << "link3:\n" << link3 << "\n";
        std::cout << "link4:\n" << link4 << "\n";
    }

    Matrix const plaquette = link1 * link2 * link3.adjoint() * link4.adjoint();
    assert(is_unitary(plaquette));
    return plaquette;
}

std::complex<double> get_plaquette_trace_sum(Configuration const &links) {
    double real = 0.0;
    double imag = 0.0;

#pragma omp parallel for reduction(+ : real, imag) collapse(4)
    for (int n1 = 0; n1 < links.length_time; ++n1) {
        for (int n2 = 0; n2 < links.length_space; ++n2) {
            for (int n3 = 0; n3 < links.length_space; ++n3) {
                for (int n4 = 0; n4 < links.length_space; ++n4) {
                    for (int mu = 0; mu < 4; ++mu) {
                        for (int nu = 0; nu < mu; ++nu) {
                            Matrix const plaquette =
                                get_plaquette(n1, n2, n3, n4, mu, nu, links);
                            auto const summand = plaquette.trace();
                            assert(std::isfinite(summand.real()));
                            assert(std::isfinite(summand.imag()));
                            real += summand.real();
                            imag += summand.imag();
                        }
                    }
                }
            }
        }
    }

    return std::complex<double>{real, imag};
}

std::complex<double> get_plaquette_trace_average(Configuration const &links) {
    double const summands = links.get_volume() * (3 + 2 + 1) * number_of_colors;
    return get_plaquette_trace_sum(links) / summands;
}

double get_link_energy(Configuration const &links, double const beta) {
    double links_part = 0.0;
    links_part += links.get_volume() * (3 + 2 + 1);
    links_part -= get_plaquette_trace_sum(links).real();
    return links_part * beta;
}

double get_momentum_energy(Configuration const &momenta, double const beta) {
    double momentum_part = 0.0;
#pragma omp parallel for reduction(+ : momentum_part) collapse(4)
    for (int n1 = 0; n1 < momenta.length_time; ++n1) {
        for (int n2 = 0; n2 < momenta.length_space; ++n2) {
            for (int n3 = 0; n3 < momenta.length_space; ++n3) {
                for (int n4 = 0; n4 < momenta.length_space; ++n4) {
                    for (int mu = 0; mu < 4; ++mu) {
                        Matrix const &momentum = momenta(n1, n2, n3, n4, mu);
                        auto const trace = (momentum * momentum).trace();
                        assert(is_real(trace));
                        auto const summand = trace.real();
                        assert(std::isfinite(summand));
                        momentum_part += summand;
                    }
                }
            }
        }
    }
     return 0.5 * momentum_part;
}

double
get_energy(Configuration const &links, Configuration const &momenta, double const beta) {
    return get_link_energy(links, beta) + get_momentum_energy(momenta, beta);
}
