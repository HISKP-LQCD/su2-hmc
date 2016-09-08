// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "hybrid-monte-carlo.hpp"

#include "pauli-matrices.hpp"
#include "sanity-checks.hpp"

#include <unsupported/Eigen/MatrixFunctions>

#include <cassert>
#include <iostream>
#include <random>

std::complex<double> constexpr imag_unit{0, 1};

void global_gauge_transformation(Matrix const &transformation, Configuration &links) {
    for (int i = 0; i < links.get_size(); ++i) {
        links[i] *= transformation;
    }
}

Matrix random_from_algebra(std::mt19937 &engine, std::normal_distribution<double> &dist) {
    std::vector<double> coefficients = {dist(engine), dist(engine), dist(engine)};
    auto const pauli_matrices = PauliMatrices::get_instance();
    Matrix algebra_element;
    for (int i = 0; i < 3; ++i) {
        Matrix scaled_generator = coefficients[i] * pauli_matrices.get(i);
        algebra_element += scaled_generator;
    }

    assert(is_hermitian(algebra_element));
    assert(is_traceless(algebra_element));
    return algebra_element;
}

Matrix group_from_algebra(Matrix const &algebra_element) {
    auto const exponent = imag_unit * algebra_element;
    auto const group_element = exponent.exp();
    assert(is_unitary(group_element));
    return group_element;
}

Matrix random_from_group(std::mt19937 &engine, std::normal_distribution<double> &dist) {
    auto const algebra_element = random_from_algebra(engine, dist);
    auto const group_element = group_from_algebra(algebra_element);
    assert(is_unitary(group_element));
    return group_element;
}

Configuration make_hot_start(int const length_space,
                             int const length_time,
                             double const std,
                             int const seed) {
    std::mt19937 engine(seed);
    std::normal_distribution<double> dist(0, std);

    Configuration links(length_space, length_time);
    randomize_group(links, engine, dist);
    return links;
}

void randomize_algebra(Configuration &config,
                       std::mt19937 &engine,
                       std::normal_distribution<double> &dist) {
    for (int i = 0; i < config.get_size(); ++i) {
        auto const next = random_from_algebra(engine, dist);
        config[i] = next;
    }
}

void randomize_group(Configuration &config,
                     std::mt19937 &engine,
                     std::normal_distribution<double> &dist) {
    for (int i = 0; i < config.get_size(); ++i) {
        auto const next = random_from_group(engine, dist);
        config[i] = next;
    }
}

void md_step(Configuration &links,
             Configuration &momenta,
             Configuration &momenta_half,
             std::mt19937 &engine,
             std::normal_distribution<double> &dist,
             double const time_step,
             double const beta) {
// Update `momenta_half`.
#pragma omp parallel for
    for (int n1 = 0; n1 < links.length_time; ++n1) {
        for (int n2 = 0; n2 < links.length_space; ++n2) {
            for (int n3 = 0; n3 < links.length_space; ++n3) {
                for (int n4 = 0; n4 < links.length_space; ++n4) {
                    for (int mu = 0; mu < 4; ++mu) {
                        momenta_half(n1, n2, n3, n4, mu) = compute_new_momentum(
                            n1, n2, n3, n4, mu, links, momenta, time_step, beta);
                    }
                }
            }
        }
    }
// Update `links`.
#pragma omp parallel for
    for (int n1 = 0; n1 < links.length_time; ++n1) {
        for (int n2 = 0; n2 < links.length_space; ++n2) {
            for (int n3 = 0; n3 < links.length_space; ++n3) {
                for (int n4 = 0; n4 < links.length_space; ++n4) {
                    for (int mu = 0; mu < 4; ++mu) {
                        links(n1, n2, n3, n4, mu) = compute_new_link(
                            n1, n2, n3, n4, mu, links, momenta_half, time_step);
                    }
                }
            }
        }
    }
// Update `momenta`.
#pragma omp parallel for
    for (int n1 = 0; n1 < links.length_time; ++n1) {
        for (int n2 = 0; n2 < links.length_space; ++n2) {
            for (int n3 = 0; n3 < links.length_space; ++n3) {
                for (int n4 = 0; n4 < links.length_space; ++n4) {
                    for (int mu = 0; mu < 4; ++mu) {
                        momenta(n1, n2, n3, n4, mu) = compute_new_momentum(
                            n1, n2, n3, n4, mu, links, momenta_half, time_step, beta);
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
    result +=
        time_step / 2 * compute_momentum_derivative(n1, n2, n3, n4, mu, links, beta);
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

        auto &link3 = links(coords, nu);
        ++coords[mu];
        auto &link1 = links(coords, nu);
        --coords[mu];
        ++coords[nu];
        auto &link2 = links(coords, mu);

        staples += link1 * link2.adjoint() * link3.adjoint();

        coords = old_coords;

        --coords[nu];
        auto &link6 = links(coords, nu);
        auto &link5 = links(coords, mu);
        ++coords[mu];
        auto &link4 = links(coords, nu);

        staples += link4.adjoint() * link5.adjoint() * link6;
    }
    return staples;
}

Matrix compute_momentum_derivative(int const n1,
                                             int const n2,
                                             int const n3,
                                             int const n4,
                                             int const mu,
                                             Configuration const &links,
                                             double const beta) {
    auto staples = get_staples(n1, n2, n3, n4, mu, links);
    Matrix const links_staples = links(n1, n2, n3, n4, mu) * staples;
    Matrix const minus_adjoint = links_staples - links_staples.adjoint().eval();
    Matrix const derivative = imag_unit * beta / 6.0 * minus_adjoint;

#ifndef NDEBUG
    if (!is_traceless(derivative)) {
        std::cerr << "Momentum is not traceless!\n";
        std::cerr << n1 << ", " << n2 << ", " << n3 << ", " << n4 << "; " << mu << "\n";
        std::cerr << "UV:\n";
        std::cerr << links_staples << "\n";
        std::cerr << "UV - (UV)^\\dagger:\n";
        std::cerr << minus_adjoint << "\n";
        std::cerr << "\\dot H:\n";
        std::cerr << derivative << std::endl;
    }
#endif

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
    auto const exponent = imag_unit * time_step * momenta_half(n1, n2, n3, n4, mu);
    auto const rotation = exponent.exp();
    assert(is_unitary(rotation));

    auto const new_link = rotation * links(n1, n2, n3, n4, mu);
    assert(is_unitary(new_link));

    return new_link;
}

Matrix get_plaquette(int const n1,
                               int const n2,
                               int const n3,
                               int const n4,
                               int const mu,
                               int const nu,
                               Configuration const &links) {
    std::vector<int> coords{n1, n2, n3, n4};

    auto &link1 = links(coords, mu);
    auto &link4 = links(coords, nu);
    ++coords[mu];
    auto &link2 = links(coords, nu);
    --coords[mu];
    ++coords[nu];
    auto &link3 = links(coords, mu);

    Matrix const plaquette = link1 * link2 * link3.adjoint() * link4.adjoint();
    assert(is_unitary(plaquette));
    return plaquette;
}

double get_plaquette_trace_real(Configuration const &links) {
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (int n1 = 0; n1 < links.length_time; ++n1) {
        for (int n2 = 0; n2 < links.length_space; ++n2) {
            for (int n3 = 0; n3 < links.length_space; ++n3) {
                for (int n4 = 0; n4 < links.length_space; ++n4) {
                    for (int mu = 0; mu < 4; ++mu) {
                        for (int nu = 0; nu < 4; ++nu) {
                            auto const plaquette =
                                get_plaquette(n1, n2, n3, n4, mu, nu, links);
                            double const summand = plaquette.trace().real();
                            assert(std::isfinite(summand));
                            sum += summand;
                        }
                    }
                }
            }
        }
    }

    return sum;
}

double get_link_energy(Configuration const &links) {
    double links_part = 0.0;
    links_part += links.get_volume() * 4 * 4;
    links_part -= get_plaquette_trace_real(links);
    return links_part;
}

double get_momentum_energy(Configuration const &momenta) {
    double momentum_part = 0.0;
#pragma omp parallel for reduction(+ : momentum_part)
    for (int n1 = 0; n1 < momenta.length_time; ++n1) {
        for (int n2 = 0; n2 < momenta.length_space; ++n2) {
            for (int n3 = 0; n3 < momenta.length_space; ++n3) {
                for (int n4 = 0; n4 < momenta.length_space; ++n4) {
                    for (int mu = 0; mu < 4; ++mu) {
                        auto const &momentum = momenta(n1, n2, n3, n4, mu);
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

double get_energy(Configuration const &links, Configuration const &momenta) {
    return get_link_energy(links) + get_momentum_energy(momenta);
}
