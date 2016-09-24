// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

/// \file

#pragma once

#include "configuration.hpp"
#include "matrix.hpp"

#include <random>

double md_evolution(Configuration &links,
                    std::mt19937 &engine,
                    std::normal_distribution<double> &dist,
                    double const time_step,
                    int const md_steps,
                    double const beta);

void md_momentum_half_step(Configuration &links,
                           Configuration &momenta,
                           std::mt19937 &engine,
                           std::normal_distribution<double> &dist,
                           double const time_step,
                           double const beta);

void md_link_step(Configuration &links,
                  Configuration &momenta,
                  std::mt19937 &engine,
                  std::normal_distribution<double> &dist,
                  double const time_step,
                  double const beta);

void md_momentum_step(Configuration &links,
                      Configuration &momenta,
                      std::mt19937 &engine,
                      std::normal_distribution<double> &dist,
                      double const time_step,
                      double const beta);

/**
  Compute the new momentum half a timestep further.
  */
Matrix compute_new_momentum(int const n1,
                            int const n2,
                            int const n3,
                            int const n4,
                            int const mu,
                            Configuration const &links,
                            Configuration const &momenta,
                            double const time_step,
                            double const beta);

/**
  Compute the staples for a given link.
  */
Matrix get_staples(int const n1,
                   int const n2,
                   int const n3,
                   int const n4,
                   int const mu,
                   Configuration const &links);

Matrix compute_momentum_derivative(int const n1,
                                   int const n2,
                                   int const n3,
                                   int const n4,
                                   int const mu,
                                   Configuration const &links,
                                   double const beta);

Matrix compute_new_link(int const n1,
                        int const n2,
                        int const n3,
                        int const n4,
                        int const mu,
                        Configuration const &links,
                        Configuration const &momenta_half,
                        double const time_step);

Matrix get_plaquette(int const n1,
                     int const n2,
                     int const n3,
                     int const n4,
                     int const mu,
                     int const nu,
                     Configuration const &links,
                     bool const debug = false);

std::complex<double> get_plaquette_trace_sum(Configuration const &links);
std::complex<double> get_plaquette_trace_average(Configuration const &links);

double get_energy(Configuration const &links, Configuration const &momenta, double const beta);
double get_link_energy(Configuration const &links, double const beta);
double get_momentum_energy(Configuration const &momenta, double const beta);
