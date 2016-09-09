// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

/// \file

#pragma once

#include "configuration.hpp"
#include "matrix.hpp"

#include <random>

void global_gauge_transformation(Matrix const &transformation, Configuration &links);

Matrix random_from_algebra(std::mt19937 &engine, std::normal_distribution<double> &dist);
Matrix random_from_group(std::mt19937 &engine, std::normal_distribution<double> &dist);
Matrix group_from_algebra(Matrix const &algebra_element);

/**
  Creates a SU(2) random configuration.
  */
Configuration make_hot_start(int const length_space,
                             int const length_time,
                             double const std,
                             int const seed);

/**
  Randomizes the whole lattice with su(2) matrices.
  */
void randomize_algebra(Configuration &configuration,
                       std::mt19937 &engine,
                       std::normal_distribution<double> &dist);

/**
  Randomizes the whole lattice with SU(2) matrices.
  */
void randomize_group(Configuration &configuration,
                     std::mt19937 &engine,
                     std::normal_distribution<double> &dist);

/**
  Perform a single molecular dynamics step.
  */
void md_step(Configuration &links,
             Configuration &momenta,
             Configuration &momenta_half,
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

double get_plaquette_trace_real(Configuration const &links);
std::complex<double> get_average_plaquette(Configuration const &links);

double get_energy(Configuration const &links, Configuration const &momenta);
double get_link_energy(Configuration const &links);
double get_momentum_energy(Configuration const &momenta);
