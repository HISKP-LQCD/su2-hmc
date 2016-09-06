// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#pragma once

#include "configuration.hpp"

#include <random>

SU2Matrix generate_from_gaussian(std::mt19937 &engine,
                                 std::normal_distribution<double> &dist);
