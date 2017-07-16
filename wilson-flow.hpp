// Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

#pragma once

#include "configuration.hpp"

Configuration exp(Configuration const &links);

Configuration flow_step(Configuration const &initial, double const time_step);

Configuration get_global_stout_exponential(Configuration const &links);
