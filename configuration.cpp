// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "configuration.hpp"

Configuration::Configuration(int const length_space, int const length_time)
    : length_space(length_space),
      length_time(length_time),
      spacing_n4(4),
      spacing_n3(length_space * spacing_n4),
      spacing_n2(length_space * spacing_n3),
      spacing_n1(length_space * spacing_n2),
      volume(length_space * length_space * length_space * length_time) {}
