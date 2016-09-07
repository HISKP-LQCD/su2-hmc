// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "configuration.hpp"

#include <fstream>
#include <iostream>

Configuration::Configuration(int const length_space, int const length_time)
    : length_space(length_space),
      length_time(length_time),
      spacing_n4(4),
      spacing_n3(length_space * spacing_n4),
      spacing_n2(length_space * spacing_n3),
      spacing_n1(length_space * spacing_n2),
      volume(length_space * length_space * length_space * length_time * 4),
      data(volume) {}

void Configuration::save(std::string const &path) {
    std::ofstream os(path, std::ios::out | std::ios::binary);
    os.write(reinterpret_cast<char *>(data.data()), storage_size());
}
