// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "hybrid-monte-carlo.hpp"

#include <iostream>
#include <random>

int main() {
    std::mt19937 engine;
    std::normal_distribution<double> dist(0, 1);

    std::cout << "Start" << std::endl;

    int const length_time = 10;
    int const length_space = 5;

    auto links = make_hot_start(length_space, length_time, 1, 0);

    std::cout << "Element: ";
    std::cout << links(0, 0, 0, 0, 0)(0, 0) << std::endl;

    auto momenta = make_hot_start(length_space, length_time, 1, 0);
    auto momenta_half = make_hot_start(length_space, length_time, 1, 0);

    const double time_step = 0.1;
    const double beta = 1.0;

    md_step(links, momenta, momenta_half, engine, dist, time_step, beta);

}
