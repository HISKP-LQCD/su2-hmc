// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "hybrid-monte-carlo.hpp"

#include <iostream>
#include <random>

int main() {
    std::mt19937 engine;
    std::normal_distribution<double> dist(0, 1);

    std::cout << "Start" << std::endl;

    auto links = make_hot_start(5, 10, 1, 0);

    std::cout << "Element: ";
    std::cout << links(0, 0, 0, 0, 0)(0, 0) << std::endl;
}
