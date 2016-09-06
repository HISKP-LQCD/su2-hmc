// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "hybrid-monte-carlo.hpp"

#include <iostream>
#include <random>

int main() {
    std::mt19937 engine;
    std::normal_distribution<double> dist(0, 1);

    generate_from_gaussian(engine, dist);
}
