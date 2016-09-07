// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "hybrid-monte-carlo.hpp"

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <random>

namespace ptree = boost::property_tree;

int main() {
    std::mt19937 engine;
    std::normal_distribution<double> dist(0, 1);

    ptree::ptree config;
    try {
    ptree::read_ini("hmc.ini", config);
    }
    catch (ptree::ini_parser::ini_parser_error e) {
        std::cerr << e.what() << std::endl;
        abort();
    }

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

    std::cout << "Element: ";
    std::cout << links(0, 0, 0, 0, 0)(0, 0) << std::endl;

}
