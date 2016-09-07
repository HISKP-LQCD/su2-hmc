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
    } catch (ptree::ini_parser::ini_parser_error e) {
        std::cerr << e.what() << std::endl;
        abort();
    }

    std::cout << "Start" << std::endl;

    int const length_time = config.get<int>("lattice.length_time");
    int const length_space = config.get<int>("lattice.length_space");

    auto links = make_hot_start(length_space, length_time,
                                config.get<double>("init.hot_start_std"),
                                config.get<int>("init.seed"));

    std::cout << "Element: ";
    std::cout << links(0, 0, 0, 0, 0)(0, 0) << std::endl;

    auto momenta = make_hot_start(length_space, length_time, 1, 0);
    auto momenta_half = make_hot_start(length_space, length_time, 1, 0);

    const double time_step = config.get<double>("md.time_step");
    const double beta = config.get<double>("md.beta");

    md_step(links, momenta, momenta_half, engine, dist, time_step, beta);

    std::cout << "Element: ";
    std::cout << links(0, 0, 0, 0, 0)(0, 0) << std::endl;

    links.save("links.bin");
}
