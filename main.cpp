// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "hybrid-monte-carlo.hpp"

#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <random>

namespace ptree = boost::property_tree;

int main() {
    std::mt19937 engine;
    std::normal_distribution<double> dist(0, 1);

    std::cout << "sizeof(value_type): " << sizeof(Configuration::value_type)
              << std::endl;

    boost::format config_filename_format("gauge-links-%04d.bin");

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

    const double time_step = config.get<double>("md.time_step");
    const double beta = config.get<double>("md.beta");

    const int chain_total = config.get<int>("chain.total");
    const int chain_skip = config.get<int>("chain.skip");
    const int md_steps = config.get<int>("md.steps");

    Configuration momenta(length_space, length_time);
    Configuration momenta_half(length_space, length_time);

    int configs_stored = 0;
    int configs_computed = 0;

    while (configs_stored < chain_total) {
        Configuration const old_links = links;
        randomize(momenta, engine, dist);
        double const old_energy = get_energy(links, momenta);
        for (int md_step_idx = 0; md_step_idx != md_steps; ++md_step_idx) {
            md_step(links, momenta, momenta_half, engine, dist, time_step, beta);
        }
        double const new_energy = get_energy(links, momenta);

        std::cout << old_energy << "\t" << new_energy << std::endl;

        // TODO Accept-Reject.

        ++configs_computed;

        if (chain_skip == 0 || configs_computed % chain_skip == 0) {
            std::string filename = (config_filename_format % configs_stored).str();
            links.save(filename);
            ++configs_stored;
        }
    }


    std::cout << "Element: ";
    std::cout << links(0, 0, 0, 0, 0)(0, 0) << std::endl;

    links.save("links.bin");
}
