// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "hybrid-monte-carlo.hpp"
#include "sanity-checks.hpp"

#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <random>

#define OUTPUT

namespace ptree = boost::property_tree;

int main() {
    std::mt19937 engine;
    std::normal_distribution<double> dist(0, 1);
    std::uniform_real_distribution<double> uniform(0, 1);

    std::cout << "sizeof(value_type): " << sizeof(Configuration::value_type) << std::endl;

    boost::format config_filename_format("gauge-links-%04d.bin");

    ptree::ptree config;
    try {
        ptree::read_ini("hmc.ini", config);
    } catch (ptree::ini_parser::ini_parser_error e) {
        std::cerr << e.what() << std::endl;
        abort();
    }

    int const length_time = config.get<int>("lattice.length_time");
    int const length_space = config.get<int>("lattice.length_space");

    bool const do_write_config = config.get<bool>("output.links", false);

    auto links = make_hot_start(length_space, length_time,
                                config.get<double>("init.hot_start_std"),
                                config.get<int>("init.seed"));

    const double time_step = config.get<double>("md.time_step");
    const double beta = config.get<double>("md.beta");

    const int chain_total = config.get<int>("chain.total");
    const int chain_skip = config.get<int>("chain.skip");
    const int md_steps = config.get<int>("md.steps");

    Configuration momenta(length_space, length_time);
    Configuration momenta_half(length_space, length_time);

    std::ofstream ofs_energy("energy.tsv");
    std::ofstream ofs_plaquette("plaquette.tsv");
    std::ofstream ofs_energy_reject("energy-reject.tsv");
    std::ofstream ofs_plaquette_reject("plaquette-reject.tsv");
    std::ofstream ofs_boltzmann("boltzmann.tsv");
    std::ofstream ofs_accept("accept.tsv");

    int number_stored = 0;
    int number_computed = 0;
    int number_accepted = 0;

    while (number_computed < chain_total) {
        Configuration const old_links = links;

        // FIXME The standard deviation here should depend on the time step.
        randomize_algebra(momenta, engine, dist);

        double factor_sum = 0.0;
        int factor_count = 0;

        double const old_energy = get_energy(links, momenta, beta);
        for (int md_step_idx = 0; md_step_idx != md_steps; ++md_step_idx) {
            double const old_links_energy = get_link_energy(links, beta);
            double const old_momentum_energy = get_momentum_energy(momenta, beta);
            md_step(links, momenta, momenta_half, engine, dist, time_step, beta);
            double const new_links_energy = get_link_energy(links, beta);
            double const new_momentum_energy = get_momentum_energy(momenta, beta);

            double const links_energy_difference = new_links_energy - old_links_energy;
            double const momentum_energy_difference = new_momentum_energy - old_momentum_energy;

            double const factor = -links_energy_difference / momentum_energy_difference;

#ifdef OUTPUT
            double const energy_difference =
                links_energy_difference + momentum_energy_difference;

            std::cout << "ΔMD Energy: Links = " << links_energy_difference
                      << ", momentum = " << momentum_energy_difference
                      << ", total = " << energy_difference << ", ratio = " << factor
                      << std::endl;
#endif

            factor_sum += factor;
            ++factor_count;

            if (factor_count > 5 && std::abs(factor_sum / factor_count - 1) > 0.1) {
                std::cerr
                    << "WARNUNG: Link and momentum energy transfer does not match up, factor is "
                    << factor << std::endl;
                //abort();
            }
        }

        double const new_energy = get_energy(links, momenta, beta);
        double const energy_difference = new_energy - old_energy;

        ofs_boltzmann << energy_difference << std::endl;

        std::cout << "HMD Energy: " << old_energy << " → " << new_energy
                  << "; ΔE = " << energy_difference << std::endl;

        double const average_plaquette = get_plaquette_trace_average(links).real();

#ifdef OUTPUT
        auto const transformation = random_from_group(engine, dist);
        global_gauge_transformation(transformation, links);
        double const average_plaquette_2 = get_plaquette_trace_average(links).real();
        double const plaquette_diference = average_plaquette_2 - average_plaquette;
        std::cout << "Global transform: " << average_plaquette << " → "
                  << average_plaquette_2 << "; Δ = " << plaquette_diference << std::endl;
        assert(is_equal(average_plaquette, average_plaquette_2));
#endif

        ++number_computed;


        // Accept-Reject.
        const bool accepted =
            energy_difference <= 0 || std::exp(-energy_difference) >= uniform(engine);
        if (accepted) {
            std::cout << "Accepted.\n";

            ofs_energy << number_computed << "\t" << (new_energy / links.get_volume())
                       << std::endl;
            ofs_plaquette << number_computed << "\t" << average_plaquette << std::endl;

            ++number_accepted;
        } else {
            std::cout << "Rejected.\n";
            links = old_links;

            ofs_energy_reject << number_computed << "\t"
                              << (new_energy / links.get_volume()) << std::endl;
            ofs_plaquette_reject << number_computed << "\t" << average_plaquette
                                 << std::endl;
        }

        if (do_write_config &&
            (chain_skip == 0 || number_computed % chain_skip == 0)) {
            std::string filename = (config_filename_format %
           number_stored).str();
            links.save(filename);
            ++number_stored;
        }

        ofs_accept << (accepted ? 1 : 0) << std::endl;

#ifdef OUTPUT
        std::cout << "Plaquette: " << average_plaquette << std::endl;
        std::cout << "Plaquette after global SU(2) transformation: "
                  << average_plaquette_2 << std::endl;
#endif

        auto const acceptance_rate = static_cast<double>(number_accepted) / number_computed;

        std::cout << "Acceptance rate: " << number_accepted << " / " << number_computed << " = "
                  << acceptance_rate << "\n"
                  << std::endl;
    }
}
