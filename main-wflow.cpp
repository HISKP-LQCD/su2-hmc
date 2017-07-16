// Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

#include "wilson-flow.hpp"
#include "hybrid-monte-carlo.hpp"

#include <iomanip>
#include <iostream>
#include <string>

int main(int const argc, char const **argv) {

    if (argc != 7) {
        std::cout << "Please supply five arguments:\n"
                     "- Lattice extent in space [int]\n"
                     "- Lattice extent in time [int]\n"
                     "- Input path [string]\n"
                     "- Output path [string]\n"
                     "- Time step [double]\n"
                     "- Number of integration steps [int]"
                  << std::endl;
        std::exit(1);
    }

    int const length_space = std::stoi(argv[1]);
    int const length_time = std::stoi(argv[2]);
    std::string path_in(argv[3]);
    std::string path_out(argv[4]);
    double const time_step = std::stod(argv[5]);
    int const integration_steps = std::stoi(argv[6]);

    Configuration initial(length_space, length_time);
    initial.load(path_in);

    auto const plaquette_initial = get_plaquette_trace_average(initial);
    std::cout << "Initial plaquette: " << plaquette_initial.real() << std::endl;

    Configuration flowed = initial;

    int const column_width = 15;

    std::cout << std::setw(column_width) << "Flow step" << std::setw(column_width) << "Plaquette" << std::endl;
    for (int integration_step = 0;
         integration_step < integration_steps; ++integration_step) {


        flowed = flow_step(flowed, time_step);

        auto const plaquette = get_plaquette_trace_average(flowed);

        std::cout << std::setw(column_width) << integration_step
                  << std::setw(column_width) << plaquette.real() << std::endl;
    }

    flowed(0, 0, 0, 0, 0)(0, 0).real(1.0);
    flowed(0, 0, 0, 0, 0)(0, 0).imag(1.5);
    flowed(0, 0, 0, 0, 0)(1, 0).real(2.0);
    flowed(0, 0, 0, 0, 0)(1, 0).imag(2.5);
    flowed(0, 0, 0, 0, 0)(0, 1).real(3.0);
    flowed(0, 0, 0, 0, 0)(0, 1).imag(3.5);
    flowed(0, 0, 0, 0, 0)(1, 1).real(4.0);
    flowed(0, 0, 0, 0, 0)(1, 1).imag(4.5);

    flowed.save(path_out);
}
