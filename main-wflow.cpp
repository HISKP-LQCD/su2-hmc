// Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

#include "wilson-flow.hpp"

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

    Configuration flowed = initial;

    std::cout << "Flow step" << std::endl;
    for (int integration_step = 0;
         integration_step < integration_steps; ++integration_step) {
        std::cout << "  " << std::setw(5) << integration_step << std::endl;

        flowed = flow_step(flowed, time_step);
    }

    flowed.save(path_out);
}
