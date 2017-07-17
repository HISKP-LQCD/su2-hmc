// Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

#include "configuration.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
    void check_config(Configuration &cfg) {
        for (int n1 = 0; n1 < cfg.length_time; ++n1)
            for (int n2 = 0; n2 < cfg.length_space; ++n2)
                for (int n3 = 0; n3 < cfg.length_space; ++n3)
                    for (int n4 = 0; n4 < cfg.length_space; ++n4) {
                        int coords[] = {n1, n2, n3, n4};
                        for (int mu = 0; mu < 4; ++mu) {
                            auto const actual = cfg(n1, n2, n3, n4, mu)(0, 0).real();
                            double const expected = coords[mu];
                            if (fabs(actual - expected) > 1e-10) {
                                std::cout << "cfg(t=" << n1 << ", x=" << n2
                                          << ", y=" << n3 << ", z=" << n4 << ", mu=" << mu
                                          << ") = " << actual << " but should be "
                                          << expected << std::endl;
                            }
                        }
                    }
    }

    void read_test(Configuration &cfg, std::string const &path) {
        cfg.load(path);
        check_config(cfg);
    }

    void write_test(Configuration &cfg, std::string const &path) {
        for (int n1 = 0; n1 < cfg.length_time; ++n1)
            for (int n2 = 0; n2 < cfg.length_space; ++n2)
                for (int n3 = 0; n3 < cfg.length_space; ++n3)
                    for (int n4 = 0; n4 < cfg.length_space; ++n4) {
                        int coords[] = {n1, n2, n3, n4};
                        for (int mu = 0; mu < 4; ++mu) {
                            double const expected = coords[mu];
                            cfg(n1, n2, n3, n4, mu)(0, 0).real(expected);
                        }
                    }

        check_config(cfg);
        cfg.save(path);
    }
}

int main(int argc, char **argv) {
    if (argc != 5) {
        std::cout << "Please supply four arguments:\n"
                     "- Lattice extent in space [int]\n"
                     "- Lattice extent in time [int]\n"
                     "- Path [string]\n"
                     "- Action [string, either 'read' or 'write']\n"
                  << std::endl;
        std::exit(1);
    }

    int const length_space = std::stoi(argv[1]);
    int const length_time = std::stoi(argv[2]);
    std::string path(argv[3]);
    std::string action(argv[4]);

    Configuration cfg(length_space, length_time);

    bool read;

    if (action == "read")
        read = true;
    else if (action == "write")
        read = false;
    else
        throw std::domain_error("The given action is not supported");

    if (read)
        read_test(cfg, path);
    else
        write_test(cfg, path);
}
