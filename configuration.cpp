// Copyright © 2016-2017 Martin Ueding <dev@martin-ueding.de>

#include "configuration.hpp"

#include <fstream>
#include <iostream>

Configuration::Configuration(int const length_space, int const length_time)
    : length_space(length_space),
      length_time(length_time),
      spacing_n4(4),
      spacing_n3(length_space * spacing_n4),
      spacing_n2(length_space * spacing_n3),
      spacing_n1(length_space * spacing_n2),
      volume(length_space * length_space * length_space * length_time),
      data(volume * 4) {}

void Configuration::save(std::string const &path) const {
    std::ofstream os(path, std::ios::out | std::ios::binary);
    os.write(reinterpret_cast<char const *>(data.data()), storage_size());
}

void Configuration::load(std::string const &path) {
    std::ifstream ifs(path, std::ios::in | std::ios::binary);
    ifs.read(reinterpret_cast<char *>(data.data()), storage_size());
}

void global_gauge_transformation(Matrix const &transformation, Configuration &links) {
    Matrix const adjoint = transformation.adjoint();
    for (int i = 0; i < links.get_size(); ++i) {
        links[i] = transformation * links[i] * adjoint;
    }
}

Configuration operator+(Configuration const &left, Configuration const &right) {
    assert(left.length_space == right.length_space);
    assert(left.length_time == right.length_time);

    Configuration result(right.length_space, right.length_time);
#pragma omp parallel for collapse(4)
    for (int n1 = 0; n1 < right.length_time; ++n1)
        for (int n2 = 0; n2 < right.length_space; ++n2)
            for (int n3 = 0; n3 < right.length_space; ++n3)
                for (int n4 = 0; n4 < right.length_space; ++n4)
                    for (int mu = 0; mu < 4; ++mu) {
                        result(n1, n2, n3, n4, mu) =
                            left(n1, n2, n3, n4, mu) + right(n1, n2, n3, n4, mu);
                    }

    return result;
}

Configuration operator-(Configuration const &left, Configuration const &right) {
    assert(left.length_space == right.length_space);
    assert(left.length_time == right.length_time);

    Configuration result(right.length_space, right.length_time);
#pragma omp parallel for collapse(4)
    for (int n1 = 0; n1 < right.length_time; ++n1)
        for (int n2 = 0; n2 < right.length_space; ++n2)
            for (int n3 = 0; n3 < right.length_space; ++n3)
                for (int n4 = 0; n4 < right.length_space; ++n4)
                    for (int mu = 0; mu < 4; ++mu) {
                        result(n1, n2, n3, n4, mu) =
                            left(n1, n2, n3, n4, mu) - right(n1, n2, n3, n4, mu);
                    }

    return result;
}

Configuration operator*(Configuration const &left, Configuration const &right) {
    assert(left.length_space == right.length_space);
    assert(left.length_time == right.length_time);

    Configuration result(right.length_space, right.length_time);
#pragma omp parallel for collapse(4)
    for (int n1 = 0; n1 < right.length_time; ++n1)
        for (int n2 = 0; n2 < right.length_space; ++n2)
            for (int n3 = 0; n3 < right.length_space; ++n3)
                for (int n4 = 0; n4 < right.length_space; ++n4)
                    for (int mu = 0; mu < 4; ++mu) {
                        result(n1, n2, n3, n4, mu) =
                            left(n1, n2, n3, n4, mu) * right(n1, n2, n3, n4, mu);
                    }

    return result;
}

Configuration operator*(double const left, Configuration const &right) {
    Configuration result(right.length_space, right.length_time);
#pragma omp parallel for collapse(4)
    for (int n1 = 0; n1 < right.length_time; ++n1)
        for (int n2 = 0; n2 < right.length_space; ++n2)
            for (int n3 = 0; n3 < right.length_space; ++n3)
                for (int n4 = 0; n4 < right.length_space; ++n4)
                    for (int mu = 0; mu < 4; ++mu) {
                        result(n1, n2, n3, n4, mu) = left * right(n1, n2, n3, n4, mu);
                    }

    return result;
}

Configuration make_hot_start(int const length_space,
                             int const length_time,
                             double const std,
                             int const seed) {
    std::mt19937 engine(seed);
    std::normal_distribution<double> dist(0, std);

    Configuration links(length_space, length_time);
    randomize_group(links, engine, dist);
    return links;
}

void randomize_algebra(Configuration &config,
                       std::mt19937 &engine,
                       std::normal_distribution<double> &dist) {
    for (int i = 0; i < config.get_size(); ++i) {
        Matrix const next = random_from_algebra(engine, dist);
        config[i] = next;
    }
}

void randomize_group(Configuration &config,
                     std::mt19937 &engine,
                     std::normal_distribution<double> &dist) {
    for (int i = 0; i < config.get_size(); ++i) {
        Matrix const next = random_from_group(engine, dist);
        config[i] = next;
    }
}
