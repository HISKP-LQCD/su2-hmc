// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#pragma once

#include <unsupported/Eigen/MatrixFunctions>

#include <cassert>
#include <vector>

class Configuration {
  public:
    using value_type = Eigen::Matrix2cd;

    Configuration(int const length_space, int const length_time);

    value_type &operator()(
        int const n1, int const n2, int const n3, int const n4, int const mu) {
        return data[get_index(n1, n2, n3, n4, mu)];
    }

    const value_type &operator()(
        int const n1, int const n2, int const n3, int const n4, int const mu) const {
        return data[get_index(n1, n2, n3, n4, mu)];
    }

    value_type &operator()(std::vector<int> const &coords, int const mu) {
        return data[get_index(coords[0], coords[1], coords[2], coords[3], mu)];
    }

    const value_type &operator()(std::vector<int> const &coords, int const mu) const {
        return data[get_index(coords[0], coords[1], coords[2], coords[3], mu)];
    }

    int length_space, length_time;

    size_t storage_size() const { return data.size() * sizeof(value_type); };

    void save(std::string const &path) const;

    int get_volume() const { return volume; }

  private:
    int spacing_n4, spacing_n3, spacing_n2, spacing_n1;
    int volume;

    std::vector<value_type> data;

    size_t get_index(int const n1, int const n2, int const n3, int const n4, int const mu) const {
        assert(-1 <= n1 && n1 <= length_time);
        assert(-1 <= n2 && n2 <= length_space);
        assert(-1 <= n3 && n3 <= length_space);
        assert(-1 <= n4 && n4 <= length_space);

        // Periodic boundary conditions.
        int const n1_p = (n1 + length_time) % length_time;
        int const n2_p = (n2 + length_space) % length_space;
        int const n3_p = (n3 + length_space) % length_space;
        int const n4_p = (n4 + length_space) % length_space;

        assert(0 <= n1_p && n1_p < length_time);
        assert(0 <= n2_p && n2_p < length_space);
        assert(0 <= n3_p && n3_p < length_space);
        assert(0 <= n4_p && n4_p < length_space);

        int const index = n1_p * spacing_n1 + n2_p * spacing_n2 +
                          n3_p * spacing_n3 + n4_p * spacing_n4 + mu;

        assert(0 <= index && index < volume);

        return index;
    }
};
