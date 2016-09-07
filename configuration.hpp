// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#pragma once

#include <unsupported/Eigen/MatrixFunctions>

#include <cassert>
#include <vector>

class Configuration {
    using value_type = Eigen::Matrix2cd;

  public:
    Configuration(int const length_space, int const length_time);

    value_type &operator()(
        int const n1, int const n2, int const n3, int const n4, int const mu) {
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

        return data[index];
    }

    value_type &operator()(std::vector<int> const &coords, int const mu) {
        return operator()(coords[0], coords[1], coords[2], coords[3], mu);
    }

    int const length_space, length_time;

  private:
    int const spacing_n4, spacing_n3, spacing_n2, spacing_n1;
    int const volume;

    std::vector<value_type> data;
};
