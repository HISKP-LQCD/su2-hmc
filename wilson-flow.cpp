// Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

#include "wilson-flow.hpp"

#include "hybrid-monte-carlo.hpp"

Configuration exp(Configuration const &links) {
    Configuration result(links.length_space, links.length_time);
#pragma omp parallel for collapse(4)
    for (int n1 = 0; n1 < links.length_time; ++n1)
        for (int n2 = 0; n2 < links.length_space; ++n2)
            for (int n3 = 0; n3 < links.length_space; ++n3)
                for (int n4 = 0; n4 < links.length_space; ++n4)
                    for (int mu = 0; mu < 4; ++mu) {
                        result(n1, n2, n3, n4, mu) = links(n1, n2, n3, n4, mu).exp();
                    }

    return result;
}

Configuration flow_step(Configuration const &initial, double const time_step) {
    Configuration const &W0 = initial;
    Configuration Z0 = time_step * get_global_stout_exponential(W0);

    Configuration W1 = exp(1.0 / 4.0 * Z0) * W0;
    Configuration Z1 = time_step * get_global_stout_exponential(W1);

    Configuration W2 = exp(8.0 / 9.0 * Z1 - 17.0 / 36.0 * Z0) * W1;
    Configuration Z2 = time_step * get_global_stout_exponential(W2);

    Configuration W3 = exp(3.0 / 4.0 * Z2 - 8.0 / 9.0 * Z1 + 17.0 / 36.0 * Z0) * W2;

    return W3;
}

Configuration get_global_stout_exponential(Configuration const &links) {
    Configuration result(links.length_space, links.length_time);
#pragma omp parallel for collapse(4)
    for (int n1 = 0; n1 < links.length_time; ++n1)
        for (int n2 = 0; n2 < links.length_space; ++n2)
            for (int n3 = 0; n3 < links.length_space; ++n3)
                for (int n4 = 0; n4 < links.length_space; ++n4)
                    for (int mu = 0; mu < 4; ++mu) {
                        result(n1, n2, n3, n4, mu) =
                            - get_stout_exponential(n1, n2, n3, n4, mu, links);
                    }

    return result;
}
