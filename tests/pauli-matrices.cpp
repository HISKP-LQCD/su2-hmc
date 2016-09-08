// Copyright © 2016 Martin Ueding <dev@martin-ueding.de>

#include "../pauli-matrices.hpp"
#include "../sanity-checks.hpp"

#include <gtest/gtest.h>

class PauliMatricesTest : public ::testing::Test {
  protected:
    PauliMatricesTest() : pauli_matrices(PauliMatrices::get_instance()) {}

    PauliMatrices pauli_matrices;
};

TEST_F(PauliMatricesTest, squaring) {
    for (int i = 0; i < 3; ++i) {
        ASSERT_TRUE(is_unity(pauli_matrices.get(i) * pauli_matrices.get(i)));
    }
}

TEST_F(PauliMatricesTest, hermitian) {
    for (int i = 0; i < 3; ++i) {
        ASSERT_TRUE(is_hermitian(pauli_matrices.get(i)));
    }
}

TEST_F(PauliMatricesTest, traceless) {
    for (int i = 0; i < 3; ++i) {
        ASSERT_TRUE(is_traceless(pauli_matrices.get(i)));
    }
}
