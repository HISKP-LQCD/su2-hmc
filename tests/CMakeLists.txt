add_executable(tests

    ../configuration.cpp
    ../hybrid-monte-carlo.cpp
    ../hybrid-monte-carlo.cpp
    ../pauli-matrices.cpp
    ../sanity-checks.cpp
    ../su2matrix.cpp
    hybrid-monte-carlo.cpp
    main.cpp
    sanity-checks.cpp

)

target_link_libraries(tests gtest)
