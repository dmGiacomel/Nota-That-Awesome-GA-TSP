#include <iostream>
#include <vector>
#include "TSPInstance.hpp"
#include "Utils.hpp"
//arquivo de teste gerado pelo chat gpt
//compile: g++ -std=c++20 test_tsp_instance.cpp TSPInstance.cpp -O2 -o test

void test_flattened_coordinates() {
    std::cout << "Testing Utils::flattenedMatrixCoordinates...\n";

    bool ok = true;

    if (Utils::flattenedMatrixCoordinates(0, 0, 3) != 0) ok = false;
    if (Utils::flattenedMatrixCoordinates(0, 1, 3) != 1) ok = false;
    if (Utils::flattenedMatrixCoordinates(1, 0, 3) != 3) ok = false;
    if (Utils::flattenedMatrixCoordinates(2, 2, 3) != 8) ok = false;

    if (ok) std::cout << "OK.\n";
    else    std::cout << "FAILED.\n";
}

void test_random_instance() {
    std::cout << "Testing random TSPInstance generation...\n";

    TSPInstance tsp(4, 1.0, 2.0, 123);

    bool ok = true;

    // diagonal must be infinity
    for (int i = 0; i < 4; i++) {
        double d = tsp.getEdgeCost(i, i);
        if (!std::isinf(d)) ok = false;
    }

    // symmetric
    for (int i = 0; i < 4; i++) {
        for (int j = i+1; j < 4; j++) {
            if (tsp.getEdgeCost(i,j) != tsp.getEdgeCost(j,i))
                ok = false;
        }
    }

    if (ok) std::cout << "OK.\n";
    else    std::cout << "FAILED.\n";
}

void test_path_cost() {
    std::cout << "Testing path cost...\n";

    // a fixed 4x4 matrix:
    // 0  1  2  3
    // 1  0  4  5
    // 2  4  0  6
    // 3  5  6  0

    std::vector<double> adj = {
        // row 0
        std::numeric_limits<double>::infinity(), 1, 2, 3,
        // row 1
        1, std::numeric_limits<double>::infinity(), 4, 5,
        // row 2
        2, 4, std::numeric_limits<double>::infinity(), 6,
        // row 3
        3, 5, 6, std::numeric_limits<double>::infinity()
    };

    TSPInstance tsp(adj, 4);

    std::vector<int> path = {0, 1, 3, 2};

    // cost:
    // 0→1 = 1
    // 1→3 = 5
    // 3→2 = 6
    // 2→0 = 2
    // total = 14

    double cost = tsp.getPathCost(path);

    if (cost == 14.0) std::cout << "OK.\n";
    else std::cout << "FAILED: got " << cost << "\n";
}

int main() {
    test_flattened_coordinates();
    test_random_instance();
    test_path_cost();
    return 0;
}
