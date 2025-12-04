#include <iostream>
#include <vector>
#include "TSPInstance.hpp"
#include "Utils.hpp"
//arquivo de teste gerado pelo chat gpt
//compile: g++ -std=c++20 test_tsp_instance.cpp TSPInstance.cpp -O2 -o test

void test_flattened_coordinates() {
    std::cout << "Testing Utils::flattenedMatrixCoordinates...\n";

    bool ok = true;

    if (Utils::flat2DIdx(0, 0, 3) != 0) ok = false;
    if (Utils::flat2DIdx(0, 1, 3) != 1) ok = false;
    if (Utils::flat2DIdx(1, 0, 3) != 3) ok = false;
    if (Utils::flat2DIdx(2, 2, 3) != 8) ok = false;

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

void test_segments_intersect() {
    std::cout << "Testing TSPInstance::segmentsIntersect\n";

    // Create a small instance with coordinates for 6 cities
    size_t n = 6;
    std::vector<std::tuple<double,double>> coords(n);
    // 0: (0,0)
    coords[0] = std::make_tuple(0.0, 0.0);
    // 1: (1,1)
    coords[1] = std::make_tuple(1.0, 1.0);
    // 2: (0,1)
    coords[2] = std::make_tuple(0.0, 1.0);
    // 3: (1,0)
    coords[3] = std::make_tuple(1.0, 0.0);
    // 4: (2,0)
    coords[4] = std::make_tuple(2.0, 0.0);
    // 5: (3,0)
    coords[5] = std::make_tuple(3.0, 0.0);

    // adjacency matrix placeholder (not used by segmentsIntersect)
    std::vector<double> adj(n * n, std::numeric_limits<double>::infinity());

    TSPInstance tsp(adj, coords, n);

    bool ok = true;

    // Case 1: crossing segments: (0,1) vs (2,3) -> (0,0)-(1,1) crosses (0,1)-(1,0)
    if (!tsp.segmentsIntersect(0,1,2,3)) ok = false;

    // Case 2: non-crossing parallel: (0,3) vs (2,4) -> (0,0)-(1,0) and (0,1)-(2,0) do not cross
    if (tsp.segmentsIntersect(0,3,2,4)) ok = false;

    // Case 3: touching at endpoint: (0,3) vs (3,4) -> touch at city 3
    if (!tsp.segmentsIntersect(0,3,3,4)) ok = false;

    // Case 4: colinear overlapping: (3,5) vs (4,5) -> (1,0)-(3,0) and (2,0)-(3,0) overlap/touch
    if (!tsp.segmentsIntersect(3,5,4,5)) ok = false;

    if (ok) std::cout << "OK.\n";
    else    std::cout << "FAILED.\n";
}



int main() {
    test_flattened_coordinates();
    test_random_instance();
    test_path_cost();
    test_segments_intersect();
    return 0;
}
