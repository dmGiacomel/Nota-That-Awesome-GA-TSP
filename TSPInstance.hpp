#ifndef TSP_INSTANCE_HPP
#define TSP_INSTANCE_HPP
#include <vector>
#include <tuple>
#include <cstddef>
#include <random>
#include <limits>

//only the symmetric case is considered so far
class TSPInstance{
public:

    //assumes the path is a hamiltonian cycle
    //path is given by an array of visited cities
    //the first city should not be repeated at the end of the array
    double getPathCost(const std::vector<int>& path) const;

    //assumes the two cities exist
    double getEdgeCost(size_t city_a, size_t city_b) const;

    std::tuple<double,double> getCityCoordinates(size_t city) const;

    size_t getTourSize() const;

    //dummy TSP instance
    TSPInstance();

    //adj_matrix should be in flattened form
    TSPInstance(const std::vector<double>& adj_matrix, size_t n_cities);

    //initialize with adj matrix and coordinates vector
    TSPInstance(const std::vector<double>& adj_matrix, const std::vector<std::tuple<double, double>>& coordinates, size_t n_cities);

    //generates random symetric TSP instance given the integer seed provided
    //edge wheights range is [a,b)
    TSPInstance(size_t n_cities, double a, double b, int random_seed);

    TSPInstance(const TSPInstance&) = default;
    TSPInstance(TSPInstance&&) noexcept = default;
    TSPInstance& operator=(const TSPInstance&) = default;
    TSPInstance& operator=(TSPInstance&&) noexcept = default;

private:
    std::vector<double> adj_matrix;
    std::vector<std::tuple<double, double>> coordinates;
    size_t n_cities;
};

#endif