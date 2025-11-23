#include "TSPInstance.hpp"
#include "Utils.hpp"

double TSPInstance::getPathCost(const std::vector<size_t>& path) const{
    double path_cost{0.0};
    for (size_t i{0}; i < n_cities - 1; i++){
        path_cost += adj_matrix[Utils::flat2DIdx(path[i],path[i+1],n_cities)];
    }
    path_cost += adj_matrix[Utils::flat2DIdx(path[n_cities-1],path[0],n_cities)];
    return path_cost;
}

double TSPInstance::getEdgeCost(size_t city_a, size_t city_b) const{
    return adj_matrix[Utils::flat2DIdx(city_a, city_b, n_cities)];
}

//assumes city exists
std::tuple<double,double> TSPInstance::getCityCoordinates(size_t city) const{
    return coordinates[city];
}

size_t TSPInstance::getTourSize() const{
    return n_cities;
}

TSPInstance::TSPInstance(const std::vector<double>& adj_matrix, const std::vector<std::tuple<double, double>>& coordinates, size_t n_cities)
    :adj_matrix(adj_matrix),coordinates(coordinates),n_cities(n_cities)
{
}

//dummy TSP instance
TSPInstance::TSPInstance(){
}

//adj_matrix should be in flattened form
TSPInstance::TSPInstance(const std::vector<double>& adj_matrix, size_t n_cities)
    :adj_matrix(adj_matrix), n_cities(n_cities)
{
}

//generates random symetric TSP instance given the integer seed provided
//edge wheights range is [a,b)
//doesn't do anything special with coordinates
TSPInstance::TSPInstance(size_t n_cities, double a, double b, int random_seed)
    :   n_cities(n_cities), 
        adj_matrix(std::vector<double>(n_cities * n_cities, std::numeric_limits<double>::infinity())),
        coordinates(std::vector<std::tuple<double, double>>(n_cities))
{
    //standart mersenne_twister_generator
    //https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution.html
    std::mt19937 random_generator(random_seed);
    std::uniform_real_distribution<> distribution(a, b);

    //randomly filling edge weights but the main diagonal (no self-loops allowed)
    for (size_t i = 0; i < n_cities; i++){
        for (size_t j = i + 1; j < n_cities; j++){
            auto random = distribution(random_generator);
            adj_matrix[Utils::flat2DIdx(i, j, n_cities)] = random;
            adj_matrix[Utils::flat2DIdx(j, i, n_cities)] = random;
        }
    }
}