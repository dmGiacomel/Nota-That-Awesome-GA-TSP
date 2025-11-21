#include "TSPInstance.hpp"
#include "Utils.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

double TSPInstance::getPathCost(const std::vector<int>& path) const{
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

TSPInstance TSPInstance::fromTSPLIB(const std::string& path){
    std::ifstream ifs(path);
    if(!ifs.is_open()){
        throw std::runtime_error("unable to open file: " + path);
    }

    std::string line;
    size_t dimension = 0;
    bool has_node_coord = false;
    bool is_euc2d = false;
    bool has_edge_weight_section = false;
    std::string edgeWeightType;

    // read header
    while(std::getline(ifs, line)){
        // trim
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        if(line.empty()) continue;
        if(line.rfind("DIMENSION", 0) == 0){
            auto pos = line.find(":");
            if(pos != std::string::npos){
                std::stringstream ss(line.substr(pos+1));
                ss >> dimension;
            }
        } else if(line.rfind("EDGE_WEIGHT_TYPE", 0) == 0){
            auto pos = line.find(":");
            if(pos != std::string::npos){
                std::stringstream ss(line.substr(pos+1));
                ss >> edgeWeightType;
                if(edgeWeightType == "EUC_2D") is_euc2d = true;
            }
        } else if(line == "NODE_COORD_SECTION"){
            has_node_coord = true;
            break;
        } else if(line == "EDGE_WEIGHT_SECTION"){
            has_edge_weight_section = true;
            break;
        } else if(line == "EOF"){
            break;
        }
    }

    if(dimension == 0){
        throw std::runtime_error("DIMENSION not found or zero in TSPLIB file");
    }

    std::vector<std::tuple<double,double>> coords(dimension);
    std::vector<double> adj;

    if(has_node_coord){
        // read coords
        size_t read = 0;
        while(read < dimension && std::getline(ifs, line)){
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            if(line.empty() || line == "EOF") break;
            std::stringstream ss(line);
            int id;
            double x,y;
            ss >> id >> x >> y;
            if(id >= 1 && static_cast<size_t>(id) <= dimension){
                coords[id-1] = std::make_tuple(x,y);
                ++read;
            }
        }
        // require EUC_2D coordinates: if not declared as EUC_2D, reject
        if(!is_euc2d){
            throw std::runtime_error("TSPLIB file must declare EDGE_WEIGHT_TYPE: EUC_2D when providing NODE_COORD_SECTION");
        }

        // build adj matrix from euclidean distances
        adj.assign(dimension * dimension, 0.0);
        for(size_t i=0;i<dimension;i++){
            for(size_t j=0;j<dimension;j++){
                if(i==j){
                    adj[Utils::flat2DIdx(i,j,dimension)] = 0.0;
                } else {
                    double xi = std::get<0>(coords[i]);
                    double yi = std::get<1>(coords[i]);
                    double xj = std::get<0>(coords[j]);
                    double yj = std::get<1>(coords[j]);
                    double d = std::sqrt((xi-xj)*(xi-xj) + (yi-yj)*(yi-yj));
                    adj[Utils::flat2DIdx(i,j,dimension)] = d;
                }
            }
        }

    } else if(has_edge_weight_section){
        // read full matrix values (assume flattened row-major numbers)
        adj.assign(dimension * dimension, std::numeric_limits<double>::infinity());
        size_t count = 0;
        while(std::getline(ifs, line)){
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            if(line.empty() || line == "EOF") break;
            std::stringstream ss(line);
            double v;
            while(ss >> v){
                if(count >= adj.size()) break;
                adj[count++] = v;
            }
            if(count >= adj.size()) break;
        }
        if(count < adj.size()){
            throw std::runtime_error("EDGE_WEIGHT_SECTION does not contain enough values");
        }
    } else {
        throw std::runtime_error("TSPLIB format not supported: missing NODE_COORD_SECTION or EDGE_WEIGHT_SECTION");
    }

    // construct instance
    TSPInstance inst(adj, coords, dimension);
    return inst;
}

void TSPInstance::printAdjMatrix(std::ostream& os) const{
    size_t n = n_cities;
    if(n == 0){
        os << "Empty adjacency matrix\n";
        return;
    }
    os << std::fixed << std::setprecision(2);
    size_t limit = std::min<size_t>(n, 10);

    os << "      |";
    for(size_t j=0;j<limit;j++){
        os << std::setw(7) << (j+1) << " |";
    }
    if(limit < n) os << " ...";
    os << '\n';

    os << std::string(static_cast<size_t>(limit) * 10, '-') << (limit < n ? "----" : "") << '\n';

    for(size_t i=0;i<limit;i++){
        os << std::setw(5) << (i+1) << " |";
        for(size_t j=0;j<limit;j++){
            double v = adj_matrix[Utils::flat2DIdx(i,j,n)];
            os << std::setw(7) << v << " |";
        }
        if(limit < n) os << " ...";
        os << '\n';
    }
    if(limit < n){
        os << "Matrix too large, printed first " << limit << "x" << limit << " entries\n";
    }
}