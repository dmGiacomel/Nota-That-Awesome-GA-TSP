#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <iomanip>
#include <map>
#include <chrono> 
#include "TSPInstance.hpp"
#include "GeneticAlgorithmTSP.hpp"

//estrutura da entrada
//./exec_parte_2.x 
// <arquivo/problema> 
// <pop_size_1> 
// <mutation_rate_1> 
// <crossover_rate_1> 
// <heuristics_rate_1> 
// <elitism_rate_1> 
// <max_iter_1> 
// <tournament_size_1>
// <pop_size_2> 
// <mutation_rate_2> 
// <crossover_rate_2> 
// <heuristics_rate_2> 
// <elitism_rate_2> 
// <max_iter_2> 
// <tournament_size_2>
//  grid_rows
//  grid_columns
// <seed>

struct Region{
    int region_id;
    // original IDs from the TSP file
    std::vector<size_t> global_city_ids;
};

// Compute cut points to divide data into k equal parts (quantiles)
std::vector<double> compute_cuts(std::vector<double>& values, int k){
    if(values.empty() || k <= 1){
        return {};
    }
    std::sort(values.begin(), values.end());
    std::vector<double> cuts;
    size_t n = values.size();

    for(int i = 1; i < k; ++i){
        size_t idx = (n * i) / k;
        cuts.push_back(values[idx]);
    }
    return cuts;
}

// Divide the instance into balanced regions using quantile cuts on x and y
std::vector<Region> split_instance_balanced(const TSPInstance& instance, int grid_rows, int grid_cols){
    size_t n = instance.getTourSize();
    if(n == 0){
        return {};
    }
    std::vector<double> all_x, all_y;
    all_x.reserve(n);
    all_y.reserve(n);

    for(size_t i = 0; i < n; ++i){
        auto [x, y] = instance.getCityCoordinates(i);
        all_x.push_back(x);
        all_y.push_back(y);
    }

    std::vector<double> cuts_x = compute_cuts(all_x, grid_cols);
    std::vector<double> cuts_y = compute_cuts(all_y, grid_rows);

    std::vector<Region> regions(grid_rows * grid_cols);

    for(int i = 0; i < (int)regions.size(); ++i){
        regions[i].region_id = i;
    }

    for(size_t i = 0; i < n; ++i){

        auto [x, y] = instance.getCityCoordinates(i);

        int col = 0;
        while(col < (int)cuts_x.size() && x >= cuts_x[col]){
            col++;
        }

        int row = 0;
        while(row < (int)cuts_y.size() && y >= cuts_y[row]){
            row++;
        }

        int idx = row * grid_cols + col;
        regions[idx].global_city_ids.push_back(i);
    }
    return regions;
}

// Create a smaller TSPInstance containing only the region's cities
// Returns the new instance and a mapping LocalID -> GlobalID
std::pair<TSPInstance, std::vector<size_t>> create_subinstance(const TSPInstance& global_inst, const std::vector<size_t>& region_cities){
    size_t n_sub = region_cities.size();
    
    // local coordinates vector
    std::vector<std::tuple<double, double>> sub_coords(n_sub);
    // local adjacency matrix (flattened)
    std::vector<double> sub_adj(n_sub * n_sub);

    for(size_t i = 0; i < n_sub; ++i){

        size_t global_id_i = region_cities[i];
        sub_coords[i] = global_inst.getCityCoordinates(global_id_i);
        // fill diagonal 
        sub_adj[i * n_sub + i] = 0.0;
        for(size_t j = 0; j < n_sub; ++j){

            size_t global_id_j = region_cities[j];
            // copy the edge cost from the global matrix
            // access the global edge cost via public API getEdgeCost
            double cost = global_inst.getEdgeCost(global_id_i, global_id_j);
            sub_adj[i * n_sub + j] = cost;
        }
    }

    return { TSPInstance(sub_adj, sub_coords, n_sub), region_cities };
}

std::vector<std::vector<size_t>> generate_diversified_seeds(
    const std::vector<std::vector<size_t>>& regional_solutions,
    size_t n_global,
    size_t n_seeds,
    std::mt19937& rng)
{
    std::vector<std::vector<size_t>> seeds(n_seeds);

    auto n_regions = regional_solutions.size();
    for(size_t s = 0; s < n_seeds; s++){
        auto regional_sol_copy = regional_solutions;
        for(auto &region : regional_sol_copy){
            std::rotate(
                region.begin(),
                region.begin() + (rng() % region.size()),
                region.end()
            );
        }

        std::shuffle(regional_sol_copy.begin(), regional_sol_copy.end(), rng);

        for(auto &region : regional_sol_copy){
            seeds[s].insert(
                seeds[s].end(),
                region.begin(),
                region.end()
            );
        }
    }
    return seeds;
}

int main(int argc, char** argv){

    //inicia o rand generator com a seed passada por argv
    std::mt19937 random_generator(std::atoi(argv[18]));

    std::string arquivo_tsp = argv[1];
    Parameters parameters_local{
        std::stoull(argv[2]),     //pop
        std::stod(argv[3]),     //taxa de mutatacao
        std::stod(argv[4]),     //taxa de crossover
        std::stod(argv[5]),     //taxa de heuristica
        std::stod(argv[6]),     //taxa de elitismo
        std::stoull(argv[7]),     //iteracoes maximas
        0,                           //random seed
        std::atoi(argv[8])      //tamanho do torneio
    };

    Parameters parameters_global{
        std::stoull(argv[9]),     //pop
        std::stod(argv[10]),     //taxa de mutatacao
        std::stod(argv[11]),     //taxa de crossover
        std::stod(argv[12]),     //taxa de heuristica
        std::stod(argv[13]),     //taxa de elitismo
        std::stoull(argv[14]),     //iteracoes maximas
        0,                          //random seed
        std::atoi(argv[15])      //tamanho do torneio
    };

    int grid_rows = std::atoi(argv[16]);
    int grid_columns = std::atoi(argv[17]);

    TSPInstance tsp_instance = TSPInstance::fromTSPLIB(arquivo_tsp);
    
    auto t_start_local = std::chrono::steady_clock::now();

    auto regions = split_instance_balanced(tsp_instance, grid_rows, grid_columns);
    std::vector<std::vector<size_t>> assembled_global_solution;

    for(const auto& reg : regions){
        if(reg.global_city_ids.empty()){
            continue;
        }
        auto [tsp_local, local_to_global] = create_subinstance(tsp_instance, reg.global_city_ids);

        parameters_local.random_seed = random_generator();
        GeneticAlgorithmTSP ga_local(tsp_local, parameters_local);
        ReturnInfo info_local = ga_local.solve();

        // Translates local route (0..k) to global IDs and adds to solution
        auto solution_size = info_local.best_solution.size();
        for(size_t i{0}; i < solution_size; i++){
            info_local.best_solution[i] = local_to_global[info_local.best_solution[i]];
        }
        assembled_global_solution.emplace_back(info_local.best_solution);
    }

    auto t_end_local = std::chrono::steady_clock::now();
    auto t_total_local = std::chrono::duration<double>(t_end_local - t_start_local).count();

    parameters_global.random_seed = random_generator();
    GeneticAlgorithmTSP ga(tsp_instance, parameters_global);

    auto seeds = generate_diversified_seeds(
        assembled_global_solution,
        tsp_instance.getTourSize(),
        size_t(parameters_global.pop_size * 0.8), //20% da pop inicial continua aleatória
        random_generator
    );

    ga.setInjectedSolutions(seeds);
    ReturnInfo return_info = ga.solve();

    std::cout << return_info.best_fitness << ",";
    std::cout << return_info.seconds_elapsed << ","; // t_sol_global
    std::cout << t_total_local << ",";               //t_sol_local
    std::cout << return_info.n_generations_untill_best << std::endl;
}