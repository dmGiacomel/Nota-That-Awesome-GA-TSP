#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <iomanip>
#include <map>
#include "TSPInstance.hpp"
#include "GeneticAlgorithmTSP.hpp"

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
    // Flatten the concatenated regional solution once
    std::vector<size_t> base;
    base.reserve(n_global);
    for(const auto& reg : regional_solutions)
        for(size_t x : reg)
            base.push_back(x);

    std::vector<std::vector<size_t>> seeds;
    seeds.reserve(n_seeds);

    for(size_t s = 0; s < n_seeds; ++s){
        std::vector<size_t> tour = base;

        // 1. Random cyclic shift
        size_t shift = rng() % n_global;
        std::rotate(tour.begin(), tour.begin() + shift, tour.end());

        // 2. Random segment reversal (big 2-opt)
        size_t i = rng() % n_global;
        size_t j = rng() % n_global;
        if(i > j) std::swap(i, j);
        std::reverse(tour.begin() + i, tour.begin() + j + 1);

        // 3. Random inversion of region order
        // (simulate your "random permutation of sub-solutions")
        // We rebuild but using a shuffled order
        std::vector<size_t> idx(regional_solutions.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::shuffle(idx.begin(), idx.end(), rng);

        std::vector<size_t> permuted;
        permuted.reserve(n_global);
        for(size_t r : idx)
            permuted.insert(permuted.end(),
                            regional_solutions[r].begin(),
                            regional_solutions[r].end());

        // This overwrites with a valid permuted seed
        tour = std::move(permuted);
        seeds.push_back(std::move(tour));
    }
    return seeds;
}

namespace fs = std::filesystem;

int main(int argc, char** argv){

    std::mt19937 random_generator(42);

    // Part 2 configuration
    // Usage: tsp_regions.exe <tsp_file> [grid_rows grid_cols]
    std::string target_file = "./arquivos_tsp/brd14051.tsp"; // default test instance
    int GRID_ROWS = 4; // default 2x2 grid
    int GRID_COLS = 4;

    if(argc >= 2){
        target_file = argv[1];
    }
    if(argc >= 4){
        try{
            GRID_ROWS = std::stoi(argv[2]);
            GRID_COLS = std::stoi(argv[3]);
            if(GRID_ROWS <= 0) GRID_ROWS = 2;
            if(GRID_COLS <= 0) GRID_COLS = 2;
        } catch(...){
            std::cerr << "Invalido. Usando 2x2" << std::endl;
            GRID_ROWS = 2;
            GRID_COLS = 2;
        }
    }

    std::cout << "Abordagem Regional:\n";
    std::cout << "Lendo arquivo: " << target_file << std::endl;

    try{
        // Load global instance
        TSPInstance tsp_global = TSPInstance::fromTSPLIB(target_file);
        size_t n_global = tsp_global.getTourSize();
        std::cout << "Instancia global carregada. Cidades: " << n_global << "\n\n";

        // Divide to regions using quantiles
        std::cout << "Dividindo em grade " << GRID_ROWS << "x" << GRID_COLS << " (quantis)\n";
        auto regions = split_instance_balanced(tsp_global, GRID_ROWS, GRID_COLS);

        std::vector<std::vector<size_t>> assembled_global_solution;
        double sum_regional_costs = 0.0;

        // Resolves each region independently
        for(const auto& reg : regions){
            if(reg.global_city_ids.empty()){
                continue;
            }
            std::cout << "Processando regiao " << reg.region_id << " (" << reg.global_city_ids.size() << " cidades)";
            // Creates sub-instance (mapping local 0..k -> global ID)
            auto [tsp_local, local_to_global] = create_subinstance(tsp_global, reg.global_city_ids);

            // Configures GA for the region (fewer iterations due to smaller size)
            Parameters params_local;
            params_local.pop_size = 600;
            params_local.max_iter = 1000;
            params_local.mutation_rate = 0.05;
            params_local.crossover_rate = 0.7;
            params_local.heuristics_rate = 0.04; // delete-cross
            params_local.elitism = 0.01; 
            params_local.random_seed = random_generator();
            params_local.tournament_size = 3;

            // Executes local GA
            GeneticAlgorithmTSP ga_local(tsp_local, params_local);
            ReturnInfo info_local = ga_local.solve();

            std::cout << " Melhor custo local: " << info_local.best_fitness << std::endl;
            sum_regional_costs += info_local.best_fitness;

            // Translates local route (0..k) to global IDs and adds to solution
            auto solution_size = info_local.best_solution.size();
            for(size_t i{0}; i < solution_size; i++){
                info_local.best_solution[i] = local_to_global[info_local.best_solution[i]];
            }
            assembled_global_solution.emplace_back(info_local.best_solution);
        }

        std::cout << "\nSolucoes regionais concluidas\n";
        std::cout << "Soma dos custos regionais (sem conexao): " << sum_regional_costs << "\n";

        // Final global optimization
        std::cout << "\nIniciando montagem e refinamento global\n";

        Parameters params_global;
        params_global.pop_size = 600; 
        params_global.max_iter = 2000; // final refinement
        params_global.mutation_rate = 0.07;
        params_global.crossover_rate = 0.7;
        params_global.heuristics_rate = 0.05; // delete-cross
        params_global.elitism = 0.01;
        params_global.random_seed = random_generator();
        params_global.tournament_size = 3;

        GeneticAlgorithmTSP ga_global(tsp_global, params_global);

        auto seeds = generate_diversified_seeds(
            assembled_global_solution,
            n_global/2,
            params_global.pop_size, 
            random_generator
        );

        ga_global.setInjectedSolutions(seeds);

        // Execute global GA (the injected solution will seed the initial population)
        ReturnInfo info_global = ga_global.solve();

        std::cout << "Resultado global apos agrupamento das regioes:\n";
        std::cout << "Distancia: " << info_global.best_fitness << "\n";
        std::cout << "Tempo: " << info_global.seconds_elapsed << "s\n";
        std::cout << "Geracoes ate o melhor: " << info_global.n_generations_untill_best << "\n";

    } catch (const std::exception& e){
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}