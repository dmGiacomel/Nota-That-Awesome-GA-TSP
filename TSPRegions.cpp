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


namespace fs = std::filesystem;

int main(int argc, char** argv){

    // Part 2 configuration
    // Usage: tsp_regions.exe <tsp_file> [grid_rows grid_cols]
    std::string target_file = "./arquivos_tsp/d657.tsp"; // default test instance
    int GRID_ROWS = 2; // default 2x2 grid
    int GRID_COLS = 2;

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

        std::vector<size_t> assembled_global_solution;
        assembled_global_solution.reserve(n_global);

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
            params_local.pop_size = 1000;
            params_local.max_iter = 1000;
            params_local.mutation_rate = 0.05;
            params_local.crossover_rate = 0.7;
            params_local.heuristics_rate = 0.2; // delete-cross
            params_local.elitism = 0.01; 
            params_local.random_seed = 42 + reg.region_id;
            params_local.tournament_size = 3;


            // Executes local GA
            GeneticAlgorithmTSP ga_local(tsp_local, params_local);
            ReturnInfo info_local = ga_local.solve();

            std::cout << " Melhor custo local: " << info_local.best_fitness << std::endl;
            sum_regional_costs += info_local.best_fitness;

            // Translates local route (0..k) to global IDs and adds to solution
            for(size_t id_local : info_local.best_solution){
                assembled_global_solution.push_back(local_to_global[id_local]);
            }
        }

        std::cout << "\nSolucoes regionais concluidas\n";
        std::cout << "Soma dos custos regionais (sem conexao): " << sum_regional_costs << "\n";

        // Final global optimization
        // The 'assembled_global_solution' is a simple concatenation; we need to refine the inter-region links
        // We'll use it (and mutations of it) to seed the global GA.
        
        std::cout << "\nIniciando montagem e refinamento global\n";

        Parameters params_global;
        params_global.pop_size = 1000; 
        params_global.max_iter = 1000; // final refinement
        params_global.mutation_rate = 0.05;
        params_global.crossover_rate = 0.7;
        params_global.heuristics_rate = 0.2; // delete-cross
        params_global.elitism = 0.01; // 1 individual
        params_global.random_seed = 100;
        params_global.tournament_size = 3;

        GeneticAlgorithmTSP ga_global(tsp_global, params_global);

        // Injeta a solucao regional concatenada como semente para o AG global
        std::vector<std::vector<size_t>> injected;
        if(assembled_global_solution.size() == n_global){
            injected.push_back(assembled_global_solution);
        }
        ga_global.setInjectedSolutions(injected);

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