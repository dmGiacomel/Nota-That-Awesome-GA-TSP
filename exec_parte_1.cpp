#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <iomanip>
#include <map>
#include "TSPInstance.hpp"
#include "GeneticAlgorithmTSP.hpp"

namespace fs = std::filesystem;

//compilar: g++ TSPInstance.cpp GeneticAlgorithmTSP.cpp exec_parte1.cpp -O2 -o exec_parte_1.x

//estrutura da entrada
//./exec_parte_1.x <arquivo/problema> <pop_size> <mutation_rate> <crossover_rate> <heuristics_rate> <elitism_rate> <max_iter> <random_seed> <tournament_size>
int main(int argc, char** argv){

    std::string arquivo_tsp = argv[1];
    Parameters parameters{
        std::stoull(argv[2]),     //pop
        std::stod(argv[3]),     //taxa de mutatacao
        std::stod(argv[4]),     //taxa de crossover
        std::stod(argv[5]),     //taxa de heuristica
        std::stod(argv[6]),     //taxa de elitismo
        std::stoull(argv[7]),     //iteracoes maximas
        std::atoi(argv[8]),     //random seed
        std::atoi(argv[9])      //tamanho do torneio
    };

    TSPInstance tsp_instance = TSPInstance::fromTSPLIB(arquivo_tsp);
    GeneticAlgorithmTSP ga(tsp_instance, parameters);
    ReturnInfo return_info = ga.solve();

    std::cout << return_info.best_fitness << ",";
    std::cout << return_info.seconds_elapsed << ",";
    std::cout << return_info.n_generations_untill_best << std::endl;
}