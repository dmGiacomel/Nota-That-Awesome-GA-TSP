#ifndef GENETIC_ALGORITHM_TSP_HPP
#define GENETIC_ALGORITHM_TSP_HPP
#include <vector>
#include "TSPInstance.hpp"
#include <algorithm>
#include <random>
#include <iterator>

class Parameters {
public:
    size_t pop_size;
    double mutation_rate;
    double crossover_rate;
    double elitism;
    size_t max_iter;
    int random_seed;
};

class ReturnInfo{
public:
    std::vector<size_t> best_solution;
    double best_fitness;
    double seconds_elapsed;
    double n_generations_untill_best;
};

class GeneticAlgorithmTSP {
public:
    GeneticAlgorithmTSP(const TSPInstance& instance, Parameters parameters);
    ReturnInfo solve();

private:

    void generateInitialPop();

    //returns index of selected individual based on current rank
    size_t selection();

    void updatePopFitness();

    void updatePopRank();

    void updatePopAuxFitness();

    void updatePopAuxRank();

    size_t generateOffspring();

    Parameters parameters;
    const TSPInstance& tsp_instance;

    std::vector<std::vector<size_t>> pop;
    std::vector<double> pop_fitness;
    std::vector<size_t> pop_rank;
    std::vector<std::vector<size_t>> pop_aux;
    std::vector<double> pop_aux_fitness;
    std::vector<double> pop_aux_rank;

    std::vector<size_t> global_best;
    double global_best_fitness;
    size_t best_iter;
    size_t current_iter;
    size_t n_elitist_sol;

    std::mt19937 random_generator;
};
#endif