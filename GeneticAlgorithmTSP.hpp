#ifndef GENETIC_ALGORITHM_TSP_HPP
#define GENETIC_ALGORITHM_TSP_HPP
#include <vector>
#include "TSPInstance.hpp"
#include <algorithm>
#include <tuple>
#include <random>
#include <iterator>
#include <numeric>
#include <omp.h>
#include <thread>


class Parameters {
public:
    size_t pop_size;
    double mutation_rate;
    double crossover_rate;
    double heuristics_rate; // probability to apply local heuristics (delete-cross)
    double elitism;
    size_t max_iter;
    int random_seed;
    int tournament_size;
};

class ReturnInfo{
public:
    std::vector<size_t> best_solution;
    std::vector<double> best_fitness_per_iter;
    double best_fitness;
    double seconds_elapsed;
    size_t n_generations_untill_best;
};

class GeneticAlgorithmTSP {
public:
    GeneticAlgorithmTSP(const TSPInstance& instance, Parameters parameters);
    // Allow injecting one or more initial solutions to seed the population
    void setInjectedSolutions(const std::vector<std::vector<size_t>>& sols);
    ReturnInfo solve();

private:

    void generateInitialPop();

    size_t selection();

    void updatePopFitness();

    std::tuple<std::vector<size_t>, std::vector<size_t>>  crossover(size_t individual_a, size_t individual_b);

    void updatePopRank();

    void generateOffspring();

    void applyHeuristics();

    void deleteCross(size_t individual);

    void applyElitism();

    void applyMutation();

    void mutate(size_t individual);

    Parameters parameters;
    const TSPInstance& tsp_instance;

    std::vector<std::vector<size_t>> pop;
    std::vector<double> pop_fitness;
    std::vector<size_t> pop_rank;
    std::vector<std::vector<size_t>> pop_aux;
    std::vector<double> pop_aux_fitness;
    std::vector<size_t> pop_aux_rank;

    std::vector<size_t> global_best;
    std::vector<double> best_fitness_per_iter;
    double global_best_fitness;
    size_t best_iter;
    size_t current_iter;
    size_t n_elitist_sol;
    
    size_t n_cities;

    std::mt19937 random_generator;
    // Optional injected solutions to seed initial population
    std::vector<std::vector<size_t>> injected_solutions;
};
#endif