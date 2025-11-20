#include "GeneticAlgorithmTSP.hpp"


GeneticAlgorithmTSP::GeneticAlgorithmTSP(const TSPInstance& instance, Parameters parameters)
    :parameters(parameters), 
    tsp_instance(instance),
    n_elitist_sol(parameters.elitism * parameters.pop_size), 
    best_iter(0),
    current_iter(0),
    global_best_fitness(0),
    random_generator(parameters.random_seed)
{

    pop = std::vector<std::vector<size_t>>(
        this->parameters.pop_size, 
        std::vector<size_t>(this->tsp_instance.getTourSize())
    );

    pop_aux = std::vector<std::vector<size_t>>(
        this->parameters.pop_size, 
        std::vector<size_t>(this->tsp_instance.getTourSize())
    );
}