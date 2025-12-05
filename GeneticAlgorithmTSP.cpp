#include "GeneticAlgorithmTSP.hpp"
#include <chrono>

GeneticAlgorithmTSP::GeneticAlgorithmTSP(const TSPInstance& instance, Parameters parameters)
    :parameters(parameters),
    tsp_instance(instance),
    n_elitist_sol(static_cast<size_t>(parameters.elitism * parameters.pop_size)), 
    best_iter(0),
    current_iter(0),
    global_best_fitness(std::numeric_limits<double>::infinity()),
    n_cities(instance.getTourSize()),
    random_generator(parameters.random_seed),
    pop(std::vector<std::vector<size_t>>(
        parameters.pop_size, 
        std::vector<size_t>(instance.getTourSize())
    )),
    pop_aux(std::vector<std::vector<size_t>>(
        parameters.pop_size, 
        std::vector<size_t>(instance.getTourSize())
    )),
    best_fitness_per_iter(std::vector<double>(parameters.max_iter, 0.0)),
    pop_rank(parameters.pop_size),
    pop_aux_rank(parameters.pop_size),
    pop_fitness(parameters.pop_size, 0.0),
    pop_aux_fitness(parameters.pop_size, 0.0)
{
}

ReturnInfo GeneticAlgorithmTSP::solve(){

    auto t_start = std::chrono::high_resolution_clock::now();
    omp_set_dynamic(0);
    omp_set_num_threads(std::thread::hardware_concurrency());

    generateInitialPop();
    updatePopFitness();
    updatePopRank();

    best_fitness_per_iter[0] = pop_aux_fitness[pop_aux_rank[0]];
    global_best_fitness = best_fitness_per_iter[0];
    global_best = pop_aux[pop_aux_rank[0]];
    best_iter = 0; 

    std::swap(pop, pop_aux);
    std::swap(pop_rank, pop_aux_rank);
    std::swap(pop_fitness, pop_aux_fitness);

    for(current_iter = 1; current_iter < parameters.max_iter; current_iter++){
        applyElitism();
        generateOffspring();
        applyMutation();
        applyHeuristics();
        updatePopFitness();
        updatePopRank();

        double best_iter_fitness = pop_aux_fitness[pop_aux_rank[0]];
        best_fitness_per_iter[current_iter] = best_iter_fitness;
        if (best_iter_fitness < global_best_fitness){
            global_best_fitness = best_iter_fitness;
            best_iter = current_iter;
            global_best = pop_aux[pop_aux_rank[0]];
        }

        std::swap(pop, pop_aux);
        std::swap(pop_rank, pop_aux_rank);
        std::swap(pop_fitness, pop_aux_fitness);
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    ReturnInfo return_info = {
        global_best,
        best_fitness_per_iter,
        global_best_fitness,
        elapsed,
        best_iter
    };

    return return_info;
}

void GeneticAlgorithmTSP::generateInitialPop(){
    std::vector<size_t> raw_individual(n_cities);
    std::iota(raw_individual.begin(), raw_individual.end(), 0);
    // If user injected solutions, seed the first slots with them (if size matches)
    size_t injected = 0;
    if(!injected_solutions.empty()){
        for(const auto &sol : injected_solutions){
            if(injected >= parameters.pop_size) break;
            if(sol.size() != n_cities) continue; // only accept full-length solutions
            pop_aux[injected] = sol;
            ++injected;
        }
    }

    // Fill remaining population with random permutations
    for(size_t idx = injected; idx < parameters.pop_size; ++idx){
        auto &individual = pop_aux[idx];
        individual = raw_individual;
        std::shuffle(individual.begin(), individual.end(), random_generator);
    }
}

void GeneticAlgorithmTSP::setInjectedSolutions(const std::vector<std::vector<size_t>>& sols){
    injected_solutions.clear();
    for(const auto &s : sols){
        if(s.size() == n_cities) injected_solutions.push_back(s);
    }
}

size_t GeneticAlgorithmTSP::selection(){
    std::uniform_int_distribution<size_t> random_individual(0, parameters.pop_size - 1);
    size_t winner = random_individual(random_generator);
    double best_fitness = pop_fitness[winner];

    for(size_t i = 0; i < parameters.tournament_size - 1; i++){
        size_t candidate = random_individual(random_generator);
        if (pop_fitness[candidate] < best_fitness){
            winner = candidate;
            best_fitness = pop_fitness[candidate];
        }
    }

    return winner;
}

void GeneticAlgorithmTSP::applyMutation(){
    std::uniform_real_distribution will_mutate(0.0, 1.0);
    for (size_t i = 0; i < parameters.pop_size; i++){
        if(will_mutate(random_generator) < parameters.mutation_rate){
            mutate(i);
        }
    }
}

void GeneticAlgorithmTSP::applyHeuristics(){
    std::uniform_real_distribution will_apply(0.0, 1.0);
    for (size_t i = 0; i < parameters.pop_size; i++){
        if(will_apply(random_generator) < parameters.heuristics_rate){
            deleteCross(i);
        }
    }
}

void GeneticAlgorithmTSP::mutate(size_t individual){
    std::uniform_int_distribution<size_t> random_city(0, n_cities - 1);
    std::swap(pop_aux[random_city(random_generator)], pop_aux[random_city(random_generator)]);
}

void GeneticAlgorithmTSP::deleteCross(size_t individual){
    auto &tour = pop_aux[individual];

    for (size_t i{0}; i < n_cities - 3; i++){
        for (size_t j = 0; j < n_cities - 1; j++){
            if (j <= i + 1) continue;
            size_t a1 = tour[i];
            size_t b1 = tour[i+1];
            size_t a2 = tour[j];
            size_t b2 = tour[j+1];

            if (tsp_instance.segmentsIntersect(a1, b1, a2, b2)){
                size_t left = i + 1;
                size_t right = j;

                std::reverse(tour.begin() + left, tour.begin() + right + 1);
                return;
            }
        }
    }
}

void GeneticAlgorithmTSP::updatePopFitness(){
    #pragma omp parallel for schedule(static)
    for(size_t i = 0; i < parameters.pop_size; i++){
        pop_aux_fitness[i] = tsp_instance.getPathCost(pop_aux[i]);
    }
}

void GeneticAlgorithmTSP::updatePopRank(){
    std::iota(pop_aux_rank.begin(), pop_aux_rank.end(), 0);

    std::sort(pop_aux_rank.begin(),pop_aux_rank.end(),
        [&](size_t a, size_t b){
            return pop_aux_fitness[a] < pop_aux_fitness[b];
        }
    );
}

void GeneticAlgorithmTSP::generateOffspring(){

    std::uniform_real_distribution crossover_dist(0.0, 1.0);
    for (size_t i{n_elitist_sol}; i < parameters.pop_size; i += 2){
        size_t parent_a = selection();
        size_t parent_b;
        do
            parent_b = selection();
        while (parent_a == parent_b);
        
        if(crossover_dist(random_generator) < parameters.crossover_rate){
            auto [son_a, son_b] = crossover(parent_a, parent_b);
            pop_aux[i] = son_a;
            if (i + 1 < parameters.pop_size) 
                pop_aux[i + 1] = son_b;
        }else{
            pop_aux[i] = pop[parent_a];
            if (i + 1 < parameters.pop_size) 
                pop_aux[i+1] = pop[parent_b];
        }
    }
}

std::tuple<std::vector<size_t>, std::vector<size_t>> 
GeneticAlgorithmTSP::crossover(size_t parent_a, size_t parent_b){
    const std::vector<size_t>& A = pop[parent_a];
    const std::vector<size_t>& B = pop[parent_b];
    std::vector<size_t> son_a(n_cities);
    std::vector<size_t> son_b(n_cities);

    std::uniform_int_distribution<size_t> random_cut_dist(0, n_cities - 1);
    size_t first_cut = random_cut_dist(random_generator);
    size_t second_cut = random_cut_dist(random_generator);
    if (first_cut > second_cut){
        std::swap(first_cut, second_cut);
    }

    std::copy(
        A.begin() + first_cut, 
        A.begin() + second_cut + 1, 
        son_b.begin() + first_cut
    );
    std::copy(
        B.begin() + first_cut,
        B.begin() + second_cut + 1,
        son_a.begin() + first_cut
    );

    std::vector<char> already_in_a(n_cities, 0);
    std::vector<char> already_in_b(n_cities, 0);
    for(size_t i = first_cut; i <= second_cut; i++){
        already_in_a[son_a[i]] = 1;
        already_in_b[son_b[i]] = 1;
    }

    size_t idx_a = (second_cut + 1) % n_cities;
    size_t idx_b = idx_a;

    for (size_t i = 0; i < n_cities; i++){
        //filling son_a (has the segment from B)
        size_t city_in_A = A[(second_cut + 1 + i) % n_cities];
        if (!already_in_a[city_in_A]){
            son_a[idx_a] = city_in_A;
            idx_a = (idx_a + 1) % n_cities;
        }
        //filling son_b (has the segment from A)
        size_t city_in_B = B[(second_cut + 1 + i) % n_cities];
        if(!already_in_b[city_in_B]){
            son_b[idx_b] = city_in_B;
            idx_b = (idx_b + 1) % n_cities;
        }
    }

    return std::make_tuple(son_a, son_b);
}

void GeneticAlgorithmTSP::applyElitism(){
    for (size_t i{0}; i < n_elitist_sol; i++){
        pop_aux[i] = pop[pop_rank[i]];
    }
}