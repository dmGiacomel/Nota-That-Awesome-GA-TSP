#include <iostream>
#include <string>
#include <filesystem>
#include <random>
#include <vector>
#include <tuple>
#include <algorithm>
#include "TSPInstance.hpp"

namespace fs = std::filesystem;

int main(){
    std::string pasta_arquivos = "./arquivos_tsp/";

    std::cout << "Starting file processing" << std::endl;

    // Parameters for sampling segment pairs per instance
    const std::size_t MAX_SAMPLES = 1000; // cap number of random segment-pair tests
    const std::size_t MAX_EXAMPLES = 10; // how many example intersecting pairs to print

    try{
        for(const fs::directory_entry& entry : fs::directory_iterator(pasta_arquivos)){
            if(!entry.is_regular_file()) continue;
            if(entry.path().extension() != ".tsp") continue;

            std::string nome_arquivo = entry.path().string();
            std::cout << "Processing: " << nome_arquivo << std::endl;

            try{
                TSPInstance tsp = TSPInstance::fromTSPLIB(nome_arquivo);

                std::size_t n = tsp.getTourSize();
                if(n < 4){
                    std::cout << "  Skipping (need at least 4 distinct points)." << std::endl;
                    continue;
                }

                // Deterministic RNG seeded with filename length and n to keep results reproducible
                std::seed_seq seed{static_cast<unsigned int>(nome_arquivo.size()), static_cast<unsigned int>(n)};
                std::mt19937 rng(seed);
                std::uniform_int_distribution<std::size_t> dist(0, n-1);

                std::size_t samples = std::min<std::size_t>(MAX_SAMPLES, (n*(n-1)/2));
                std::size_t intersections_found = 0;
                std::vector<std::tuple<std::size_t,std::size_t,std::size_t,std::size_t>> examples;

                for(std::size_t s = 0; s < samples; ++s){
                    // pick two distinct segments AB and CD with all endpoints distinct
                    std::size_t a=0,b=0,c=0,d=0;
                    // try to sample distinct endpoints; limit attempts to avoid infinite loops
                    const int MAX_ATTEMPTS = 50;
                    int attempts = 0;
                    do{
                        a = dist(rng);
                        b = dist(rng);
                        c = dist(rng);
                        d = dist(rng);
                        ++attempts;
                    } while((a==b || c==d || a==c || a==d || b==c || b==d) && attempts < MAX_ATTEMPTS);

                    if(a==b || c==d || a==c || a==d || b==c || b==d){
                        // fallback: skip this sample if couldn't find distinct endpoints
                        continue;
                    }

                    bool inter = tsp.segmentsIntersect(a,b,c,d);
                    if(inter){
                        ++intersections_found;
                        if(examples.size() < MAX_EXAMPLES) examples.emplace_back(a,b,c,d);
                    }
                }

                std::cout << "  Points: " << n << " | Samples tested: " << samples << " | Intersections found: " << intersections_found << std::endl;
                if(!examples.empty()){
                    std::cout << "  Example intersections (a,b,c,d): ";
                    for(const auto &t : examples){
                        std::cout << "(" << std::get<0>(t) << "," << std::get<1>(t) << "," << std::get<2>(t) << "," << std::get<3>(t) << ") ";
                    }
                    std::cout << std::endl;
                }

                std::cout << "File processed successfully: " << nome_arquivo << std::endl;
            } catch(const std::exception& e){
                std::cerr << "Error processing " << nome_arquivo << ": " << e.what() << std::endl;
            }
        }
    } catch(const fs::filesystem_error& e){
        std::cerr << "Error: Unable to open directory: " << pasta_arquivos << std::endl;
        return 1;
    }

    std::cout << "All files processed" << std::endl;
    return 0;
}
