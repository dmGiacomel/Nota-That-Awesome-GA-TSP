#include <iostream>
#include <string>
#include <filesystem>
#include "TSPInstance.hpp"

namespace fs = std::filesystem;

int main(){
    std::string pasta_arquivos = "./arquivos_tsp/";

    std::cout << "Starting file processing" << std::endl;

    try{
        for(const auto& entry : fs::directory_iterator(pasta_arquivos)){
            if(!entry.is_regular_file()) continue;
            if(entry.path().extension() != ".tsp") continue;

            std::string nome_arquivo = entry.path().string();
            std::cout << "Processing: " << nome_arquivo << std::endl;

            try{
                TSPInstance tsp = TSPInstance::fromTSPLIB(nome_arquivo);
                // Printing of adjacency matrix
                tsp.printAdjMatrix(std::cout);
                std::cout << "File loaded successfully: " << nome_arquivo << std::endl;
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
