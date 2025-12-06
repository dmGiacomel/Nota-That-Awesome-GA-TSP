#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <iomanip>
#include "TSPInstance.hpp"
#include <algorithm> // Para sort
#include <iomanip>   // Para setw

namespace fs = std::filesystem;

// Estrutura para representar a região
struct Regiao{
    int id_regiao;
    std::vector<size_t> cidades_na_regiao;
};

// Função auxiliar para calcular cortes (percentis)
std::vector<double> calcular_cortes(std::vector<double>& valores, int k){
    if (valores.empty() || k <= 1) return{};
    std::sort(valores.begin(), valores.end());
    std::vector<double> cortes;
    size_t n = valores.size();
    for (int i = 1; i < k; ++i){
        size_t idx = (n * i) / k;
        cortes.push_back(valores[idx]);
    }
    return cortes;
}

// Função de divisão balanceada (Quantis)
std::vector<Regiao> dividir_instancia_balanceada(const TSPInstance& instancia, int grid_rows, int grid_cols){
    size_t n = instancia.getTourSize();
    if (n == 0) return{};

    // Extrair todas as coordenadas
    std::vector<double> todos_x; 
    std::vector<double> todos_y;
    todos_x.reserve(n);
    todos_y.reserve(n);

    for (size_t i = 0; i < n; ++i){
        auto [x, y] = instancia.getCityCoordinates(i);
        todos_x.push_back(x);
        todos_y.push_back(y);
    }

    // Calcular cortes
    std::vector<double> cortes_x = calcular_cortes(todos_x, grid_cols);
    std::vector<double> cortes_y = calcular_cortes(todos_y, grid_rows);

    // nicializar regiões
    std::vector<Regiao> regioes(grid_rows * grid_cols);
    for(int i=0; i < regioes.size(); ++i) regioes[i].id_regiao = i;

    // Distribuir cidades
    for (size_t i = 0; i < n; ++i){
        auto [x, y] = instancia.getCityCoordinates(i);

        // Acha coluna baseada nos cortes de X
        int col = 0;
        while (col < cortes_x.size() && x >= cortes_x[col]) col++;

        // Acha linha baseada nos cortes de Y
        int row = 0;
        while (row < cortes_y.size() && y >= cortes_y[row]) row++;

        int idx = row * grid_cols + col;
        regioes[idx].cidades_na_regiao.push_back(i);
    }

    return regioes;
}

int main(){
    std::string pasta_arquivos = "./arquivos_tsp/";

    if (!fs::exists(pasta_arquivos)){
        std::cerr << "Erro: A pasta " << pasta_arquivos << " nao existe" << std::endl;
        return 1;
    }

    std::cout << "Iniciando processamento de instancias reais" << std::endl;

    for (const auto& entry : fs::directory_iterator(pasta_arquivos)){
        if (entry.is_regular_file() && entry.path().extension() == ".tsp"){
            
            std::string caminho = entry.path().string();
            std::string nome_arquivo = entry.path().filename().string();
            
            try{
                // Carrega a instância
                TSPInstance tsp = TSPInstance::fromTSPLIB(caminho);
                size_t n = tsp.getTourSize();

                std::cout << "\nArquivo: " << nome_arquivo << " | Cidades: " << n << std::endl;

                // Teste 1: Grade 2x2
                std::cout << "  -> Dividindo em 2x2 (4 regioes):" << std::endl;
                auto regioes_2x2 = dividir_instancia_balanceada(tsp, 2, 2);
                
                for (const auto& r : regioes_2x2){
                    std::cout << "     Regiao " << r.id_regiao << ": " 
                              << r.cidades_na_regiao.size() << " cidades." << std::endl;
                }

                // Teste 2: Grade 4x4 (Só se tiver cidades suficientes)
                if (n > 50){
                    std::cout << "  -> Dividindo em 4x4 (16 regioes):" << std::endl;
                    auto regioes_4x4 = dividir_instancia_balanceada(tsp, 4, 4);
                    
                    // Imprime estatisticas resumidas para nao poluir
                    size_t min_c = n, max_c = 0;
                    for (const auto& r : regioes_4x4){
                        if (r.cidades_na_regiao.size() < min_c) min_c = r.cidades_na_regiao.size();
                        if (r.cidades_na_regiao.size() > max_c) max_c = r.cidades_na_regiao.size();
                    }
                    std::cout << "     (Min cidades/regiao: " << min_c << ", Max: " << max_c << ")" << std::endl;
                }

            } catch (const std::exception& e){
                // Ignora arquivos que nao sejam EUC_2D ou estejam corrompidos
                // std::cerr << "Pulei " << nome_arquivo << ": " << e.what() << std::endl;
            }
        }
    }

    return 0;
}