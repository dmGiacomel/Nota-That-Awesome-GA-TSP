#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <map>
#include "TSPInstance.hpp"
#include "GeneticAlgorithmTSP.hpp"

// ==========================================
// 1. UTILITÁRIOS DE DIVISÃO (Quantis)
// ==========================================

struct Regiao {
    int id_regiao;
    std::vector<size_t> cidades_globais;
};

// Função auxiliar para calcular percentis
std::vector<double> calcular_cortes(std::vector<double>& valores, int k) {
    if (valores.empty() || k <= 1) return {};
    std::sort(valores.begin(), valores.end());
    std::vector<double> cortes;
    size_t n = valores.size();
    for (int i = 1; i < k; ++i) {
        size_t idx = (n * i) / k;
        cortes.push_back(valores[idx]);
    }
    return cortes;
}

// Divide a instância em regiões equilibradas
std::vector<Regiao> dividir_instancia(const TSPInstance& instancia, int grid_rows, int grid_cols) {
    size_t n = instancia.getTourSize();
    if (n == 0) return {};

    std::vector<double> todos_x, todos_y;
    todos_x.reserve(n);
    todos_y.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        auto [x, y] = instancia.getCityCoordinates(i);
        todos_x.push_back(x);
        todos_y.push_back(y);
    }

    std::vector<double> cortes_x = calcular_cortes(todos_x, grid_cols);
    std::vector<double> cortes_y = calcular_cortes(todos_y, grid_rows);

    std::vector<Regiao> regioes(grid_rows * grid_cols);
    for(int i=0; i < (int)regioes.size(); ++i) regioes[i].id_regiao = i;

    for (size_t i = 0; i < n; ++i) {
        auto [x, y] = instancia.getCityCoordinates(i);
        int col = 0;
        while (col < (int)cortes_x.size() && x >= cortes_x[col]) col++;
        int row = 0;
        while (row < (int)cortes_y.size() && y >= cortes_y[row]) row++;
        int idx = row * grid_cols + col;
        regioes[idx].cidades_globais.push_back(i);
    }
    return regioes;
}

// Cria uma sub-instância apenas com as cidades daquela região
std::pair<TSPInstance, std::vector<size_t>> criar_sub_instancia(const TSPInstance& global_inst, const std::vector<size_t>& cidades_regiao) {
    size_t n_sub = cidades_regiao.size();
    std::vector<std::tuple<double, double>> sub_coords(n_sub);
    std::vector<double> sub_adj(n_sub * n_sub);

    for(size_t i = 0; i < n_sub; ++i) {
        size_t global_id_i = cidades_regiao[i];
        sub_coords[i] = global_inst.getCityCoordinates(global_id_i);
        for(size_t j = 0; j < n_sub; ++j) {
            sub_adj[i * n_sub + j] = global_inst.getEdgeCost(global_id_i, cidades_regiao[j]);
        }
    }
    return {TSPInstance(sub_adj, sub_coords, n_sub), cidades_regiao};
}

// ==========================================
// 2. FUNÇÃO DE REFINAMENTO (Manual)
// ==========================================
// Como não podemos injetar a população no AG do seu amigo, usamos esta função
// para aplicar o 1-delete-cross na rota costurada final.
void refinar_rota_costurada(std::vector<size_t>& rota, const TSPInstance& tsp, int max_iteracoes) {
    size_t n = rota.size();
    for (int iter = 0; iter < max_iteracoes; ++iter) {
        bool melhorou = false;
        // Procura cruzamentos
        for (size_t i = 0; i < n - 1; ++i) {
            for (size_t j = i + 2; j < n; ++j) {
                if (i == 0 && j == n - 1) continue; // Adjacente no ciclo

                size_t u = rota[i];
                size_t v = rota[i+1];
                size_t x = rota[j];
                size_t y = rota[(j+1) % n];

                if (tsp.segmentsIntersect(u, v, x, y)) {
                    // Remove cruzamento (2-opt)
                    std::reverse(rota.begin() + i + 1, rota.begin() + j + 1);
                    melhorou = true;
                    // 1-delete-cross: para após o primeiro e reinicia o loop externo
                    goto proxima_iteracao; 
                }
            }
        }
        if (!melhorou) break; // Sem cruzamentos, sai
        proxima_iteracao:;
    }
}

// ==========================================
// 3. CENÁRIOS DE EXECUÇÃO
// ==========================================

// --- PARTE 1: AG GLOBAL (Cenário A) ---
ReturnInfo rodar_cenario_global(const TSPInstance& tsp) {
    std::cout << "\n>>> Executando Cenário A: AG Global..." << std::endl;
    
    Parameters params;
    params.pop_size = std::max((size_t)50, tsp.getTourSize()); 
    params.max_iter = 2000;
    params.mutation_rate = 0.01;
    params.crossover_rate = 0.9;
    params.heuristics_rate = 0.1; // 1-delete-cross
    params.elitism = 1.0 / params.pop_size;
    params.random_seed = 42;
    params.tournament_size = 5;

    GeneticAlgorithmTSP ga(tsp, params);
    return ga.solve();
}

// --- PARTE 2: REGIONAL + COSTURA (Cenário B) ---
double rodar_cenario_regional(const TSPInstance& tsp) {
    std::cout << "\n>>> Executando Cenário B: Regional + Juncao..." << std::endl;
    auto inicio = std::chrono::high_resolution_clock::now();

    // 1. Dividir
    int linhas = 2, colunas = 2; // Grade 2x2
    auto regioes = dividir_instancia(tsp, linhas, colunas);
    std::vector<size_t> rota_costurada;

    // 2. Resolver Regiões
    for (const auto& reg : regioes) {
        if (reg.cidades_globais.empty()) continue;
        
        auto [tsp_local, mapa_ids] = criar_sub_instancia(tsp, reg.cidades_globais);
        
        Parameters params;
        params.pop_size = std::max((size_t)20, tsp_local.getTourSize()); 
        params.max_iter = 500; // Menos iterações
        params.mutation_rate = 0.05;
        params.crossover_rate = 0.9;
        params.heuristics_rate = 0.2; 
        params.elitism = 0.1; 
        params.random_seed = 100 + reg.id_regiao;
        params.tournament_size = 3;

        GeneticAlgorithmTSP ga(tsp_local, params);
        ReturnInfo info = ga.solve();

        // Concatena na rota final (Costura simples)
        for (size_t id_local : info.best_solution) {
            rota_costurada.push_back(mapa_ids[id_local]);
        }
    }

    // 3. Otimizar a Junção (Rodar AG/Heurística na rota completa)
    // Como não podemos injetar no AG do amigo, rodamos o delete-cross manual aqui
    std::cout << "    Refinando a juncao das regioes..." << std::endl;
    refinar_rota_costurada(rota_costurada, tsp, 1000);

    auto fim = std::chrono::high_resolution_clock::now();
    double tempo = std::chrono::duration<double>(fim - inicio).count();
    
    double custo_final = tsp.getPathCost(rota_costurada);
    std::cout << "    Tempo Cenario B: " << tempo << "s" << std::endl;
    
    return custo_final;
}

// ==========================================
// MAIN
// ==========================================
int main() {
    std::string arquivo = "./arquivos_tsp/eil101.tsp"; 

    try {
        TSPInstance tsp = TSPInstance::fromTSPLIB(arquivo);
        std::cout << "Instancia: " << arquivo << " (" << tsp.getTourSize() << " cidades)\n";

        // Roda Parte 1
        ReturnInfo resultado_A = rodar_cenario_global(tsp);
        
        // Roda Parte 2
        double resultado_B_custo = rodar_cenario_regional(tsp);

        // Comparação Final
        std::cout << "\n============================================\n";
        std::cout << "             RELATORIO FINAL                \n";
        std::cout << "============================================\n";
        std::cout << "Resultado A (Global Puro): " << resultado_A.best_fitness << "\n";
        std::cout << "Resultado B (Regional):    " << resultado_B_custo << "\n";
        
        double diferenca = ((resultado_B_custo - resultado_A.best_fitness) / resultado_A.best_fitness) * 100.0;
        
        std::cout << "\nComparativo:\n";
        if (resultado_B_custo < resultado_A.best_fitness) {
            std::cout << ">> O metodo Regional venceu por " << std::fixed << std::setprecision(2) << -diferenca << "%!\n";
        } else {
            std::cout << ">> O metodo Global venceu por " << std::fixed << std::setprecision(2) << diferenca << "%.\n";
        }
        std::cout << "============================================\n";

    } catch (const std::exception& e) {
        std::cerr << "Erro: " << e.what() << std::endl;
    }
    return 0;
}