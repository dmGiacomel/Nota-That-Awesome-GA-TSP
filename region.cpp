#include "TSPInstance.hpp"
#include <vector>
#include <algorithm>
#include <cmath>

struct Regiao{
    int id_regiao;
    // Bounds (for visual/debug reference only)
    double x_min, x_max, y_min, y_max;
    std::vector<size_t> cidades_na_regiao;
};

/**
 * @brief Calculates cut values (thresholds) to split a vector into 'k' equal parts
 */
std::vector<double> calcular_cortes(std::vector<double>& valores, int k){
    if (valores.empty() || k <= 1) return{};
    
    // Sorts to find percentiles
    std::sort(valores.begin(), valores.end());
    
    std::vector<double> cortes;
    size_t n = valores.size();
    
    for (int i = 1; i < k; ++i){
        // Index of the percentile (1/3, 2/3, etc)
        size_t idx = (n * i) / k;
        cortes.push_back(valores[idx]);
    }
    return cortes;
}

/**
 * @brief Split the instance based on the distribution of points (quantiles)
 * Ensures load balancing across regions
 */
std::vector<Regiao> dividir_instancia_balanceada(const TSPInstance& instancia, int grid_rows, int grid_cols){
    
    size_t n = instancia.getTourSize();
    if (n == 0) return{};

    // Extract all X and Y coordinates
    std::vector<double> todos_x;
    std::vector<double> todos_y;
    todos_x.reserve(n);
    todos_y.reserve(n);

    for (size_t i = 0; i < n; ++i){
        auto [x, y] = instancia.getCityCoordinates(i);
        todos_x.push_back(x);
        todos_y.push_back(y);
    }

    // Calculate cuts (thresholds) for X and Y
    // for 2 columns, returns [median_x]
    // for 3 columns, returns [tertile_1_x, tertile_2_x]
    std::vector<double> cortes_x = calcular_cortes(todos_x, grid_cols);
    std::vector<double> cortes_y = calcular_cortes(todos_y, grid_rows);

    // Initialize regions
    std::vector<Regiao> regioes;
    regioes.resize(grid_rows * grid_cols);

    // Set IDs and approximate bounds (optional, for visualization)
    for (int r = 0; r < grid_rows; ++r){
        for (int c = 0; c < grid_cols; ++c){
            int idx = r * grid_cols + c;
            regioes[idx].id_regiao = idx;
            
            // Define bounds based on cuts (informational only)
            regioes[idx].x_min = (c == 0) ? -1e9 : cortes_x[c-1];
            regioes[idx].x_max = (c == grid_cols - 1) ? 1e9 : cortes_x[c];
            regioes[idx].y_min = (r == 0) ? -1e9 : cortes_y[r-1];
            regioes[idx].y_max = (r == grid_rows - 1) ? 1e9 : cortes_y[r];
        }
    }

    // Distribute the cities
    for (size_t i = 0; i < n; ++i){
      
        auto [x, y] = instancia.getCityCoordinates(i);

        // Find the column based on X cuts
        int col = 0;
        while (col < cortes_x.size() && x >= cortes_x[col]){
            col++;
        }

        // Find the row based on Y cuts
        int row = 0;
        while (row < cortes_y.size() && y >= cortes_y[row]){
            row++;
        }

        // Add to the correct region
        int idx = row * grid_cols + col;
        regioes[idx].cidades_na_regiao.push_back(i);
    }

    return regioes;
}