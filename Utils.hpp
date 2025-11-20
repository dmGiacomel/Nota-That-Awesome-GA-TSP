#ifndef UTILS_HPP
#define UTILS_HPP 
#include <cstddef>

namespace Utils{
    inline size_t flat2DIdx(size_t row, size_t column, size_t row_size){
        return row * row_size + column;
    }

    template <typename datatype>
    inline datatype** allocateMatrix(size_t rows, size_t columns) {

        datatype** m = new datatype*[rows];
        m[0] = new datatype[rows * columns];

        for (size_t i = 1; i < rows; i++) {
            m[i] = m[0] + i * columns;
        }

        return m;
    }

    template <typename datatype>
    inline void freeMatrix(datatype** matrix) {

        if (matrix == nullptr)
            return;

        delete[] matrix[0];
        delete[] matrix;
    }


}

#endif