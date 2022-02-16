#include "hw1.h"

using Matrix = std::vector<std::vector<double>>;

namespace algebra {
Matrix zeros(size_t n, size_t m)
{

    Matrix M(n, std::vector<double>(m));
    for (size_t i {}; i < n; i++) {
        for (size_t j {}; j < m; j++)
            M[i][j] = 0;
    }
    return M;
}

Matrix ones(size_t n, size_t m)
{
    Matrix M(n, std::vector<double>(m));
    for (size_t i {}; i < n; i++) {
        for (size_t j {}; j < m; j++)
            M[i][j] = 1;
    }
    return M;
}

Matrix random(size_t n, size_t m, double min, double max)
{
    if (min > max)
        throw std::logic_error("min should be smaller than max");
    else {
        // Efficient pseudo-random generator for random number in [min,max)
        std::random_device seeder;
        std::default_random_engine generator { seeder() };
        std::uniform_real_distribution<double> distribution(min, max);

        Matrix M(n, std::vector<double>(m));
        for (size_t i {}; i < n; i++) {
            for (size_t j {}; j < m; j++)
                M[i][j] = distribution(generator);
        }
        return M;
    }
}

void show(const Matrix& matrix)
{
    if (matrix.empty()) {
        throw std::logic_error("Matrix is empty");
    } else {
        double row { static_cast<double>(matrix.size()) };
        double col { static_cast<double>(matrix[0].size()) };
        for (size_t i {}; i < row; i++) {
            for (size_t j {}; j < col; j++)
                std::cout << std::setprecision(3) << std::fixed << std::setw(7) << matrix[i][j];
            std::cout << std::endl;
        }
    }
}

Matrix multiply(const Matrix& matrix, double c)
{
    if (matrix.empty()) {
        Matrix M {};
        return M;
    } else {
        double row { static_cast<double>(matrix.size()) };
        double col { static_cast<double>(matrix[0].size()) };
        Matrix M(row, std::vector<double>(col));
        M = matrix;
        for (size_t i {}; i < row; i++) {
            for (size_t j {}; j < col; j++)
                M[i][j] *= c;
        }
        return M;
    }
}

Matrix multiply(const Matrix& matrix1, const Matrix& matrix2)
{

    if (matrix1.empty() && matrix2.empty()) {
        Matrix M {};
        return M;
    } else {
        double row1 { static_cast<double>(matrix1.size()) };
        double col1 { static_cast<double>(matrix1[0].size()) };
        double row2 { static_cast<double>(matrix2.size()) };
        double col2 { static_cast<double>(matrix2[0].size()) };
        if (row2 == col1) {
            Matrix M(row1, std::vector<double>(col2));
            for (size_t i {}; i < row1; i++) {
                for (size_t j {}; j < col2; j++) {
                    for (size_t k {}; k < col1; k++)
                        M[i][j] += (matrix1[i][k] * matrix2[k][j]);
                }
            }
            return M;
        } else {
            throw std::logic_error("Matrices sizes don't match");
        }
    }
}

Matrix sum(const Matrix& matrix, double c)
{
    if (matrix.empty()) {
        Matrix M {};
        return M;
    } else {
        double row { static_cast<double>(matrix.size()) };
        double col { static_cast<double>(matrix[0].size()) };
        Matrix M(row, std::vector<double>(col));
        for (size_t i {}; i < row; i++) {
            for (size_t j {}; j < col; j++)
                M[i][j] = c + matrix[i][j];
        }
        return M;
    }
}

Matrix sum(const Matrix& matrix1, const Matrix& matrix2)
{
    if (matrix1.empty() && matrix2.empty()) {
        Matrix M {};
        return M;
    } else {
        double row1 { static_cast<double>(matrix1.size()) };
        double col1 { static_cast<double>(matrix1[0].size()) };
        double row2 { static_cast<double>(matrix2.size()) };
        double col2 { static_cast<double>(matrix2[0].size()) };
        if (row1 == row2 && col1 == col2) {
            Matrix M(row1, std::vector<double>(col1));
            for (size_t i {}; i < row1; i++) {
                for (size_t j {}; j < col1; j++)
                    M[i][j] = (matrix1[i][j] + matrix2[i][j]);
            }
            return M;
        } else {
            throw std::logic_error("Matrices sizes don't match");
        }
    }
}
}