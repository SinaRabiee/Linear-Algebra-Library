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

Matrix transpose(const Matrix& matrix)
{
    if (matrix.empty()) {
        Matrix M {};
        return M;
    } else {
        double row { static_cast<double>(matrix.size()) };
        double col { static_cast<double>(matrix[0].size()) };
        Matrix M(col, std::vector<double>(row));
        for (size_t i {}; i < row; i++) {
            for (size_t j {}; j < col; j++)
                M[j][i] = matrix[i][j];
        }
        return M;
    }
}

Matrix minor(const Matrix& matrix, size_t n, size_t m)
{
    if (matrix.empty()) {
        throw std::logic_error("Matrix is empty");
    } else {
        double row { static_cast<double>(matrix.size()) };
        double col { static_cast<double>(matrix[0].size()) };
        if (n < row && m < col) {
            double minor_row {};
            double minor_col {};
            Matrix temp(row, std::vector<double>(col));
            Matrix M(row - 1, std::vector<double>(col - 1));

            for (size_t i {}; i < row; i++) {
                if (i != n) {
                    for (size_t j {}; j < col; j++) {
                        if (j != m) {
                            temp[minor_row][minor_col] = matrix[i][j];
                            minor_col++;
                        }
                    }
                    minor_col = 0;
                    minor_row++;
                }
            }
            for (size_t i {}; i < row - 1; i++) {
                for (size_t j {}; j < col - 1; j++)
                    M[i][j] = temp[i][j];
            }
            return M;
        } else
            throw std::logic_error("Invalid minor size");
    }
}

double determinant(const Matrix& matrix)
{

    if (matrix.empty())
        return 1.0;
    else {
        double row { static_cast<double>(matrix.size()) };
        double col { static_cast<double>(matrix[0].size()) };
        double det {};
        if (row == col) {
            if (row == 1)
                det = matrix[0][0];
            else if (row == 2)
                det = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
            else {
                for (size_t i {}; i < row; i++)
                    det += pow(-1, i) * matrix[0][i] * determinant(minor(matrix, 0, i));
            }
            return det;
        } else
            throw std::logic_error("Matrix is not square");
    }
}

Matrix inverse(const Matrix& matrix)
{
    if (matrix.empty()) {
        Matrix M {};
        return M;
    } else {
        double row { static_cast<double>(matrix.size()) };
        double col { static_cast<double>(matrix[0].size()) };
        if (row == col) {
            if (determinant(matrix) == 0)
                throw std::logic_error("Matrix is singular (Determinant = 0)");
            else {
                double det { determinant(matrix) };
                Matrix M(row, std::vector<double>(col));
                for (size_t i {}; i < row; i++) {
                    for (size_t j {}; j < col; j++)
                        M[i][j] = pow(-1, (i + j)) * determinant(minor(matrix, i, j));
                }
                M = multiply(transpose(M), 1 / det);
                return M;
            }
        } else
            throw std::logic_error("Matrix is not square");
    }
}

Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis)
{
    if (matrix1.empty() || matrix2.empty()) {
        if (matrix1.empty() && matrix2.empty()) {
            Matrix M {};
            return M;
        } else
            throw std::logic_error("An empty matrix is is not concatable with a non-empty matirix");
    } else {
        double row1 { static_cast<double>(matrix1.size()) };
        double col1 { static_cast<double>(matrix1[0].size()) };
        double row2 { static_cast<double>(matrix2.size()) };
        double col2 { static_cast<double>(matrix2[0].size()) };

        if (axis == 0) {
            if (col1 != col2) {
                throw std::logic_error("columns of matrices don't match");
            } else {
                Matrix M(row1 + row2, std::vector<double>(col1));
                for (size_t i {}; i < row1; i++) {
                    for (size_t j {}; j < col1; j++)
                        M[i][j] = matrix1[i][j];
                }
                for (size_t i {}; i < row2; i++) {
                    for (size_t j {}; j < col2; j++)
                        M[i + static_cast<size_t>(row1)][j] = matrix2[i][j];
                }
                return M;
            }
        } else if (axis == 1) {
            if (row1 != row2) {
                throw std::logic_error("rows of matrices don't match");
            } else {
                Matrix M(row1, std::vector<double>(col1 + col2));
                for (size_t i {}; i < row1; i++) {
                    for (size_t j {}; j < col1; j++)
                        M[i][j] = matrix1[i][j];
                }
                for (size_t i {}; i < row2; i++) {
                    for (size_t j {}; j < col2; j++)
                        M[i][j + static_cast<size_t>(col1)] = matrix2[i][j];
                }
                return M;
            }
        } else {
            throw std::logic_error("Inavalid axis value");
        }
    }
}

Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2)
{
    if (matrix.empty()) {
        throw std::logic_error("Matrix is empty");
    } else {
        double row { static_cast<double>(matrix.size()) };
        double col { static_cast<double>(matrix[0].size()) };
        if (r1 < row && r2 < row) {
            Matrix temp(1, std::vector<double>(col));
            Matrix M { matrix };
            for (size_t j {}; j < col; j++) {
                temp[0][j] = matrix[r2][j];
                M[r2][j] = matrix[r1][j];
                M[r1][j] = temp[0][j];
            }
            return M;
        } else
            throw std::logic_error("Invalid values for desired rows (Out of range)");
    }
}

Matrix ero_multiply(const Matrix& matrix, size_t r, double c)
{
    if (matrix.empty()) {
        throw std::logic_error("Matrix is empty");
    } else {
        double row { static_cast<double>(matrix.size()) };
        double col { static_cast<double>(matrix[0].size()) };
        if (r < row) {
            Matrix M { matrix };
            for (size_t j {}; j < col; j++)
                M[r][j] = c * matrix[r][j];
            return M;
        } else
            throw std::logic_error("Invalid values for desired row");
    }
}

Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2)
{
    if (matrix.empty()) {
        throw std::logic_error("Matrix is empty");
    } else {
        double row { static_cast<double>(matrix.size()) };
        double col { static_cast<double>(matrix[0].size()) };
        if (r1 < row && r2 < row) {
            Matrix temp(1, std::vector<double>(col));
            Matrix M { matrix };
            for (size_t j {}; j < col; j++)
                M[r2][j] += c * matrix[r1][j];
            return M;
        } else
            throw std::logic_error("Invalid values for desired rows");
    }
}

Matrix upper_triangular(const Matrix& matrix)
{

    if (matrix.empty()) {
        Matrix M {};
        return M;
    } else {
        double row { static_cast<double>(matrix.size()) };
        double col { static_cast<double>(matrix[0].size()) };
        if (row == col) {
            double temp {};
            Matrix M { matrix };
            /*for (size_t i {}; i < row - 1; i++)
                if (M[i][i] == 0)
                    M = algebra::ero_swap(M, i, i + 1);*/
            for (size_t i {}; i < row; i++) {
                if (M[i][i] == 0)
                    M = algebra::ero_swap(M, i, i + 1);
                for (size_t k { i + 1 }; k < row; k++) {
                    temp = -(M[k][i] / M[i][i]);
                    M = algebra::ero_sum(M, i, temp, k);
                }
            }
            return M;
        } else
            throw std::logic_error("Matrix is not square");
    }
}
}