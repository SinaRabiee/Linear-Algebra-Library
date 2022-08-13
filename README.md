# AP HW1 (Linear algebra library in C++)
In this project, we implement different functions for working with vectors and matrices.

## zeros 
Creats a `n x m` matrix of zeros
```cpp
Matrix zeros(size_t n, size_t m)
```
## ones
Creats a `n x m` matrix of ones
```cpp
Matrix ones(size_t n, size_t m)
```
## random
Produces a `n x m` matrix of random numbers ranging from `min` to `max`
```cpp
Matrix random(size_t n, size_t m, double min, double max)
```
## show
Diplays the input of the function in a proper way
```cpp
void show(const Matrix& matrix)
```
## mutiply (scalar number and matrix)
Multiplies all the elements of `matrix1` with scalar number `c`
```cpp
Matrix multiply(const Matrix& matrix, double c)
```
## mutiply (two matrices)
Multiplies `matrix1` with `matrix2` in that order
```cpp
Matrix multiply(const Matrix& matrix1, const Matrix& matrix2)
```
## sum (scalar number and matrix)
Addes scalar number `c` to all the elements of `matrix1`
```cpp
Matrix sum(const Matrix& matrix, double c)
```
## sum (two matrices)
Addes `matrix1` to `matrix2`
```cpp
Matrix sum(const Matrix& matrix1, const Matrix& matrix2)
```
## transpose
Creates the transpose of the input matrix
```cpp
Matrix transpose(const Matrix& matrix)
```
## minor
Creates the minor of the input matrix with respect to nth row and mth column. A minor of a matrix is the same matrix with desired row and column removed.
```cpp
Matrix minor(const Matrix& matrix, size_t n, size_t m)
```
## determinant
Calculates the determinant of the input matrix
```cpp
double determinant(const Matrix& matrix)
```
## inverse
Calculates the inverse of the input matrix
```cpp
Matrix inverse(const Matrix& matrix)
```
## concatenate
Concatenates the input matrices along the specified axis (axis=0 alignes them vertically and axis=1 alignes them horizontally)
```cpp
Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis=0)
```
## ero_swap
Part of the elementary row operations. It swaps the `r1` row with `r2`.
```cpp
Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2)
```
## ero_multiply
Part of the elementary row operations. It multiplies the `r` row with constant `c`.
```cpp
Matrix ero_multiply(const Matrix& matrix, size_t r, double c)
```
## ero_sum
Part of the elementary row operations. It multiplies the `r1` row with constant `c` and adds it to `r2` row.
```cpp
Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2)
```
## upper_triangular
This functions uses elementary row operations to calculate the upper triangle form of the input matrix .
```cpp
Matrix upper_triangular(const Matrix& matrix)
```








