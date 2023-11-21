#include <vector>
#include <iostream>
#include <cmath>
#include "QRDecomposition.h"

double residualTraceNorm(std::vector<std::vector<double>>& A, std::vector<double>& eigenvalues) {
    double result = 0;
    for (size_t i = 0; i < A.size(); ++i) {
        result += A[i][i] - eigenvalues[i];
    }
    return std::abs(result);
}

double residualLengthNorm(std::vector<std::vector<double>>& A, std::vector<double>& eigenvalues) {
    double l2A = 0, l2eigenvalues = 0;
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A.size(); ++j) {
            l2A += A[i][j] * A[i][j];
        }
        l2eigenvalues += eigenvalues[i] * eigenvalues[i];
    }
    l2A = sqrt(l2A);
    l2eigenvalues = sqrt(l2eigenvalues);
    return std::abs(l2A - l2eigenvalues);
}

bool getEigenvalues(size_t n, std::vector<std::vector<double>>& A,
                    std::vector<double>& eigenvalues, const double EPS) {
    
    return true;
}

void printResult(const std::vector<double>& x, size_t m) {
    for (size_t i = 0; i < std::min(x.size(), m); ++i) {
        std::cout << x[i] << std::endl;
    }
}