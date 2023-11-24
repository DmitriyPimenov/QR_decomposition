#include <vector>
#include <iostream>
#include <cmath>
#include "QRDecomposition.h"

void leftMultiplication(std::vector<std::vector<double>>& A, const double cosPhi, const double sinPhi,
                        const size_t i, const size_t j, const size_t k) {
    double xi = A[i][k], xj = A[j][k];
    A[i][k] = xi * cosPhi - xj * sinPhi;
    A[j][k] = xi * sinPhi + xj * cosPhi;
}

void rightMultiplication(std::vector<std::vector<double>>& A, const double cosPhi, const double sinPhi,
                        const size_t i, const size_t j, const size_t k) {
    double xi = A[k][i], xj = A[k][j];
    A[k][i] = xi * cosPhi - xj * sinPhi;
    A[k][j] = xi * sinPhi + xj * cosPhi;
}

// reducing the matrix to an almost upper triangular form using the rotation method
void toAlmTrnForm(size_t n, std::vector<std::vector<double>>& A) {
    for (size_t j = 0; j < n - 2; ++j) {
        for (size_t i = j + 2; i < n; ++i) {
            double x = A[j + 1][j], y = A[i][j];
            if (std::abs(x) + std::abs(y) < 1e-10) {
                continue;
            }
            double cosPhi = x / sqrt(x*x + y*y); // Lemma 2 p.43
            double sinPhi = -y / sqrt(x*x + y*y); // Lemma 2 p.43
            for (size_t k = j; k < n; ++k) { // Lemma 5 p.45
                leftMultiplication(A, cosPhi, sinPhi, j + 1, i, k); // Lemma 4 p.45
            }
            for (size_t k = 0; k < n; ++k) {
                rightMultiplication(A, cosPhi, sinPhi, j + 1, i, k);
            }
        }
    }
}

// QR matrix decomposition using the reflection method
void QRbyReflMtd(std::vector<std::vector<double>>& A, size_t n, std::vector<std::vector<double>>& U) {
    for (size_t k = 0; k < n; ++k) {
        double sk = 0.;
        if (k < n - 1) {
            sk = A[k + 1][k] * A[k + 1][k]; // (15) p.118
        }
        std::vector<double> x(2, 0.);
        x[0] = A[k][k] - sqrt(A[k][k] * A[k][k] + sk); // (17) p.118
        if (k < n - 1) {
            x[1] = A[k + 1][k]; // (17) p.118
        }
        double xNorm = sqrt(x[0] * x[0] + sk); // (18) p.119
        if (xNorm < 1e-10) { // if x coordinates are too small, in particular if x = (0, 0)
            x[0] = A[k][k] < 0 ? -1. : 1.;
        }
        else { // (19) p.119
            x[0] /= xNorm;
            x[1] /= xNorm;
        }

        U[k] = x;
        if (k < n - 1) { // 3. p.119-120
            double u11 = 1 - 2*x[0]*x[0], u12 = -2*x[0]*x[1], u22 = 1 - 2*x[1]*x[1];
            for (size_t i = k; i < n; ++i) {
                double A1 = A[k][i], A2 = A[k + 1][i];
                A[k][i] = u11 * A1 + u12 * A2;
                A[k + 1][i] = u12 * A1 + u22 * A2;
            }
        }
        else {
            A[k][k] *= x[0];
        }
    }
}

void prodRU(size_t n, std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& U) {
    for (size_t j = 0; j < n - 1; ++j) {
        for (size_t i = 0; i <= j + 1; ++i) {
            double Rij = A[i][j], Rij1 = A[i][j + 1];
            double u11 = 1 - 2*U[j][0]*U[j][0], u12 = -2*U[j][0]*U[j][1], u22 = 1 - 2*U[j][1]*U[j][1];
            A[i][j] = Rij * u11 + Rij1 * u12;
            A[i][j + 1] = Rij * u12 + Rij1 * u22;
        }
    }
    for (size_t i = 0; i < n; ++i) {
        A[i][n - 1] *= U[n - 1][0];
    }
}

double matrixError(size_t n, std::vector<std::vector<double>>& A) {
    double error = 0;
    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            error = std::max(error, std::abs(A[i][j]));
        }
    }
    return error;
}

bool getEigenvalues(size_t n, std::vector<std::vector<double>>& A, std::vector<double>& eigenvalues, const double EPS) {
    std::vector<std::vector<double>> U(n, std::vector<double>(2));

    toAlmTrnForm(n, A);
    while (matrixError(n, A) > EPS) {
        QRbyReflMtd(A, n, U);
        prodRU(n, A, U);
    }

    for (size_t i = 0; i < n; ++i) {
        eigenvalues[i] = A[i][i];
    }

    return true;
}

double residualTraceNorm(std::vector<std::vector<double>>& A, std::vector<double>& eigenvalues) {
    double result = 0;
    for (size_t i = 0; i < A.size(); ++i) {
        result += A[i][i] - eigenvalues[i];
    }
    return std::abs(result);
}

double Norm2(std::vector<std::vector<double>>& A) {
    double l2A = 0;
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A.size(); ++j) {
            l2A += A[i][j] * A[i][j];
        }
    }
    l2A = sqrt(l2A);
    return l2A;
}

void printResult(const std::vector<double>& x, size_t m) {
    for (size_t i = 0; i < std::min(x.size(), m); ++i) {
        std::cout << x[i] << std::endl;
    }
}