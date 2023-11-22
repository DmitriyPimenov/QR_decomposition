#include <vector>
#include <iostream>
#include <cmath>
#include "QRDecomposition.h"

void pr(std::vector<std::vector<double>>& A) { // Не забыть удалить
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A.size(); ++j)
            std::cout << (std::abs(A[i][j]) < 1e-10 ? 0. : A[i][j]) << ' ';
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void QR(std::vector<std::vector<double>>& A, size_t n, size_t ii) {
    std::vector<std::vector<double>> U(n, std::vector<double> (n, 0.));
    for (size_t i = 0; i < n; ++i)
        U[i][i] = 1.;
    for (size_t k = 0; k < n; ++k) {
        std::vector<double> a(n - k);
        for (size_t i = k; i < n; ++i)
            a[i - k] = A[i][k];

        double aNorm = 0;
        for (double i : a)
            aNorm += i * i;
        aNorm = sqrt(aNorm);

        std::vector<double> x = a;
        x[0] -= aNorm;
        double xNorm = 0;
        for (double i : x)
            xNorm += i * i;
        xNorm = sqrt(xNorm);
        for (double& i : x)
            i /= xNorm;
        if (xNorm < 1e-10) {
            x[0] = A[k][k] > 0 ? 1. : -1.;
            for (size_t i = 1; i < n - k; ++i) {
                x[i] = 0.;
            }
        }

        std::vector<std::vector<double>> u(n - k, std::vector<double> (n - k, 0.));
        for (size_t i = 0; i < n - k; ++i)
            u[i][i] = 1.;
        for (size_t i = 0; i < n - k; ++i) {
            for (size_t j = 0; j < n - k; ++j) {
                u[i][j] -= 2 * x[i] * x[j];
            }
        }
        std::vector<std::vector<double>> tmp(n - k, std::vector<double> (n - k, 0.));
        for (size_t i = 0; i < n - k; ++i) {
            for (size_t j = 0; j < n - k; ++j) {
                for (size_t t = 0; t < n - k; ++t) {
                    tmp[i][j] += u[i][t] * A[k + t][k + j];
                }
            }
        }
        for (size_t i = 0; i < n - k; ++i) {
            for (size_t j = 0; j < n - k; ++j) {
                A[k + i][k + j] = tmp[i][j];
            }
        }

        std::vector<std::vector<double>> ut(n, std::vector<double> (n, 0.));
        for (size_t i = 0; i < n; ++i) {
            ut[i][i] = 1.;
        }
        for (size_t i = 0; i < n - k; ++i) {
            for (size_t j = 0; j < n - k; ++j) {
                ut[k + i][k + j] = u[i][j];
            }
        }
        tmp = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                for (size_t t = 0; t < n; ++t) {
                    tmp[i][j] += ut[i][t] * U[t][j];
                }
            }
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                U[i][j] = tmp[i][j];
            }
        }
        if (ii == 0) {
            std::cout << "x" << std::endl;
            for (double i : x)
                std::cout << i << ' ';
            std::cout << std::endl;
            std::cout << "U" << std::endl;
            pr(U);
            std::cout << "A" << std::endl;
            pr(A);
        }
    }

    std::vector<std::vector<double>> R = A, Q(n, std::vector<double>(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            Q[i][j] = U[j][i]; // Q = U^t
        }
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            A[i][j] = 0.;
            for (size_t t = 0; t < n; ++t) {
                A[i][j] += R[i][t] * Q[t][j];
            }
        }
    }

    std::vector<std::vector<double>> checkA(n, std::vector<double> (n, 0.));
    std::vector<std::vector<double>> QU = checkA;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t t = 0; t < n; ++t) {
                checkA[i][j] += Q[i][t] * R[t][j];
                QU[i][j] += Q[i][t] * U[t][j];
            }
        }
    }
    if (ii == 0) {
        std::cout << "R" << std::endl;
        pr(R);
        std::cout << "Q" << std::endl;
        pr(Q);
        std::cout << "checkA" << std::endl;
        pr(checkA);
        std::cout << "QU" << std::endl;
        pr(QU);
        std::cout << "A" << std::endl;
        pr(A);
    }
}

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

bool getEigenvalues(size_t n, std::vector<std::vector<double>>& A,
                    std::vector<double>& eigenvalues, const double EPS) {
    pr(A);
    for (size_t j = 0; j < n - 2; ++j) {
        for (size_t i = j + 2; i < n; ++i) {
            double x = A[j + 1][j], y = A[i][j];
            if (std::abs(x) + std::abs(y) < EPS) {
                continue;
            }
            double cosPhi = x / sqrt(x*x + y*y); // Lemma 2
            double sinPhi = -y / sqrt(x*x + y*y); // Lemma 2
            for (size_t k = j; k < n; ++k) { // Lemma 5
                leftMultiplication(A, cosPhi, sinPhi, j + 1, i, k); // Lemma 4
            }
            for (size_t k = 0; k < n; ++k) {
                rightMultiplication(A, cosPhi, sinPhi, j + 1, i, k);
            }
        }
    }

    pr(A);
    for (size_t i = 0; i < 1000; ++i)
        QR(A, n, i);
    pr(A);

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

void printResult(const std::vector<double>& x, size_t m) {
    for (size_t i = 0; i < std::min(x.size(), m); ++i) {
        std::cout << x[i] << std::endl;
    }
}