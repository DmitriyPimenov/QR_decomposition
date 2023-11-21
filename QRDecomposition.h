#ifndef QRDECOMPOSITION_QRDECOMPOSITION_H
#define QRDECOMPOSITION_QRDECOMPOSITION_H

double residualTraceNorm(std::vector<std::vector<double>>& A, std::vector<double>& eigenvalues);

double residualLengthNorm(std::vector<std::vector<double>>& A, std::vector<double>& eigenvalues);

bool getEigenvalues(size_t n, std::vector<std::vector<double>>& A,
                    std::vector<double>& eigenvalues, double EPS);

void printResult(const std::vector<double>& x, size_t m);

#endif //QRDECOMPOSITION_QRDECOMPOSITION_H