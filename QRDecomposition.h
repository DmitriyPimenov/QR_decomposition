#ifndef QRDECOMPOSITION_QRDECOMPOSITION_H
#define QRDECOMPOSITION_QRDECOMPOSITION_H

bool getEigenvalues(size_t n, std::vector<std::vector<double>>& A,
                    std::vector<double>& eigenvalues, double EPS);

double residualTraceNorm(std::vector<std::vector<double>>& A, std::vector<double>& eigenvalues);

double Norm2(std::vector<std::vector<double>>& A);

void printResult(const std::vector<double>& x, size_t m);

#endif //QRDECOMPOSITION_QRDECOMPOSITION_H