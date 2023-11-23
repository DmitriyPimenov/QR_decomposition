#include <iostream>
#include <vector>
#include <ctime>
#include "dataCreation.h"
#include "QRDecomposition.h"

int main(int argc, char* argv[]) {
    std::cout.precision(10);

    if (argc <= 1 || argc >= 4) {
        std::cerr << "Invalid number of arguments." << std::endl;
        exit(1);
    }

    size_t n;
    n = strtoul(argv[1], nullptr, 10);
    const size_t m = n;
    const double EPS = 1e-7;
    if (n == 0) {
        std::cerr << "The first argument is not a positive integer." << std::endl;
        exit(1);
    }

    std::vector<std::vector<double>> A(n, std::vector<double> (n));
    std::vector<double> eigenvalues(n);
    if (argc == 2) {
        createMatrixByFormula(A);
    }
    else {
        readMatrix(argv[2], A);
    }

    size_t startTime = clock();
    getEigenvalues(n, A, eigenvalues, EPS);
    size_t endTime = clock();

    std::cout << "p.6 Print result:" << std::endl;
    printResult(eigenvalues, m);
    std::cout << "p.7.1 Error trace: " << residualTraceNorm(A, eigenvalues) << std::endl;
    std::cout << "p.7.2 Error length: " << residualLengthNorm(A, eigenvalues) << std::endl;
    std::cout << "p.9 System solution time: " << endTime - startTime << std::endl;

    return 0;
}