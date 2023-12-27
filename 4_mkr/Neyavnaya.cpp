#include <cstdio>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include "inmost.h"
#include <cmath>
using namespace INMOST;


double source(double x, double y){
    return 0;
}

double u(double x, double y){
    if ((x - 0.5) * (x-0.5) + (y - 0.5) *(y-0.5) <= (1 / 8) *(1/8)) return 1;
    else return 0;
}

double g(double x, double y){
    return 0;
}

int main(int argc, char**argv){
    Solver::Initialize(&argc, &argv);
    int n = 500, N = n * n;
    double h = 1. / n, dx = 10.0, dy = 1.0;
    double delta = 0.01;
    std::vector<double> full_sol(N);

    Sparse::Matrix A;
    Sparse::Vector b;

    A.Load("./matrix.mtx");
    b.Load("./rhs.mtx");

    Sparse::Vector ans = b;

    Solver Sol(Solver::INNER_ILU2);
    Sol.SetParameter("absolute_tolerance", "1e-12");
    Sol.SetParameter("verbosity", "2");
    Sol.SetParameter("relative_tolerance", "1e-7");

    for (int i = 0; i < 10; i++) {
        Sol.SetMatrix(A);
        bool solved = Sol.Solve(b, ans);

        for (int j = 0; j < n; j++) {
            A[j][j] += 1/delta;
            b[j] += ans[j]/delta;
        }

        if(!solved){
            std::cout << "Linear solver failure!" << std::endl;
            std::cout << "Reason: " << Sol.ReturnReason() << std::endl;
        }
        ans.Save("./solution{i}.mtx");
        std::cout << "Number of iterations: " << Sol.Iterations() << '\n';
        std::cout << "Residual: " << Sol.Residual() << '\n';
    }

    Solver::Finalize();
}
