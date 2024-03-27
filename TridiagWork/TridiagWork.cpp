#include <iostream>
#include "TriDiagSolve.h"
using namespace TriDiagSolve;

int main()
{
    TriDiagSolver s{};
    std::cout << s.DefaultSolveNorm();
}
