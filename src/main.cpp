#include <iostream>
#include "elements.hpp"
#include "mesh.hpp"
#include "mat_utilities.hpp"
#include "spectralFE.hpp"
#include "solver.hpp"

using namespace FETools;
using namespace Functions;


int main(int /*argc*/, char * /*argv*/[]){
    

    std::cout << "=================================================================" << std::endl;
    
    // Define some const variables to set some paramethers of the problems
    
    //Dimension of the problem
    constexpr unsigned DIM = 2;
    // Parameter to control propragation of f.p. error
    constexpr double tol= 1.e-9;
    // Degree of the FE space used by the solver
    constexpr unsigned r = 2;

    if constexpr(DIM != 2)
    {
        std::cout << "The solver has been built for DIM == 2 solutions up till now. Aborting...\n" << std::endl;
        return 1; 
    }
    // And call the solver for the considered problem

    // serialSolver solver(r);

    // solver.setup(1, "0.1000");
    // solver.assemble();
    // solver.solve();
    // solver.process("0.1000");

    // solver._errorH1();

    serialSolver solver(r);

    solver.convergence();


    
    return 0;
}