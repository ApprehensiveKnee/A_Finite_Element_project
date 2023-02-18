#include <iostream>
#include "elements.hpp"
#include "mesh.hpp"
#include "mat_utilities.hpp"
#include "quadrature.hpp"
#include "spectralFE.hpp"
#include "functions.hpp"
#include "solver.hpp"
//#include "psolver.hpp"

using namespace FETools;


int main(int /*argc*/, char * /*argv*/[]){
    

    std::cout << "=================================================================" << std::endl;

    constexpr unsigned int r = 2;

    // if constexpr(DIM != 2)
    // {
    //     std::cout << "The solver has been built for DIM == 2 solutions up till now. Aborting...\n" << std::endl;
    //     return 1; 
    // }
    // And call the solver for the considered problem


    serialSolver solver(r);

    // solver.setup(1, "0.2000");
    // solver.assemble();
    // solver.solve();
    // solver.process("0.2000");

    solver.convergence();



    
    return 0;
}