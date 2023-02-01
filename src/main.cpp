#include <iostream>
#include "elements.hpp"
#include "mesh.hpp"
#include "mat_utilities.hpp"
#include "solver.hpp"

using namespace FETools;
using namespace Functions;


int main(int /*argc*/, char * /*argv*/[]){

    
    

    std::cout << "===========================================================================================" << std::endl;
    
    if constexpr(DIM != 2)
    {
        std::cout << "The solver has been built for DIM == 2 solutions up till now. Aborting...\n" << std::endl;
        return 1; 
    }
    // And call the solver for the considered problem

    serialSolver solver(r);

    solver.setup(1);
    solver.assemble();
    solver.solve();
    solver.process();


    
    return 0;
}