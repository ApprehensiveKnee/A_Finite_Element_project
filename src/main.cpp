#include <iostream>
#include "elements.hpp"
#include "mesh.hpp"
#include "mat_utilities.hpp"
#include "quadrature.hpp"
#include "spectralFE.hpp"
#include "functions.hpp"
#include "solver.hpp"
#include "psolver.hpp"
#include "psolver2.hpp"

// Preprocessor macro to check for out of bounds indexing
#define EIGEN_RUNTIME_NO_MALLOC


using namespace FETools;


int main(int /*argc*/, char * /*argv*/[]){
    

    std::cout << "=================================================================" << std::endl;

    // Define the degree of the FE space to solve the problem

    constexpr unsigned int r = 4;

    // And call the solver for the considered problem


    serialSolver solver(r);

    solver.setup("0.2000", true);
    solver.assemble();
    solver.solve();
    // std::cout << solver.getMat() << std::endl;
    // std::cout << solver.getRHS() << std::endl;
    // std::cout << solver.getSol() << std::endl;
    solver.process("0.2000", true);

    std::cout << "\n\nProcess ended correctly\n"<< std::endl;

    // solver.convergence();

    // ---------------------------------------------------------------------------------------------

    // ---------------------------------------------------------------------------------------------
    

    // Initialize Eigen for parallel environment

    // Eigen::initParallel();
    
    // parallelSolver psolver(r);
    // psolver.setup(true, "0.5000");
    // psolver.assemble();
    // std::cout << psolver.getMat() << std::endl;
    // psolver.solve();

    // parallelSolver2 psolver2(r);
    // std::cout << psolver2.getMat() << std::endl;
    // psolver2.setup(1, "0.5000");
    // psolver2.assemble();
    // std::cout << " COMMA " << std::endl;
    // psolver2.getMat().eval();
    // std::cout << psolver2.getMat().coeff(0,0) << std::endl;

    



    
    return 0;
}