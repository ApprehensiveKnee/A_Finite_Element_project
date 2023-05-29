// Preprocessor macro to check for out of bounds indexing
#define EIGEN_RUNTIME_NO_MALLOC
// Macro for the coloring algorithm


#include <iostream>
#include "variables.hpp"
#include "elements.hpp"
#include "mesh.hpp"
#include "mat_utilities.hpp"
#include "quadrature.hpp"
#include "spectralFE.hpp"
#include "functions.hpp"
#include "coloring.hpp"
#include "solver.hpp"
#include "psolver.hpp"
#include "psolver2.hpp"
#include "psolverColoring.hpp"


using namespace FETools;


int main(int /*argc*/, char * /*argv*/[]){
    

    std::cout << "=================================================================" << std::endl;

    // Define the degree of the FE space to solve the problem

    constexpr unsigned int r = 2;

    // And call the solver for the considered problem


    serialSolver solver(r);

    solver.setup("0.2000", true);


    solver.assemble();
    solver.solve();

    // ================ WRITE MATRIX TO A FILE =================
    // const static Eigen::IOFormat CSVFormat(StreamPrecision, DontAlignCols, ",", "\n");
    // Eigen::MatrixXd dense_mat1 = solver.getMat().toDense();
    // std::ofstream file1("system_mat_serial.csv");
    // file1 <<std::setprecision(5)<< dense_mat1.format(CSVFormat);

    // std::ofstream file01("system_rhs_serial.csv");
    // file01 <<std::setprecision(5)<< solver.getRHS().format(CSVFormat);
    // ================ WRITE MATRIX TO A FILE =================
    // std::cout << solver.getRHS() << std::endl;
    // std::cout << solver.getSol() << std::endl;
    solver.process("0.2000", true);

    // std::cout << "\n\nProcess ended correctly\n"<< std::endl;

    // solver.convergence();

    // ---------------------------------------------------------------------------------------------

    // ---------------------------------------------------------------------------------------------
    

    // Initialize Eigen for parallel environment

    //Eigen::initParallel();
    
    // parallelSolver psolver(r);
    // psolver.setup( "0.0500", true);
    // psolver.assemble();
    // // ================ WRITE MATRIX TO A FILE =================
    // // Eigen::MatrixXd dense_mat2 = psolver.getMat().toDense();
    // // std::ofstream file2("system_mat_parallel.csv");
    // // file2 <<std::setprecision(5)<< dense_mat2.format(CSVFormat);

    // // std::ofstream file02("system_rhs_parallel.csv");
    // // file02 <<std::setprecision(5)<< psolver.getRHS().format(CSVFormat);
    // // ================ WRITE MATRIX TO A FILE =================
    // psolver.solve();
    // psolver.process("0.0500", true);

    // parallelSolver2 psolver2(r);
    // std::cout << psolver2.getMat() << std::endl;
    // psolver2.setup(1, "0.5000");
    // psolver2.assemble();
    // std::cout << " COMMA " << std::endl;
    // psolver2.getMat().eval();
    // std::cout << psolver2.getMat().coeff(0,0) << std::endl;

    // ---------------------------------------------------------------------------------------------

    // ---------------------------------------------------------------------------------------------

    // parallelSolverColoring psolverColoring(r);
    // psolverColoring.setup("0.2000", true);
    // psolverColoring.assemble();
    // // ================ WRITE MATRIX TO A FILE =================
    // Eigen::MatrixXd dense_mat2 = psolverColoring.getMat().toDense();
    // std::ofstream file2("system_mat_parallel.csv");
    // file2 <<std::setprecision(5)<< dense_mat2.format(CSVFormat);

    // std::ofstream file02("system_rhs_parallel.csv");
    // file02 <<std::setprecision(5)<< psolverColoring.getRHS().format(CSVFormat);
    // // ================ WRITE MATRIX TO A FILE =================
    // psolverColoring.solve();
    // psolverColoring.process("0.2000", true);

    



    
    return 0;
}