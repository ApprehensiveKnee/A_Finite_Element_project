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

#include <boost/program_options.hpp>

using namespace FETools;

namespace po = boost::program_options;

bool isStringAllowed(const std::string &input, const std::vector<std::string> &allowedStrings)
{
    for (const std::string &allowedString : allowedStrings)
    {
        if (input == allowedString)
        {
            return true;
        }
    }
    return false;
}

int main(int argc, char *argv[])
{

    std::vector<std::string> allowedMeshes = {"0.5000", "0.2000", "0.1000", "0.0500", "0.0250", "0.0125"};

    po::options_description desc("Allowed options");
    desc.add_options()("help", "Produce this help message.\n")("serial", "Run the serial routine: this option also requires to specify a mesh to use by including the corresponding option when launching the code.\n")("mesh", po::value<std::string>(), "Specify the mesh to be used for the computation: \n the supported meshes are 0.5000, 0.2000, 0.1000, 0.0500, 0.0250, 0.0125. \n ")("parallel", po::value<int>(), "Run the parallel routine with the specified number of threads (to be included after -- parallel): this option also requires to specify the mesh and parallel version to be used with the corresponding options.\n")("version", po::value<unsigned short>()->default_value(0), "Specify the version of the parallel routine: \n *   0 -> Unordered Maps \n *   1 -> OMP Locks \n *   2 -> Coloring Algorithm \n")("performance_comparison", po::value<int>(), "Run the performance comparison between the serial and parallel routine with the specified number of threads: this options requires to specify the mesh and parallel version to use including the corresponging options when launching the code.\n")("convergence", "Run the convergence test in the serial setting.\n");

    try
    {
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("serial"))
        {
            // See if the mesh option is included
            if (!vm.count("mesh"))
            {
                std::cout << "The mesh option is required when running the serial routine. Try again." << std::endl;
                return 2;
            }
            std::string meshFile = vm["mesh"].as<std::string>();
            // If the mesh selected is different from the supported ones, exit and print error message
            if (!isStringAllowed(meshFile, allowedMeshes))
            {
                std::cout << "The mesh selected is not supported. Try again." << std::endl;
                return 2;
            }
            // Execute serial routine
            serialSolver solver(r);
            solver.setup(meshFile);
            solver.assemble();
            solver.solve();
            solver.process(meshFile, true);
        }
        else if (vm.count("parallel"))
        {
            THREADS = vm["parallel"].as<int>();
            unsigned short int version = vm["version"].as<unsigned short>();

            if (!vm.count("mesh"))
            {
                std::cout << "The mesh option is required when running the parallel routine. Try again." << std::endl;
                return 2;
            }
            // See if the version option is included
            if (!vm.count("version"))
            {
                std::cout << "The version option is required when running the parallel routine. Try again." << std::endl;
                return 3;
            }

            std::string meshFile = vm["mesh"].as<std::string>();
            if (!isStringAllowed(meshFile, allowedMeshes))
            {
                std::cout << "The mesh selected is not supported. Try again." << std::endl;
                return 2;
            }

            if (version == 0)
            {
                // Start psolver routine
                parallelSolver solver(r);
                solver.setup(meshFile);
                solver.assemble();
                solver.solve();
                solver.process(meshFile, true);
            }
            else if (version == 1)
            {
                // Start psolver2 routine
                parallelSolver2 solver(r);
                solver.setup(meshFile);
                solver.assemble();
                solver.solve();
                solver.process(meshFile, true);
            }
            else if (version == 2)
            {
                // Start psolverColoring routine
                parallelSolverColoring solver(r);
                solver.setup(meshFile);
                solver.assemble();
                solver.solve();
                solver.process(meshFile, true);
            }
            else
            {
                std::cout << "The version selected is not supported. Try again." << std::endl;
                return 3;
            }
        }
        else if (vm.count("convergence"))
        {
            // Execute convergence test in serial setting
            serialSolver solver(r);
            solver.convergence();
        }
        else if (vm.count("performance_comparison"))
        {
            PERFORMANCE_COMPARISON = true;
            THREADS = vm["performance_comparison"].as<int>();

            if (!vm.count("mesh"))
            {
                std::cout << "The mesh option is required when running the performance comparison. Try again." << std::endl;
                return 2;
            }
            if (!vm.count("version"))
            {
                std::cout << "The version option is required when running the performance comparison. Try again." << std::endl;
                return 3;
            }

            std::string meshFile = vm["mesh"].as<std::string>();
            if (!isStringAllowed(meshFile, allowedMeshes))
            {
                std::cout << "The mesh selected is not supported. Try again." << std::endl;
                return 2;
            }

            // Define the serial solver and time the assembly and solve routines
            serialSolver solver(r);
            solver.setup(meshFile);
            auto start = std::chrono::high_resolution_clock::now();
            solver.assemble();
            solver.solve();
            auto end = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cout << "Serial assembly and solve time: " << elapsed.count() << " ms" << std::endl;
            solver.process(meshFile, true);

            unsigned short int version = vm["version"].as<unsigned short>();

            if (version == 0)
            {
                // Define the chosen parallel solver and time the assembly and solve routines
                parallelSolver solver(r);
                solver.setup(meshFile);
                start = std::chrono::high_resolution_clock::now();
                solver.assemble();
                solver.solve();
                end = std::chrono::high_resolution_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
                std::cout << "Parallel assembly and solve time with VERSION == 0: " << elapsed.count() << " ms" << std::endl;
                solver.process(meshFile, true);
            }
            else if (version == 1)
            {
                parallelSolver2 solver(r);
                solver.setup(meshFile);
                start = std::chrono::high_resolution_clock::now();
                solver.assemble();
                solver.solve();
                end = std::chrono::high_resolution_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
                std::cout << "Parallel assembly and solve time with VERSION == 1: " << elapsed.count() << " ms" << std::endl;
                solver.process(meshFile, true);
            }
            else if (version == 2)
            {
                parallelSolverColoring solver(r);
                solver.setup(meshFile);
                start = std::chrono::high_resolution_clock::now();
                solver.assemble();
                solver.solve();
                end = std::chrono::high_resolution_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
                std::cout << "Parallel assembly and solve time with VERSION == 2: " << elapsed.count() << " ms" << std::endl;
                solver.process(meshFile, true);
            }
            else
            {
                std::cout << "The version selected is not supported. Try again." << std::endl;
                return 3;
            }

            // Execute performance comparison between serial and parallel routine
        }
        else if (vm.count("help"))
        {
            std::cout << desc << std::endl;
        }
        else
        {
            std::cout << "No valid options specified. Try again." << std::endl;
            return 1;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << ". Try again" << std::endl;
        return 1;
    }

    return 0;
}