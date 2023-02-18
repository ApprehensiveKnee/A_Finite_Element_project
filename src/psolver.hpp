
//===========================================================
//         HEADER FILE FOR THE CLASS PARALLEL SOLVER 
//===========================================================


// OSS: the definitions will be mainly the same, with some slight modifications for support parallel execution with OpenMP

#ifndef PSOL
#define PSOL

#include "solver.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


//================================

using namespace Eigen;
using namespace FETools;

// The parallel solver class will actively be inheriting form the sertial counterpart.
// The general structure will be the same, yet some methods will be organised "in a parallel setting"

// class parallelSolver : public serialSolver
// {
// public:
//     //default contructors
//     parallelSolver() = default;

//     parallelSolver(const unsigned short &dg)
//         :serialSolver(dg)
//         {}
    
//     virtual void setup(const unsigned short &option, const std::string& = 0) override;
//     virtual void assemble() override;
//     virtual void process(const std::string&);

//     ~parallelSolver() = default;

// protected:
//     virtual void _apply_boundary() override;
//     virtual void _computeStiff() override;
//     virtual void _computeMass() override;
//     virtual void _computeRHS() override;
//     virtual void _LocMass() override;
//     virtual void _LocRHS() override;


// private:
//     // as additional member of the parallel solver class, define a 2D array of locks to 
//     // represent the access rights to a certain subblock of the global matrix when computing the local
//     // matrixes and adding their contribute to the local one. It would be much more efficient to
//     // use a lock for every element of the global matrix, but this option is just not feasible.
//     // Instead we use a "coarse grained" access rule.

//     std::map<unsigned int, std::map<unsigned int, omp_lock_t>> locks;
    
    
// };

#endif



