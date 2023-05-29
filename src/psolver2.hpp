//===========================================================
//         HEADER FILE FOR THE CLASS PARALLEL SOLVER 
//===========================================================


// OSS: the definitions will be mainly the same, with some slight modifications for support parallel execution with OpenMP

#ifndef PSOL2
#define PSOL2

//#include <unordered_map>
#include "solver.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


#define NBLOCKS 4

//================================

using namespace Eigen;
using namespace FETools;

// ============================================================================================================= 

// The parallel solver class will actively be inheriting form the serial counterpart.
// The general structure will be the same, yet some methods will be organised "in a parallel setting"

class parallelSolver2
{
private:
  static constexpr unsigned short DIM = _DIM;
  static constexpr unsigned short THREADS = 4;
  //first of all, the mesh that will be saved as a member of solver,
  Mesh<DIM> _mesh;
  // a fe class solver object, to deal with the single elements of the mesh member
  // PLEASE NOTE: for the parallel implementation of the solver class, the _fe element, which is used mainly
  // in the assembly method, is substituted by private copies of this object inside the pragma omp prallel directive structured block
  // (to be found inside the assembly method)
  //SpectralFE<DIM> _fe;
  // a DoF Handler object, to generate the global mesh starting form r
  // PLEASE NOTE: the dof handler is used while iteration on the mesh elements by different threads.
  // Still none is required to modify it (on of its members), differently for the _fe object.
  // As such, we define a common _dof object to be used (queried) by the different threads
  DoFHandler<DIM> _dof;
  //some data functions
  const DiffusionCoefficient<DIM> _mu;
  const ReactionCoefficient<DIM> _sigma;
  const ForcingTerm<DIM> _f;
  const ExactSolution<DIM> _e;
  const functionZero<DIM> _g;

  //then some possible basic global matrixes to help us compute the cosen problem
  //defined with Eigen support:

  // a matrix to store the system matrix
  // In this case, we use the Eigen Matrix Type appositly created for multithreading processing environment
  SparseMatrix<double> _system_mat;
  // a vector to store the right hand side of the system
  VectorXd _rhs;
  // a vector to store the solution
  VectorXd _sol;
      
public:
  //default contructors
  parallelSolver2() = default;

  parallelSolver2(const unsigned short &dg)
      :_mesh(),
      _dof(dg,dg),
      _mu(),
      _sigma(),
      _f(),
      _e(),
      _g(),
      _system_mat(),
      _rhs(),
      _sol()
      {}
  
  void setup(const std::string&, const bool& option = true) ;
  void assemble();
  void solve(const bool& print = false, std::ostream& out= std::cout);
  void process(const std::string&, const bool& mesh_option= false);

  // standard getters
  const SparseMatrix<double>& getMat() const;
  const VectorXd& getRHS() const;
  const VectorXd& getSol() const;

  ~parallelSolver2() = default;

protected:
  void _apply_boundary(std::array<omp_lock_t,NBLOCKS>& rhs_locks);
  void _computeStiff(FETools::SpectralFE<DIM>&, std::array<omp_lock_t,NBLOCKS>& mat_locks);
  void _computeMass(FETools::SpectralFE<DIM>&, std::array<omp_lock_t,NBLOCKS>& mat_locks);
  void _computeRHS(FETools::SpectralFE<DIM>&, std::array<omp_lock_t,NBLOCKS>& rhs_locks);
  
  // a method to export the solution and mesh on a VTK file
    void _export(const std::string&, const bool& mesh_option = false) const;

    
};

#endif