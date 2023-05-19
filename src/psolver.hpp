
//===========================================================
//         HEADER FILE FOR THE CLASS PARALLEL SOLVER 
//===========================================================


// OSS: the definitions will be mainly the same, with some slight modifications for support parallel execution with OpenMP

#ifndef PSOL
#define PSOL

//#include <unordered_map>
#include "solver.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


//================================

using namespace Eigen;
using namespace FETools;

// =================== DEFINITION OF TEMPLATED HASH FUNCTIONS FOR UNORDERED MAPS OF STD. PAIRS =================

// Here the idea for the hash function has been taken from the Boost libraries

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}


namespace std
{
  template<typename S, typename T> 
  struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };

}


// ============================================================================================================= 

// The parallel solver class will actively be inheriting form the serial counterpart.
// The general structure will be the same, yet some methods will be organised "in a parallel setting"

class parallelSolver
{
private:
  static constexpr unsigned short DIM = 2;
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
  SparseMatrix<double> _system_mat;
  // a vector to store the right hand side of the system
  VectorXd _rhs;
  // a vector to store the solution
  VectorXd _sol;
  //finally the eigen vector that will eventually contain the solution
      
public:
  //default contructors
  parallelSolver() = default;

  parallelSolver(const unsigned short &dg)
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
  
  void setup(const std::string&, const bool &option = true) ;
  void assemble();
  void solve(const bool& print = false, std::ostream& out= std::cout);
  void process(const std::string&);

  // standard getters
  const SparseMatrix<double>& getMat() const;
  const VectorXd& getRHS() const;
  const VectorXd& getSol() const;

  ~parallelSolver() = default;

protected:
  void _apply_boundary(const Element<DIM>&,std::unordered_map<unsigned int, double>&, const std::map<unsigned int, const Function<DIM> *>&);
  void _localStiff(FETools::SpectralFE<DIM>&, std::unordered_map<std::pair<unsigned int, unsigned int>, double>&);
  void _localMass(FETools::SpectralFE<DIM>&, std::unordered_map<std::pair<unsigned int, unsigned int>, double>&);
  void _localRHS(FETools::SpectralFE<DIM>&, std::unordered_map<unsigned int, double>&);

  // a method to export the solution and mesh on a VTK file
  void _export( const std::string&) const;

    
};

#endif



