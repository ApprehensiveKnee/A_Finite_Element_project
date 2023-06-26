//===========================================================
//         HEADER FILE FOR THE PARALLEL SOLVER CLASS (3)
//===========================================================
// PARALLEL IMPLEMENTATION BASED ON THE USE OF THE OPENMP LIBRARY

#ifndef PSOLCOLOR
#define PSOLCOLOR

// Define the macro for the coloring algorithm
#define COLORING
// Include other dependencies and libraries
#include "solver.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

//================================

using namespace Eigen;
using namespace FETools;

// As for the other parallel solvers, here too the solver class is pretty much the same as the others:
// what changes is the definitions of the methods

class parallelSolverColoring
{
private:
  static constexpr unsigned short MAX_COLORS = 4; // Please note that, for a simple cases such as the ones considered, for a rectangular mesh,
  // the maximum number of possible colors has been chosen considering the maximum number of elements
  // sharing the same node, which is 4 in the worst case...
  std::vector<std::vector<unsigned int>> _colorGroups; // A vector of vectors to store the elements of the mesh grouped by color...

  // the mesh that will be saved as a member of solver...
  Mesh<DIM> _mesh;
  // a DoF Handler object, to generate the global mesh starting form r
  DoFHandler<DIM> _dof;
  // some data functions
  const DiffusionCoefficient<DIM> _mu;
  const ReactionCoefficient<DIM> _sigma;
  const ForcingTerm<DIM> _f;
  const ExactSolution<DIM> _e;
  const functionZero<DIM> _g;

  // then some possible basic global matrixes to help us compute the cosen problem
  // defined with Eigen support:

  // a matrix to store the system matrix
  SparseMatrix<double> _system_mat;
  // a vector to store the right hand side of the system
  VectorXd _rhs;
  // a vector to store the solution
  VectorXd _sol;
  // finally the eigen vector that will eventually contain the solution

public:
  // default contructors
  parallelSolverColoring() = default;

  parallelSolverColoring(const unsigned short &dg)
      : _colorGroups(MAX_COLORS),
        _mesh(),
        _dof(dg, dg),
        _mu(),
        _sigma(),
        _f(),
        _e(),
        _g(),
        _system_mat(),
        _rhs(),
        _sol()
  {
  }

  void setup(const std::string &, const bool &option = true);
  void assemble();
  void solve(const bool &print = false, std::ostream &out = std::cout);
  void process(const std::string &, const bool &mesh_option = false);

  // standard getters
  const SparseMatrix<double> &getMat() const;
  const VectorXd &getRHS() const;
  const VectorXd &getSol() const;

  ~parallelSolverColoring() = default;

private:
  void _apply_boundary(const Element<DIM> &, const std::map<unsigned int, const Function<DIM> *> &);
  void _localStiff(std::shared_ptr<FETools::SpectralFE<DIM>>);
  void _localMass(std::shared_ptr<FETools::SpectralFE<DIM>>);
  void _localRHS(std::shared_ptr<FETools::SpectralFE<DIM>>);

  // a method to export the solution and mesh on a VTK file
  void _export(const std::string &, const bool &mesh_option = false) const;
};

#endif