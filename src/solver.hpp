 #ifndef SOL
 #define SOL

#include <iostream>
#include <map>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include "mat_utilities.hpp"
#include "elements.hpp"
//#include "../gnuplot-iostream/gnuplot-iostream.hpp"

using namespace Eigen;
using namespace FETools;
using namespace Functions;

//we define a class for the local matrix:
//each element of the mesh will have a local matrix whose element will placed
//inside the global matrix using a connectivity matrix (element_vertexes)

class serialSolver
{
private:
    //first of all, the mesh that will be saved as a member of solver,
    Mesh<Element_2D> _mesh;
    // a vector to store the map: local_index --> global_index, computed using the indexMapping method of the mesh element
    std::vector<std::vector<unsigned int>> _map;
    // a fe class solver object, to deal with the single elements of the mesh member
    SpectralFE<Element_2D> _fe;
    //some data functions
    const DiffusionCoefficient _mu;
    const ReactionCoefficient _sigma;
    const ForcingTerm _f;
    const ExactSolution _e;
    const functionZero _g;

    //then some possible basic global matrixes to help us compute the cosen problem
    //defined with Eigen support:

    // a matrix to store the stiffness matrix
    SparseMatrix<double> _Stiffness;
    // a matrix to store the mass matrix
    SparseMatrix<double> _Mass;
    // a vector to store the right hand side of the system
    VectorXd _rhs;
    // a vector to store the solution
    VectorXd _sol;
    //finally the eigen vector that will eventually contain the solution
public:
    //default contructor
    serialSolver() = default;

    serialSolver(const Mesh<Element_2D> &mesh, const unsigned int &dg/*the degree of the finite element space*/)
        :_mesh(mesh),
        // TO FIX, members for the number of elements along the two directions
        _map(),
        _fe(dg),
        _mu(),
        _sigma(),
        _f(),
        _e(),
        _g(),
        _Stiffness(),
        _Mass(),
        _rhs(),
        _sol()
        {}

    //definition of a handle-member/interface with the mesh member variable to
    //decide wether to read the mesh or generate a basic one within the program
    void setup(const unsigned int &option);
    //a method to assemble the sistem
    void assemble();
    // a method to solve the sistem
    void solve();
    // a method for processising the output data and visualize the output
    void process();
    // dummy function
    //void dummy(const Element_2D & rect);
    //destructor
    ~serialSolver() = default;

private:
    //now some local solvers to compute local matrixes using the Spectral FE classes
    MatrixXd _LocStiff() const;
    MatrixXd _LocMass() const;
    VectorXd _LocRHS() const;
    //some methods to compress the local matrix calculated on the current_element onto the global matrix
    void  _calculateStiff();
    void _calculateMass();
    void _calculateRHS();
    
    // a method to construct in a consistent way with the problem considered
    // to apply the Dirichelet boundary conditions
    void _apply_boundary( SparseMatrix<double> &mat);

    // a method to plot the results using gnuplot
    //void _plot();

    // a method to export the solution in matrix market format
    void _export();

    // a method to compute the error between the computed solution and the exact solution;
    // please note that to implement such method the exact solution must be known
    void _error(std::ostream& out = std::cout);
    
};

// operator overload to construct the local matrices more easily
inline double& operator += (double & mat_elem, const MatrixXd & product)
{
    mat_elem += product(0,0);
    return mat_elem;
};

 #endif