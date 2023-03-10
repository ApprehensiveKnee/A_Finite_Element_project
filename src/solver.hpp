 #ifndef SOL
 #define SOL

#include <iostream>
#include <map>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include "mesh.hpp"
#include "mat_utilities.hpp"
#include "spectralFE.hpp"
#include "elements.hpp"
// VTK header to use in order to export the mesh ans save the solution for later visualization
//================================

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkQuad.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkDoubleArray.h>

//================================

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
    // a fe class solver object, to deal with the single elements of the mesh member
    SpectralFE<Element_2D> _fe;
    // a DoF Handler object, to generate the global mesh starting form r
    DoFHandler _dof;
    //some data functions
    const DiffusionCoefficient _mu;
    const ReactionCoefficient _sigma;
    const ForcingTerm _f;
    const ExactSolution _e;
    const functionZero _g;

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
    //default contructor
    serialSolver() = default;

    serialSolver(const unsigned short &dg)//<---the degree of the finite element space
        :_mesh(),
        _fe(dg),
        _dof(dg),
        _mu(),
        _sigma(),
        _f(),
        _e(),
        _g(),
        _system_mat(),
        _rhs(),
        _sol()
        {
            _fe.set();
        }

    //definition of a handle-member/interface with the mesh member variable to
    //decide wether to read the mesh or generate a basic one within the program
    void setup(const unsigned short &option, const std::string& = 0);
    //a method to assemble the sistem
    void assemble();
    // a method to solve the sistem
    void solve(const bool& print = false, std::ostream& out= std::cout);
    // a method to process the output data and visualize the output
    void process(const std::string&);
    // a method to define the convergence test
    void  convergence();

    // standard getters
    const SparseMatrix<double>& getMat() const;
    const VectorXd& getRHS() const;
    const VectorXd& getSol() const;

    //destructor
    ~serialSolver() = default;

//private:
    // the methods compute the local matrix for the current analyzed element
    // and compress it onto the system matrix
    // (plus a method to apply Dirichelet boundary conditions)
    void _LocStiff();
    void _LocMass();
    void _LocRHS();
    void _apply_boundary();
    //some methods to compute the global system
    void _computeStiff();
    void _computeMass();
    void _computeRHS();
    

    // // a method to export the solution and mesh on a VTK file
    void _export( const std::string&) const;

    // a method to compute the error between the computed solution and the exact solution;
    // please note that to implement such method the exact solution must be known
    double _errorL2(const bool& = false, std::ostream& out = std::cout) const;
    // a utility method to evaluate the gradient of the approximate solution using the finite differences method
    double _errorH1(const bool& = false, std::ostream& out = std::cout) const;
    
    
};

// operator overload to construct the local matrices more easily
inline double& operator += (double & mat_elem, const MatrixXd & product)
{
    mat_elem += product(0,0);
    return mat_elem;
};


 #endif