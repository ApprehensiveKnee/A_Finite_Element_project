 
//===========================================================
//             HEADER FILE FOR THE CLASS SOLVER 
//===========================================================


 #ifndef SOL
 #define SOL

#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include "mesh.hpp"
#include "functions.hpp"
#include "mat_utilities.hpp"
#include "quadrature.hpp"
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


class serialSolver
{
protected:
    static constexpr unsigned int DIM = 2;
    //first of all, the mesh that will be saved as a member of solver,
    Mesh<DIM> _mesh;
    // a fe class solver object, to deal with the single elements of the mesh member
    SpectralFE<DIM> _fe;
    // a DoF Handler object, to generate the global mesh starting form r
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
    //default contructor
    serialSolver() = default;

    serialSolver(const unsigned short &dg)//<---the degree of the finite element space
        :_mesh(),
        _fe(dg),
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

    //definition of a handle-member/interface with the mesh member variable to
    //decide wether to read the mesh or generate a basic one within the program
    virtual void setup(const unsigned short &option, const std::string& = 0);
    //a method to assemble the sistem
    virtual void assemble();
    // a method to solve the sistem
    void solve(const bool& print = false, std::ostream& out= std::cout);
    // a method to process the output data and visualize the output
    virtual void process(const std::string&);
    // a method to define the convergence test
    void  convergence();

    // standard getters
    const SparseMatrix<double>& getMat() const;
    const VectorXd& getRHS() const;
    const VectorXd& getSol() const;

    //destructor
    virtual ~serialSolver() = default;

//protected:
    // the methods compute the local matrix for the current analyzed element
    // and compress it onto the system matrix
    // (plus a method to apply Dirichelet boundary conditions)
    virtual void _LocStiff();
    virtual void _LocMass();
    virtual void _LocRHS();
    virtual void _apply_boundary();
    //some methods to compute the global system
    virtual void _computeStiff();
    virtual void _computeMass();
    virtual void _computeRHS();
    

    // // a method to export the solution and mesh on a VTK file
    void _export( const std::string&) const;

    // a method to compute the error between the computed solution and the exact solution;
    // please note that to implement such method the exact solution must be known
    double _errorL2(const bool& = false, std::ostream& out = std::cout) const;
    // a utility method to evaluate the gradient of the approximate solution using the finite differences method
    double _errorH1(const bool& = false, std::ostream& out = std::cout);
    
    
};

// operator overload to construct the local matrices more easily
inline double& operator += (double & mat_elem, const MatrixXd & product)
{
    mat_elem += product(0,0);
    return mat_elem;
};


 #endif