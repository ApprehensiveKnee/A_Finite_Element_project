#include "solver.hpp"



void serialSolver::setup(const unsigned short &option)
{
    std::cout << "\nSetting up the mesh, initializing the DoF handler and algebraic structure...\n" <<std::endl;
    // As a first step, initialize the mesh, the ID array and the DoFhandler object
    if(option == 1)
    {
        _mesh.setMesh_csv();
        _dof.genPoints(_mesh);
    }
        
        
    else if(option == 2)
    {
        _mesh.genMesh(0,1,0,1,50,50);
        _dof.genPoints(_mesh);
    }
    else
    {
        throw std::runtime_error("Something went wrong, the option chosen for the setup is not supported\n");
    }
    // Then initialise the algebraic structures
    unsigned int size = _dof.getMap()[_dof.dof_per_cell()-1][_mesh.get_nElems()-1];
    _system_mat.resize(size, size);
    _rhs = VectorXd::Zero(size);
    _sol = VectorXd::Zero(size);
    std::cout << "Solver ready..."<<std::endl;
    return;
};


void serialSolver::assemble()
{
    this->_computeStiff();
    this->_computeMass();
    this->_computeRHS();
    this->_apply_boundary();
    return;
}

void serialSolver::solve(const bool& print, std::ostream& out)
{
    // to solve the system we use the sovers provided by Eigen;
    // as reference here we use the BiCGSTAB method:
    
    BiCGSTAB<SparseMatrix<double>> solver;
    solver.compute(_system_mat);

    _sol = solver.solve(_rhs);  //solving the system
    if(print)
    {
        out << "\nThe BiCGSTAB has reached convergence in:" <<std::endl;
        out << "#iterations: " << solver.iterations() << "\n\n"<< std::endl;
        out << "The solution is:\n"  << std::endl;                                 
        out << _sol << std::endl;   
    }
    
    std::cout << "\nSolution computed in #iterations " << solver.iterations() <<".\n" << std::endl;
    
    // double relative_error = (x - xe).norm() / xe.norm(); // norm() is L2 norm
    // std::cout << "The relative error is:\n" << relative_error << std::endl;                       
    return;
}

void serialSolver::process()
{
    this->_export();
    this->_error();
    return;
}


const SparseMatrix<double>& serialSolver::getMat() const
{
    return _system_mat;
}

const VectorXd& serialSolver::getRHS() const
{
    return _rhs;
}

const VectorXd& serialSolver::getSol() const
{
    return _sol;
}



//  =========================================================
//          PRIVATE METHODS OF THE SOLVER CLASS
//  =========================================================

void serialSolver::_LocStiff()
{

    // evaluate the diffusion coefficient over the dofs == quadrature_points
    
    SparseMatrix<double> MU(_dof.dof_per_cell(),_dof.dof_per_cell());
    for (unsigned int i = 0; i < _dof.dof_per_cell(); i++)
    {
        MU.coeffRef(i,i)=_mu.value(_fe.quadrature_point(i,_dof));
    }

    // Compute the local matrix and compress it onto the global matrix using the ID array (_dof.getMap())
    MatrixXd LocStiff = ((_fe.J_cell_invT() * _fe.B_cell()).transpose()*(_fe.D_cell()*_fe.detJ())*(_fe.J_cell_invT() * _fe.B_cell())) * MU;
    for(unsigned int i = 0; i < LocStiff.rows(); i++)
    {
        for(unsigned int j = 0; j <LocStiff.cols(); j++)
        {
            _system_mat.coeffRef(_dof.getMap()[i][_fe.getCurrent().getId()-1]-1, _dof.getMap()[j][_fe.getCurrent().getId()-1]-1) += LocStiff.coeff(i,j);
        }
        
    }
    
};

void serialSolver::_LocMass()
{
    // evaluate the reaction coefficient over the dofs == quadrature_points
    
    SparseMatrix<double> SIGMA(_dof.dof_per_cell(),_dof.dof_per_cell());;
    for (unsigned int i = 0; i < _dof.dof_per_cell(); i++)
    {
        SIGMA.coeffRef(i,i) = _sigma.value(_fe.quadrature_point(i,_dof));
    }
    // Generate the matrix evaluating the basis function on the reference cell
    SparseMatrix<double> Bval(DIM*_dof.dof_per_cell(), _dof.dof_per_cell());
    for(unsigned int i = 0; i <DIM*_dof.dof_per_cell(); i++)
    {
        Bval.coeffRef(i, i/2) = 1.;
    }


    
    // Compute the local matrix and compress it onto the global matrix using the ID array (_dof.getMap())
    MatrixXd LocMass = (Bval.transpose()*(_fe.D_cell()*_fe.detJ())*Bval) * SIGMA;
    for(unsigned int i = 0; i < LocMass.rows(); i++)
    {
        for(unsigned int j = 0; j <LocMass.cols(); j++)
        {
            _system_mat.coeffRef(_dof.getMap()[i][_fe.getCurrent().getId()-1]-1, _dof.getMap()[j][_fe.getCurrent().getId()-1]-1) += LocMass.coeff(i,j);
        }
        
    }

    
};


void serialSolver::_LocRHS()
{
    // evaluate the forcing term over the dofs == quadrature_points
    
    VectorXd F(_dof.dof_per_cell());;
    for (unsigned int i = 0; i < _dof.dof_per_cell(); i++)
    {
        // Mitigate round-off error...
        F[i]=std::abs(_f.value(_fe.quadrature_point(i,_dof)))<tol?0.:_f.value(_fe.quadrature_point(i,_dof));
    }
    // Compute the local rhs and compress it onto the global rhs using the ID array (_dof.getMap())
    VectorXd LocRHS = ((_fe.getD()*_fe.detJ())*F);
    for(unsigned int i = 0; i < LocRHS.rows(); i++)
    {
        
            _rhs(_dof.getMap()[i][_fe.getCurrent().getId()-1]-1) += LocRHS.coeff(i) ;
        
    }
};

void  serialSolver::_computeStiff()
{

    //loop over the elements of the mesh
    for (Element_2D rect : _mesh.getElems())
    {
        //update the current element considered by the spectral _fe local solver
        _fe.update_current(rect);
        //and compute the local stiffeness matrix
        this->_LocStiff();

    }

    std::cout << "====================" << std::endl;
    std::cout << "Stiffness : COMPUTED" << std::endl;
    std::cout << "====================" << std::endl;

    return;
    
    
};

void  serialSolver::_computeMass()
{

    //loop over the elements of the mesh
    for (Element_2D rect : _mesh.getElems())
    {
        //update the current element considered by the spectral _fe local solver
        _fe.update_current(rect);
        //and compute the local mass matrix
        this->_LocMass();
        

    }
    std::cout << "====================" << std::endl;
    std::cout << "  Mass : COMPUTED" << std::endl;
    std::cout << "====================" << std::endl;
    return;
};


void serialSolver::_computeRHS()
{

    //loop over all the elements of the mesh
    for (Element_2D rect : _mesh.getElems())
    {
        //update the current element considered by the spectral _fe local solver
        _fe.update_current(rect);
        //and compute the local RHS vector
        this->_LocRHS();

    }
    std::cout << "====================" << std::endl;
    std::cout << "   rhs : COMPUTED" << std::endl;
    std::cout << "====================" << std::endl;
    return;
    
}


void serialSolver::_apply_boundary()
{

    // We build a map that, for each boundary tag, stores the
    // corresponding boundary function.
    std::map<unsigned int, const Function *> boundary_func;
    boundary_func[0]= &_g;
    boundary_func[1]= &_g;
    boundary_func[2]= &_g;
    boundary_func[3]= &_g;

    
    //inside this method we apply the direchelet b.c., known the value of the functions to be applied at the border
    //loop over all the nodes of the mesh given
    for(Element_2D rect : _mesh.getElems())
    {
        // get the dof for each element for the DoFHandler
        // and for each dof, see if it is on the boundary (and which boundary)
        for(unsigned int p = 0; p < _dof.dof_per_cell(); p++)
        {
            int gindex = _dof.getMap()[p][rect.getId()-1]-1;
            unsigned short bound = _dof.getPoints()[gindex].getBound();

            //for each node check that it is on the boundary (and which boundary...) and eventually substitute the value of the right hand side
            _rhs[gindex] = (bound==0)*_rhs[gindex] + 
                                (bound==1 || bound==30 ||bound == 29)*(boundary_func[0]->value(_dof.getPoints()[gindex])) + 
                                (bound==2 )*(boundary_func[1]->value(_dof.getPoints()[gindex])) + 
                                (bound==3 || bound== 28 ||bound == 27)*(boundary_func[2]->value(_dof.getPoints()[gindex])) +
                                (bound==4)*(boundary_func[3]->value(_dof.getPoints()[gindex]));
            //finally set to 0 the rows corresponging to those nodes on the boundary, apart from elements
            //that represent the solution on that specific node....
            if(bound)
            {
                _system_mat.prune([&gindex](int i, int j, float) {return i!=gindex;});
                _system_mat.coeffRef(gindex,gindex) = 1;
            }
            

        }

    } 
    std::cout << "====================" << std::endl;
    std::cout << " boundary : APPLIED" << std::endl;
    std::cout << "====================" << std::endl;
    return;
};



void serialSolver::_error(std::ostream& out)
{
    VectorXd error(_sol.rows());
    for(Element_2D elem : _mesh.getElems())
    {
        _fe.update_current(elem);
        for (unsigned int node = 0; node < _dof.dof_per_cell(); node++)
        {
            //compute the difference as the norm of the difference vector obtained between the solution vector numerically obtained with the FEM method
            //and the vector obtained valuating the exact solution on the nodes of the mesh (on which we computed the numerical solution)
            error[_dof.getMap()[node][elem.getId()-1]-1] = _e.value(_fe.quadrature_point(node, _dof)) - _sol[_dof.getMap()[node][elem.getId()-1]-1]; 
        }
        
        
    }

   out << "Norm of the error computed: \n\n" << error.norm() << std::endl;
   return;
}


void serialSolver::_export()
{
    // Define the points of the mesh
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    //Loop over the points of the global
    for(Point node: _dof.getPoints())
    {
        points->InsertNextPoint(node.getX(),node.getY(),0.);
    }
    

    // Create the mesh (connectivity of the points)
    vtkSmartPointer<vtkPolyData> mesh = vtkSmartPointer<vtkPolyData>::New();
    mesh->SetPoints(points);
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    

    // Now store the orignal mesh with the connectivity information
    for(Element_2D elem : _mesh.getElems())
    {

        vtkSmartPointer<vtkQuad> cell = vtkSmartPointer<vtkQuad>::New();
        cell->GetPointIds()->SetId(0,_dof.getMap()[0][elem.getId()-1]-1); // bottom-left corner
        cell->GetPointIds()->SetId(1,_dof.getMap()[_dof.getDeg()[0]][elem.getId()-1]-1);    // bottom-right
        cell->GetPointIds()->SetId(2,_dof.getMap()[_dof.dof_per_cell()-1][elem.getId()-1]-1); // top-right
        cell->GetPointIds()->SetId(3,_dof.getMap()[_dof.dof_per_cell()-_dof.getDeg()[DIM-1]-1][elem.getId()-1]-1); // top_left
        cells->InsertNextCell(cell);
        

    }
        
    mesh->SetPolys(cells);

    // Write the mesh to a VTK file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("./solution.vtp");
    writer->SetInputData(mesh);
    writer->Write();


    // Add the solution values as point data to the mesh
    vtkSmartPointer<vtkDoubleArray> solutionArray = 
    vtkSmartPointer<vtkDoubleArray>::New();
    solutionArray->SetName("Solution");
    solutionArray->SetNumberOfComponents(1);
    for (int i = 0; i < _dof.getMap()[_dof.dof_per_cell()-1][_mesh.get_nElems()-1]; i++)
    {
    solutionArray->InsertNextValue(_sol.coeff(i));
    }
    mesh->GetPointData()->AddArray(solutionArray);

    
    // Write the solution to the VTK file
    writer->SetFileName("./solution.vtp");
    writer->SetInputData(mesh);
    writer->Write();


    return;
}

