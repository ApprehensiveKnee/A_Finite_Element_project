#include "solver.hpp"



void serialSolver::setup(const unsigned short &option, const std::string& file_name)
{
    std::cout << "\nSetting up the mesh, initializing the DoF handler and algebraic structure...\n" <<std::endl;
    // As a first step, initialize the mesh, the ID array and the DoFhandler object
    if(option == 1)
    {
        _mesh.setMesh_csv(file_name);
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

void serialSolver::process(const std::string & file_name)
{
    this->_export(file_name);
    this->_errorL2();
    this->_errorH1();
    return;
}

void  serialSolver::convergence()
{
    for(double h = 0.10; h >= 0.0125; h= h/2)
    {
        // Loop over the meshes, call the solver and compute the error
        std::cout << "====================================" << std::endl;
        std::cout << "  Mesh -->    h = " << std::fixed << std::setprecision(4) << h <<"\n\n"<< std::endl;

        std::ostringstream stream;
        stream << std::fixed << std::setprecision(4) << h;

        setup(1,stream.str());
        assemble();
        solve();
        process(stream.str());

        double errorL2 = 1.;
        double errorH1 = 1.;

        double h_0 = 0.1;
        

        // compute the convergence order
        if(h !=0.10 )
        {   
            std::cout << "===========================================" << std::endl;
            std::cout << "  h : "<< h << std::endl;
            std::cout << "  L2 error discrepancy : "<<this->_errorL2() - errorL2 << std::endl;
            std::cout << "  convergence order for error L2 : " << std::endl;
            std::cout << std::log(this->_errorL2() - errorL2)/std::log(h - h_0) << std::endl;
            std::cout << "  H1 error discrepancy : "<< this->_errorH1() - errorH1  << std::endl;
            std::cout << "  convergence order for error H2 : " << std::endl;
            std::cout << std::log(this->_errorH1() - errorH1)/std::log(h - h_0) << std::endl;

            std::cout << "===========================================" << std::endl;
        }

        errorL2 = this->_errorL2(true);
        errorH1 = this->_errorH1();

        h_0 = h;

    }
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
    for(const Element_2D& rect : _mesh.getElems())
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



double serialSolver::_errorL2(const bool& print, std::ostream& out) const
{
    VectorXd error(_sol.rows());
    unsigned int ind(0);
    double sum = 0;
    for(const Point& p : _dof.getPoints())
    {
        error[ind] = _e.value(p) - _sol[ind]; 
        ind++;
    }
        
    if(print)
    {
        out << "Norm L2 of error computed: \n\n" << error.norm() << std::endl;
    }
        
    return error.norm();
}


void serialSolver::_export(const std::string& file_name) const
{
    // Define the points of the mesh
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    //Loop over the points of the global
    for(const Point& node: _dof.getPoints())
    {
        points->InsertNextPoint(node.getX(),node.getY(),0.);
    }
    

    // Create the mesh (connectivity of the points)
    vtkSmartPointer<vtkPolyData> mesh = vtkSmartPointer<vtkPolyData>::New();
    mesh->SetPoints(points);
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    

    // Now store the orignal mesh with the connectivity information
    for(const Element_2D& elem : _mesh.getElems())
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
    writer->SetFileName(("./solution"+file_name+".vtp").c_str());
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
    writer->SetFileName(("./solution"+file_name+".vtp").c_str());
    writer->SetInputData(mesh);
    writer->Write();


    return;
}



double serialSolver::_errorH1(const bool& print,std::ostream& out) const
{
    // first define the temp variable to store the sum of 

    double sum(0);

    // A reorganization of the solution is in order as to better "move around" while computing the evaluation of the gradient
    // For this purpouse we generate a 2D array to store the values of the solution "as they would appear of the grid": essentially
    // first all the values of the solution on the first row of points (y = 0, x = 0,1,2,... ), then all the values on the second row of points (y = 1, x= 0,1,2,...) and so on

    const unsigned int nx = (_mesh.get_nEx()+1)+(_fe.getDeg()-1)*_mesh.get_nEx();
    const unsigned int ny = (_mesh.get_nEy()+1)+(_fe.getDeg()-1)*_mesh.get_nEy();
    std::vector<std::vector<double>>sol(nx, std::vector<double>(ny));
    for(unsigned int i = 0; i < ny; ++i)
    {
        for(unsigned int j = 0; j < nx; ++j)
        {
            if(j < _fe.getDeg()+1)
            {
                sol[i][j] = _sol.coeff(i*(_fe.getDeg()+1)+j);
            }
            else
            {
                sol[i][j] =_sol.coeff( i*_fe.getDeg()+ j + ((ny-1)*(_fe.getDeg()*((j-1)/_fe.getDeg())+1)));
            }
            //std::cout << sol[i][j] << " ";
        }
        //std::cout << std::endl;
        
    }

    // Allocate the resources for the components of the gradient of the approximated solution
    double dy;
    double dx;
    double hx = (_dof.getPoints()[1].getX() - _dof.getPoints()[0].getX()); //step in x direction
    double hy = (_dof.getPoints()[_dof.getDeg()[0]+1].getY() - _dof.getPoints()[0].getY()); //step in y direction

    
    for(unsigned int i = 0; i < ny; ++i)
    {
        for(unsigned int j = 0; j < nx; ++j)
        {
            // for the points inside the 2D domain...
            if(i != 0 && j != 0 && i != ny-1 && j !=nx-1)
            {
                // ...evalute the derivative with the central finite difference
                dy = (sol[i + 1][j] - sol[i - 1][j]) / (2 * hy);
                dx = (sol[i][j + 1] - sol[i][j - 1]) / (2 * hx);
            }
            else // otherwhise use some other finite difference method
            {
                if(i == 0)
                {
                    //evaluate dy by forward finite difference method means
                    dy = (sol[i + 1][j] - sol[i][j]) / hy;
                    if(j == 0)
                    {
                        //evaluate dx by forward finite difference method means
                        dx = (sol[i][j+1] - sol[i][j]) / hx;
                    }
                    else if (j == nx-1)
                    {
                        //evaluate dx by backwards finite difference method means
                        dx = (sol[i][j] - sol[i][j-1]) / hx;
                    }
                    else
                    {
                        //evaluate dx by central finite difference method means
                        dx = (sol[i][j + 1] - sol[i][j - 1]) / (2 * hx);
                    }
                }
                else if(i == ny-1)
                {
                    //evaluate dy by backwards finite difference method means
                    dy = (sol[i][j] - sol[i-1][j]) / hy;
                    if(j == 0)
                    {
                        //evaluate dx by forward finite difference method means
                        dx = (sol[i][j+1] - sol[i][j]) / hx;
                    }
                    else if (j == nx-1)
                    {
                        //evaluate dx by backwards finite difference method means
                        dx = (sol[i][j] - sol[i][j-1]) / hx;
                    }
                    else
                    {
                        //evaluate dx by central finite difference method means
                        dx = (sol[i][j + 1] - sol[i][j - 1]) / (2 * hx);
                    }

                }
                else
                {
                    //evaluate dy by central finite difference method means
                    dy = (sol[i+1][j] - sol[i-1][j]) / (2*hy);
                    if(j == 0)
                    {
                        //evaluate dx by forward finite difference method means
                        dx = (sol[i][j+1] - sol[i][j]) / hx;
                    }
                    else if (j == nx-1)
                    {
                        //evaluate dx by backwards finite difference method means
                        dx = (sol[i][j] - sol[i][j-1]) / hx;
                    }
                    else
                    {
                        //evaluate dx by central finite difference method means
                        dx = (sol[i][j + 1] - sol[i][j - 1]) / (2 * hx);
                    }
                } 
            }

            // after computing the values of the gradient over a certain point, compute the difference
            // between the gradient of the exact solution and the approximate one
            double diff;
            double diff_grad;
            if(j < _fe.getDeg()+1)
            {
                unsigned int ind = i * (_fe.getDeg() + 1)+j;
                diff = sol[i][j] - _e.value(_dof.getPoints()[ind]);
                diff = diff*diff;
                diff_grad = (dx - _e.grad(_dof.getPoints()[ind])[0]) *
                            (dx - _e.grad(_dof.getPoints()[ind])[0]) +
                            (dy - _e.grad(_dof.getPoints()[ind])[1]) *
                            (dy - _e.grad(_dof.getPoints()[ind])[1]);
            }
            else
            {
                unsigned int ind = i*_fe.getDeg()+ j + ((ny-1)*(_fe.getDeg()*((j-1)/_fe.getDeg())+1));
                diff = sol[i][j] - _e.value(_dof.getPoints()[ind]);
                diff = diff*diff;
                diff_grad = (dx - _e.grad(_dof.getPoints()[ind])[0]) *
                            (dx - _e.grad(_dof.getPoints()[ind])[0]) +
                            (dy - _e.grad(_dof.getPoints()[ind])[1]) *
                            (dy - _e.grad(_dof.getPoints()[ind])[1]);
            }
            sum +=  diff_grad + diff;


            
        }
        
    }


    if(print)
        out << "Norm H1 of error computed :" << std::sqrt(sum) << std::endl;

    return  std::sqrt(sum);
    
 

}

