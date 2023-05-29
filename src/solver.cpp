#include "solver.hpp"


void serialSolver::setup(const std::string& file_name, const bool &option)
{
    std::cout << "\nSetting up the mesh, initializing the DoF handler and algebraic structure...\n" <<std::endl;
    // As a first step, initialize the mesh, the ID array and the DoFhandler object
    if(option == true)
    {   
        // 1D test case
        if constexpr(DIM == 1)
        {
            _mesh.setMesh_csv(file_name,_dof.getDeg()[0]+1);
            //_mesh.printMesh();
            _dof.genPoints(_mesh);
            //_dof.printPoints();
        }


        // 2D test case
        else
        {
            _mesh.setMesh_csv(file_name,_dof.getDeg()[0]+1, _dof.getDeg()[DIM-1]+1);
            //_mesh.printMesh();
            _dof.genPoints(_mesh);
            //_dof.printPoints();

            
        }
    }
        
        
    else if(option == false)
    {

        // 1D test case
        if constexpr(DIM == 1)
        {
            _mesh.genMesh({0., 1.},10,_dof.getDeg()[0]+1);
            //_mesh.printMesh();
            _dof.genPoints(_mesh);
            //_dof.printPoints();
        }

        // 2D test case
        else
        {
            _mesh.genMesh({0.,1.},{0.,1.},{10,10},{_dof.getDeg()[0]+1, _dof.getDeg()[DIM-1]+1});
            //_mesh.printMesh();
            _dof.genPoints(_mesh);
            //_dof.printPoints();
        }
    }
    else
    {
        throw std::runtime_error("Something went wrong, the option chosen for the setup is not supported\n");
    }
    // Then initialise the algebraic structures
    unsigned int size = _dof.getMap()[_dof.dof_per_cell()-1][_mesh.get_nElems()-1];
    _system_mat.resize(size, size);
    // Use a heuristic to estimate a reasonable value of non zero elements
    _system_mat.reserve(size);
    _rhs = VectorXd::Zero(size);
    _sol = VectorXd::Zero(size);
    std::cout << "Solver ready..."<<std::endl;
    return;
};


void serialSolver::assemble()
{
    this->_computeStiff();
    this->_computeMass();

    // ================ WRITE MATRIX TO A FILE =================
    // const static Eigen::IOFormat CSVFormat(StreamPrecision, DontAlignCols, ",", "\n");
    // Eigen::MatrixXd dense_mat = this->getMat().toDense();
    // std::ofstream file1("system_mat.csv");
    // file1 <<std::setprecision(5)<< dense_mat.format(CSVFormat);
    // ================ WRITE MATRIX TO A FILE =================

    this->_computeRHS();

    // ================ WRITE RHS TO A FILE =================
    // std::ofstream file2("system_rhs.csv");
    // file2 <<std::setprecision(5)<< this->getRHS().format(CSVFormat);
    // ================ WRITE RHS TO A FILE =================

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
                          
    return;
}

void serialSolver::process(const std::string & file_name, const bool& mesh_option)
{
    this->_export(file_name,mesh_option);
    this->_errorL2(true);
    this->_errorH1(true);
    return;
}

void  serialSolver::convergence()
{
    
    double errorL2 = 1.;
    double errorH1 = 1.;

    double h_0 = 0.1;

    std::ofstream convergence_file("convergence.csv");
        convergence_file << "h,eL2,eH1" << std::endl;
    std::ofstream convergence_table("convergence_table.csv");   
        convergence_table << "h rate, eL2 rate, eH1 rate" << std::endl;
    for(double h = 0.20; h >= 0.025; h= h/2)
    {
        // Loop over the meshes, call the solver and compute the error
        std::cout << "====================================" << std::endl;
        std::cout << " \n Mesh -->  h = " << std::fixed << std::setprecision(4) << h <<"\n\n"<< std::endl;

        std::ostringstream stream;
        stream << std::fixed << std::setprecision(4) << h;

        setup(stream.str(), true);
        assemble();
        solve();
        process(stream.str());
        

        // compute the convergence order
        if(h !=0.20 )
        {   
            convergence_table << h << "," <<std::log(this->_errorL2()/errorL2)/std::log(h/h_0)<< ","<< std::log(this->_errorH1()/errorH1)/std::log(h/h_0) << std::endl;
        }
        else
        {
            convergence_table << h <<std::endl;
        }
        
        errorL2 = this->_errorL2();
        errorH1 = this->_errorH1();

        h_0 = h;

        convergence_file << h << "," << errorL2 << "," << errorH1 <<std::endl;

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
//          PROTECTED METHODS OF THE SOLVER CLASS
//  =========================================================



void  serialSolver::_computeStiff()
{

    //loop over the elements of the mesh
    for (const Element<DIM>& rect : _mesh.getElems())
    {
        //update the current element considered by the spectral _fe local solver
        _fe.update_current(rect);
        //and compute the local stiffeness matrix
        
        // evaluate the diffusion coefficient over the dofs == quadrature_points

        // SparseMatrix<double> MU(_dof.dof_per_cell(),_dof.dof_per_cell());
        // for (unsigned int i = 0; i < _dof.dof_per_cell(); i++)
        // {
        //     MU.coeffRef(i,i)=_mu.value(_fe.quadrature_point(i,_dof));
        // }

        //Compute the local matrix and compress it onto the global matrix using the ID array (_dof.getMap())
        MatrixXd LocStiff = ((_fe.J_cell_invT() * _fe.B_cell()).transpose()*(_fe.D_cell()*_fe.detJ())*(_fe.J_cell_invT() * _fe.B_cell()));
        for(unsigned int i = 0; i < _dof.dof_per_cell(); i++)
        {
            for(unsigned int j = 0; j <_dof.dof_per_cell(); j++)
            {
                _system_mat.coeffRef(_dof.getMap()[i][_fe.getCurrent().getId()-1]-1, _dof.getMap()[j][_fe.getCurrent().getId()-1]-1) += LocStiff(i,j)*_mu.value(_fe.quadrature_point(j,_dof));

            }
        }

    }

    std::cout << "====================" << std::endl;
    std::cout << "Stiffness : COMPUTED" << std::endl;
    std::cout << "====================" << std::endl;

    return;
    
    
};

void  serialSolver::_computeMass()
{

    //loop over the elements of the mesh
    for (const Element<DIM>& rect : _mesh.getElems())
    {
        //update the current element considered by the spectral _fe local solver
        _fe.update_current(rect);
        //and compute the local mass matrix

        // As for the local mass, this will only bring about its contribution onto the diagonal of the system matrix.
        // We compute it in the following way...

        // Loop over the quadrature points of the current element
        if constexpr(DIM == 1)
            for(unsigned int i = 0; i < _fe.getNQ()[0]; i++)
            {
                unsigned int global_index = _dof.getMap()[i][_fe.getCurrent().getId()-1]-1;
                _system_mat.coeffRef(global_index, global_index) += _fe.getQuad()[0].getW()[i] *
                                                                    (1. / _fe.getJ().coeff(0, 0)) *
                                                                    _sigma.value(_fe.quadrature_point(i, _dof));
            }
        else
            for(unsigned int j = 0; j < _fe.getNQ()[DIM-1]; j++)
            {
                for (unsigned int i = 0; i < _fe.getNQ()[0]; i++)
                {
                    // And compute the local contributuions by means of the quadrature  weights
                    unsigned int local_ind = (j)*_fe.getNQ()[0]+i;
                    unsigned int global_index = _dof.getMap()[local_ind][_fe.getCurrent().getId()-1]-1;
                    _system_mat.coeffRef(global_index, global_index) += _fe.getQuad()[0].getW()[i] *
                                                                        (1. / _fe.getJ().coeff(0, 0)) *
                                                                        _fe.getQuad()[DIM - 1].getW()[j] *
                                                                        (1. / _fe.getJ().coeff(1, 1)) *
                                                                        _sigma.value(_fe.quadrature_point(local_ind, _dof));
                }
                
                
            }
        

    }
    std::cout << "====================" << std::endl;
    std::cout << "  Mass : COMPUTED" << std::endl;
    std::cout << "====================" << std::endl;
    return;
};


void serialSolver::_computeRHS()
{

    //loop over all the elements of the mesh
    for (const Element<DIM>& rect : _mesh.getElems())
    {
        //update the current element considered by the spectral _fe local solver
        _fe.update_current(rect);
        //and compute the local RHS vector

        // Loop over the quadrature points of the current element
        if constexpr(DIM == 1)
            for(unsigned int i = 0; i < _fe.getNQ()[0]; i++)
            {
                unsigned int global_index = _dof.getMap()[i][_fe.getCurrent().getId()-1]-1;
                _rhs(global_index) += _fe.getQuad()[0].getW()[i] *
                                    (1. / _fe.getJ().coeff(0, 0)) *
                                    _f.value(_fe.quadrature_point(i, _dof));
            }
        else
            for(unsigned int j = 0; j < _fe.getNQ()[DIM-1]; j++)
            {
                for (unsigned int i = 0; i < _fe.getNQ()[0]; i++)
                {
                    // And compute the local contribution by means of the quadrature weights
                    unsigned int local_ind = (j)*_fe.getNQ()[0] + i;
                    unsigned int global_index = _dof.getMap()[local_ind][_fe.getCurrent().getId() - 1] - 1;
                    _rhs(global_index) += _fe.getQuad()[0].getW()[i] *
                                        (1. / _fe.getJ().coeff(0, 0)) *
                                        _fe.getQuad()[DIM - 1].getW()[j] *
                                        (1. / _fe.getJ().coeff(1, 1)) *
                                        _f.value(_fe.quadrature_point(local_ind, _dof));
                }   
            }

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
    std::map<unsigned int, const Function<DIM> *> boundary_func;
    boundary_func[0]= &_g;
    boundary_func[1]= &_g;
    if constexpr(DIM == 2)
    {
        // These values of the map are used only in the 2D case
        boundary_func[2]= &_g;
        boundary_func[3]= &_g;
    }
    
    
    //inside this method we apply the direchelet b.c., known the value of the functions to be applied at the border
    //loop over all the nodes of the mesh given
    for(const Element<DIM>& rect : _mesh.getElems())
    {
        // get the dof for each element for the DoFHandler
        // and for each dof, see if it is on the boundary (and which boundary,
        // given the element considered is on the boundary too)
        for(unsigned int p = 0; p < _dof.dof_per_cell() && rect.getBound(); p++)
        {
            int gindex = _dof.getMap()[p][rect.getId()-1]-1;
            unsigned short bound = _dof.getPoints()[gindex].getBound();

            //for each node check that it is on the boundary (and which boundary...) and eventually substitute the value of the right hand side
            _rhs[gindex] = (bound == 0) * _rhs[gindex] +
                           (bound == 1 || bound == 30 || bound == 29) * (boundary_func[0]->value(_dof.getPoints()[gindex])) +
                           (bound == 2) * (boundary_func[1]->value(_dof.getPoints()[gindex]));
            if constexpr(DIM == 2)
            _rhs[gindex] += (bound == 3 || bound == 28 || bound == 27) * (boundary_func[2]->value(_dof.getPoints()[gindex])) +
                            (bound == 4) * (boundary_func[3]->value(_dof.getPoints()[gindex]));
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
    double sum(0);

    // loop over all the elements
    for(const Element<DIM>& e: _mesh.getElems())
    {
        if constexpr(DIM == 1)
            for(unsigned int i = 0; i < _fe.getNQ()[0]; i++)
            {
                // detrtmine the global index of the node currently considered
                unsigned int gindex = _dof.getMap()[i][e.getId()-1]-1;

                // compute the local weight
                double weight = (_fe.getQuad()[0].getW()[i] *
                                (1. / _fe.getJ().coeff(0, 0)));

                sum += (_e.value(_dof.getPoints()[gindex]) -
                        _sol[gindex]) *
                        (_e.value(_dof.getPoints()[gindex]) -
                        _sol[gindex]) *
                        weight;
            }
        else
        //loop over all the quadrature points of the element
            for(unsigned int j = 0; j < _fe.getNQ()[DIM-1]; j++)
            {
                for (unsigned int i = 0; i < _fe.getNQ()[0]; i++)
                {
                    // determine the index of local quadrature point currently considered
                    unsigned int q = j * _fe.getNQ()[0]+i;
                    // and the global index of the node currently considered
                    unsigned int gindex = _dof.getMap()[q][e.getId()-1]-1;

                    // compute the local weight
                    double weight = (_fe.getQuad()[0].getW()[i] *
                                    (1. / _fe.getJ().coeff(0, 0)) *
                                    _fe.getQuad()[DIM - 1].getW()[j] *
                                    (1. / _fe.getJ().coeff(1, 1)));

                    sum += (_e.value(_dof.getPoints()[gindex]) - _sol.coeff(gindex)) * (_e.value(_dof.getPoints()[gindex]) - _sol.coeff(gindex)) * weight;
                }
                
                
            }

    }

     
    if(print)
    {
        out << "Norm L2 of error computed: " << std::scientific << std::sqrt(sum) << std::endl;
    }
        
    return std::sqrt(sum);
}



double serialSolver::_errorH1(const bool& print,std::ostream& out)
{
    // first define the temp variable to store the sum of 

    double sum(0);
    

   //  Loop over each element 
   for(const Element<DIM>& elem : _mesh.getElems())
   {
        // For each element get the evaluation of the exact solution at the quadrature nodes
        // and compute the numerical derivative using the derivative matrix of the element computed by the fe object
        _fe.update_current(elem);

        if constexpr(DIM == 1)
        {
            // A simple eigen vector to store the values of the numerical solution on the local 
            // quadrature nodes, to be used to compute the derivative of the solution on the same nodes
            VectorXd un(_fe.getNQ()[0]);
            
            // Loop over the quadrature points of the element
            for(unsigned int i = 0; i < _fe.getNQ()[0]; i++)
            {
                // determine the global index of the node currently considered
                unsigned int gindex = _dof.getMap()[i][elem.getId()-1]-1;
                
                // get the exact solution on the quadrature nodes and place its values in a matrix
                un(i) = _sol[gindex];
            }
            
            // compute the derivative of the numerical solution on the quadrature nodes
            VectorXd d_un = (_fe.getJ().coeff(0, 0))*(_fe.getB()[0]* un);
            
            // loop over the quadrature points of the element
            for(unsigned int i = 0; i < _fe.getNQ()[0]; i++)
            {
                // Determine the global index of the node currently considered
                unsigned int gindex = _dof.getMap()[i][elem.getId()-1]-1;
                
                // compute the local weight
                double weight = (_fe.getQuad()[0].getW()[i] *
                                (1. / _fe.getJ().coeff(0, 0)));
                
                // compute the sum of the square of the difference between the exact and numerical derivative
                sum += (d_un(i) - _e.grad(_dof.getPoints()[gindex])[0]) *
                        (d_un(i) - _e.grad(_dof.getPoints()[gindex])[0]) *
                        weight;
            }
        }
        else
        {

            //  An eigen matrix to store the values of the numerical solution on the local quadrature nodes
            // necessary to then compute the evaluation of the derivatives
            MatrixXd un(_fe.getNQ()[0],_fe.getNQ()[DIM-1]);
            
            // Loop over the quadrature points of the element
            for(unsigned int j = 0; j < _fe.getNQ()[DIM-1]; j++)
            {
                for (unsigned int i = 0; i < _fe.getNQ()[0]; i++)
                {
                    // determine the index of local quadrature point currently considered
                    unsigned int q = (j)*_fe.getNQ()[0]+i;
                    // and the global index of the node currently considered
                    unsigned int gindex = _dof.getMap()[q][elem.getId()-1]-1;
                    
                    // get the exact solution on the quadrature nodes and place its values in a matrix
                    un(i,j) = _sol.coeff(gindex);
                    
                }
            }


            //Compute the numerical derivative
            MatrixXd ux = (_fe.getJ().coeff(0, 0))*(_fe.getB()[0]* un);
            MatrixXd uy = (_fe.getJ().coeff(1, 1))*(_fe.getB()[1]* un.transpose()).transpose();


            // Loop over the quadrature points of the element to sum all the differences over the quadrature points
            for(unsigned int j = 0; j < _fe.getNQ()[DIM-1]; j++)
            {
                for (unsigned int i = 0; i < _fe.getNQ()[0]; i++)
                {
                    // determine the index of local quadrature point currently considered
                    unsigned int q = (j)*_fe.getNQ()[0]+i;
                    
                    // sum all the contributions for the element currenlty considered
                    sum += ((un.coeff(i, j) - _e.value(_fe.quadrature_point(q, _dof))) *
                                (un.coeff(i, j) - _e.value(_fe.quadrature_point(q, _dof))) + // (U - U_n)^2
                            (ux.coeff(i, j) - _e.grad(_fe.quadrature_point(q, _dof))[0]) *
                                (ux.coeff(i, j) - _e.grad(_fe.quadrature_point(q, _dof))[0]) + // (dUx -dU_nx)^2
                            (uy.coeff(i, j) - _e.grad(_fe.quadrature_point(q, _dof))[1]) *
                                (uy.coeff(i, j) - _e.grad(_fe.quadrature_point(q, _dof))[1])) *
                        (_fe.getQuad()[0].getW()[i] *
                            (1. / _fe.getJ().coeff(0, 0)) *
                            _fe.getQuad()[DIM - 1].getW()[j] *
                            (1. / _fe.getJ().coeff(1, 1))); // (dUy -dU_ny)^2
                }
            }

        }
        
   }

   if(print)
        out << "Norm H1 of error computed: " << std::scientific <<std::sqrt(sum) << std::endl;

    return  std::sqrt(sum);
    

}



void serialSolver::_export(const std::string& file_name, const bool& mesh_option) const
{

    // If the solution directory is not already present, create it

    if (!std::filesystem::exists("./solution/"))
    {
        std::filesystem::create_directory("./solution/");
    }

    // Define the points of the mesh
    vtkSmartPointer<vtkPoints> points = 
    vtkSmartPointer<vtkPoints>::New();
    // Create the mesh (connectivity of the points)
    vtkSmartPointer<vtkPolyData> mesh = 
    vtkSmartPointer<vtkPolyData>::New();
    // Create the cells of the mesh
    vtkSmartPointer<vtkCellArray> cells = 
    vtkSmartPointer<vtkCellArray>::New();
    // Create the writer to the vtk file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    // Create the solution array
    vtkSmartPointer<vtkDoubleArray> solutionArray = 
    vtkSmartPointer<vtkDoubleArray>::New();
    // Create the exact solution array
    vtkSmartPointer<vtkDoubleArray> exactArray =
    vtkSmartPointer<vtkDoubleArray>::New();

    if constexpr(DIM == 1)
    {
        std::cout << "Nodes" << std::endl;
        //Define the points of the mesh
        for(const Point<DIM>& node: _dof.getPoints())
        {
            points->InsertNextPoint(node.getX(),0.,0.);
        }
        
        // Set the points of the mesh
        mesh->SetPoints(points);
        
        
        // Now store the orignal mesh with the connectivity informations

        vtkSmartPointer<vtkLine> cell = vtkSmartPointer<vtkLine>::New();

        // Store the original mesh, taking into consideration only the points read on the input csv file
        if(mesh_option == false)
            for(const Element<DIM>& elem : _mesh.getElems())
            {
                cell->GetPointIds()->SetId(0,_dof.getMap()[0][elem.getId()-1]-1);  // left node
                cell->GetPointIds()->SetId(1,_dof.getMap()[_dof.getDeg()[0]][elem.getId()-1]-1); // right node
                cells->InsertNextCell(cell);
            }
        // Otherwise store the refined mesh, taking into consideration even the additional points due to a degree of the
        // Finite Element method used, for r > 1
        else
            for(const Element<DIM>& elem : _mesh.getElems())
            {
                for(unsigned int i = 0; i < (_dof.getDeg()[0]); i++)
                {
                    cell->GetPointIds()->SetId(0,_dof.getMap()[i][elem.getId()-1]-1);  // left node
                    cell->GetPointIds()->SetId(1,_dof.getMap()[i+1][elem.getId()-1]-1); // right node
                    cells->InsertNextCell(cell);
                }
            }
        

        mesh->SetLines(cells);

        // Write the mesh to a VTK file
        writer->SetFileName(("./solution/solution"+file_name+".vtp").c_str());
        writer->SetInputData(mesh);
        writer->Write();

        // Now store the solution
        solutionArray->SetName("Solution");
        solutionArray->SetNumberOfComponents(1);
        for(unsigned int i = 0; i < _dof.getPoints().size(); i++)
        {
            solutionArray->InsertNextValue(_sol.coeff(i));
        }
        mesh->GetPointData()->AddArray(solutionArray);

        
        // Now store the exact solution
        exactArray->SetName("Exact Solution");
        exactArray->SetNumberOfComponents(1);
        for(unsigned int i = 0; i < _dof.getPoints().size(); i++)
        {
            exactArray->InsertNextValue(_e.value(_dof.getPoints()[i]));
        }
        mesh->GetPointData()->AddArray(exactArray);

        // Write the solution to the VTK file
        writer->SetFileName(("./solution/solution"+file_name+".vtp").c_str());
        writer->SetInputData(mesh);
        writer->Write();
    }
    else
    {
        
        //Loop over the points of the global
        for(const Point<DIM>& node: _dof.getPoints())
        {
            points->InsertNextPoint(node.getX(),node.getY(),0.);
        }
        
        // Set the points of the mesh
        mesh->SetPoints(points);
        

        // Now store the orignal mesh with the connectivity information

        // Either store the orginally analyzed mesh (the connectivity information only takes into consideration
        // the points read on the input csv file: the mesh stored on the .vtp file is the same as the one read):
        // in this case the, when visualizing the solution using the Plot Over Line (for the 1D solution) of the 
        // Warp By Scalar (for the 2D solution) on ParaView, the solution is plotted on the original mesh
        // =====> the filter does not take into consideration the solution computed on the refined mesh, while its values
        // are still visible if the Point Gaussian filter is selecter for the visualization of the data.

        vtkSmartPointer<vtkQuad> cell = vtkSmartPointer<vtkQuad>::New();

        if(mesh_option == false)
            for(const Element<DIM>& elem : _mesh.getElems())
            {
                cell->GetPointIds()->SetId(0,_dof.getMap()[0][elem.getId()-1]-1); // bottom-left corner
                cell->GetPointIds()->SetId(1,_dof.getMap()[_dof.getDeg()[0]][elem.getId()-1]-1);    // bottom-right
                cell->GetPointIds()->SetId(2,_dof.getMap()[_dof.dof_per_cell()-1][elem.getId()-1]-1); // top-right
                cell->GetPointIds()->SetId(3,_dof.getMap()[_dof.dof_per_cell()-_dof.getDeg()[DIM-1]-1][elem.getId()-1]-1); // top_left
                cells->InsertNextCell(cell);
            }
        // otherwise store the refined mesh, taking into consideration even the additional points due to a degree of the 
        // Finite Element method used, for r > 1
        else
        for(const Element<DIM>& elem : _mesh.getElems())
        {
            
            for(unsigned int i = 0; i < (_dof.getDeg()[0]+1)*_dof.getDeg()[DIM-1]; i++)
            {
                    if(i==0 || ((i + 1)%(_dof.getDeg()[0]+1)!= 0))
                    {
                        cell->GetPointIds()->SetId(0,_dof.getMap()[i][elem.getId()-1]-1); // bottom-left corner
                        cell->GetPointIds()->SetId(1,_dof.getMap()[i+1][elem.getId()-1]-1);    // bottom-right
                        cell->GetPointIds()->SetId(2,_dof.getMap()[i + _dof.getDeg()[0] + 2][elem.getId()-1]-1); // top-right
                        cell->GetPointIds()->SetId(3,_dof.getMap()[i + _dof.getDeg()[0] + 1][elem.getId()-1]-1); // top_left
                        cells->InsertNextCell(cell);

                    }
                    

            }

        }
            
        mesh->SetPolys(cells);

        // Write the mesh to a VTK file
        writer->SetFileName(("./solution/solution"+file_name+".vtp").c_str());
        writer->SetInputData(mesh);
        writer->Write();

        // Add the solution values as point data to the mesh
        solutionArray->SetName("Solution");
        solutionArray->SetNumberOfComponents(1);
        for (unsigned int i = 0; i < _dof.getPoints().size(); i++)
        {
            solutionArray->InsertNextValue(_sol.coeff(i));
        }
        mesh->GetPointData()->AddArray(solutionArray);

        //Add the exact solution values as point data to the mesh
        exactArray->SetName("Exact Solution");
        exactArray->SetNumberOfComponents(1);
        for (unsigned int i = 0; i < _dof.getPoints().size(); i++)
        {
            exactArray->InsertNextValue(_e.value(_dof.getPoints()[i]));
        }
        mesh->GetPointData()->AddArray(exactArray);

        // Write the solution to the VTK file
        writer->SetFileName(("./solution/solution"+file_name+".vtp").c_str());
        writer->SetInputData(mesh);
        writer->Write();

    }

    return;
}


