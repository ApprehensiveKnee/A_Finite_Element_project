#include "psolverColoring.hpp"


void parallelSolverColoring::setup(const std::string& file_name, const bool &option)
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

    // Generate the color group vector with the greedyColoring function
    _colorGroups = greedyColoring(_mesh, _dof, MAX_COLORS);
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


void parallelSolverColoring::assemble()
{

    std::map<unsigned int, const Function<DIM> *> boundary_func;
    boundary_func[0]= &_g;
    boundary_func[1]= &_g;
    if constexpr(DIM == 2)
    {
        // These values of the map are used only in the 2D case
        boundary_func[2]= &_g;
        boundary_func[3]= &_g;
    }

    // Outer loop over the colors (cannot be parallelized)
    for(unsigned int color = 0; color < _colorGroups.size(); color++)
    {
        std::cout << "\nAssembling the system matrix for the color " << color << "...\n" << std::endl;
        std::cout << "Elements in the color " << color << ": " << std::endl;

        for(const unsigned int& elem_id : _colorGroups[color])
        {
            std::cout << elem_id << " ";
        }
        std::cout << std::endl;
        // Inner loop over the elements of the mesh (can be parallelized)
        #pragma omp parallel num_threads(THREADS) shared(boundary_func, _system_mat, _rhs, _colorGroups, _mesh, _dof, color)
        {

            // Define the _fe object for each thread
            SpectralFE<DIM> _fe(_dof.getDeg()[0]);
            
            #pragma omp for
            for(const unsigned int& elem_id : _colorGroups[color])
            {
                // Update the current element considered by the spectral _fe local solver
                _fe.update_current(_mesh.getElems()[elem_id-1]);
                
                
                #pragma omp critical
                {

                    std::cout << "Updating the element " << elem_id << " by process "<< omp_get_thread_num() << std::endl;
                }
                    // Compute the local contribution to the global Stiffeness matrix
                    this->_localStiff(_fe);

                    std::cout << "Stiffness " << elem_id << " computed.. by process"<< omp_get_thread_num() << std::endl;
                    // Compute the local contribution to the global Mass matrix
                    this->_localMass(_fe);

                    std::cout << "Mass " << elem_id << " computed.. by process"<< omp_get_thread_num() << std::endl;
                    // Compute the local contribution to the global RHS vector
                    this->_localRHS(_fe);

                    std::cout << "RHS " << elem_id << " computed.. by process"<< omp_get_thread_num() << std::endl;
                    // Apply the boundary Dirichelet condition
                    this->_apply_boundary(_mesh.getElems()[elem_id-1], boundary_func);

                    std::cout << "Element " << elem_id << " computed.. by process"<< omp_get_thread_num() << std::endl;

                
                

            }

        }

        
        
    }

    return;
}

void parallelSolverColoring::solve(const bool& print, std::ostream& out)
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

void parallelSolverColoring::process(const std::string & file_name, const bool& mesh_option)
{
    this->_export(file_name,mesh_option);
    return;
}

const SparseMatrix<double>& parallelSolverColoring::getMat() const
{
    return _system_mat;
}

const VectorXd& parallelSolverColoring::getRHS() const
{
    return _rhs;
}

const VectorXd& parallelSolverColoring::getSol() const
{
    return _sol;
}




//  =========================================================
//          PROTECTED METHODS OF THE SOLVER CLASS
//  =========================================================


void  parallelSolverColoring::_localStiff(FETools::SpectralFE<DIM>& _fe)
{
    
    {
        MatrixXd LocStiff = ((_fe.J_cell_invT() * _fe.B_cell()).transpose()*(_fe.D_cell()*_fe.detJ())*(_fe.J_cell_invT() * _fe.B_cell()));

        for(unsigned int i = 0; i < LocStiff.rows(); i++)
        {
            // Here compute the boundary flag
            unsigned short bound = _dof.getPoints()[_dof.getMap()[i][_fe.getCurrent().getId()-1]-1].getBound();

            // Here evaluate the boundary condtion
            if(bound)
            {
                // Update the unordered map just by inserting the diagonal element equal to 1
                _system_mat.coeffRef(_dof.getMap()[i][_fe.getCurrent().getId()-1]-1, _dof.getMap()[i][_fe.getCurrent().getId()-1]-1)= 1;
                
            }
            else
            {
                
                for(unsigned int j = 0; j <LocStiff.cols(); j++)
                {
                    std::cout << "Writing: " << _dof.getMap()[i][_fe.getCurrent().getId()-1]-1 << ", " << _dof.getMap()[j][_fe.getCurrent().getId()-1]-1 << std::endl;

                    //_system_mat.coeffRef(_dof.getMap()[i][_fe.getCurrent().getId()-1]-1, _dof.getMap()[j][_fe.getCurrent().getId()-1]-1) += 1; //LocStiff(i,j)*_mu.value(_fe.quadrature_point(j,_dof));  
                        
                }   
                
            }

        }
    }


    return;
    
    
};

void  parallelSolverColoring::_localMass(FETools::SpectralFE<DIM>& _fe)
{
    // As for the stiffness code, the local mass matrix is computed and then compressed onto the global matrix
    // using the ID array (_dof.getMap()) in the same way as seen in the serial solver
    if constexpr(DIM == 1)
        for(unsigned int i =0 ; i< _fe.getNQ()[0]; i++)
        {
            // And compute the local contributuions by means of the quadrature  weights
            unsigned int global_index = _dof.getMap()[i][_fe.getCurrent().getId()-1]-1;
            // Here compute the boundary flag
            unsigned short bound = _dof.getPoints()[_dof.getMap()[i][_fe.getCurrent().getId()-1]-1].getBound();
            _system_mat.coeffRef(global_index, global_index) += (bound == 0) * // WRITE THE CONTRIBUTION ONLY IF THE DOF IS NOT ON THE BOUNDARY
                                                                _fe.getQuad()[0].getW()[i] *
                                                                (1. / _fe.getJ().coeff(0, 0)) *
                                                                _sigma.value(_fe.quadrature_point(i, _dof));
        }
    else
        for(unsigned int j = 0; j < _fe.getNQ()[DIM-1]; j++)
        {
            for (unsigned int i = 0; i < _fe.getNQ()[0]; i++)
            {
                // And compute the local contributuions by means of the quadrature  weights
                unsigned int local_ind = (j)*_fe.getNQ()[0] + i;
                unsigned int global_index = _dof.getMap()[local_ind][_fe.getCurrent().getId() - 1] - 1;
                // Here compute the boundary flag

                unsigned short bound = _dof.getPoints()[_dof.getMap()[local_ind][_fe.getCurrent().getId() - 1] - 1].getBound();


                _system_mat.coeffRef(global_index, global_index) += (bound == 0) * // WRITE THE CONTRIBUTION ONLY IF THE DOF IS NOT ON THE BOUNDARY
                                                                    _fe.getQuad()[0].getW()[i] *
                                                                    (1. / _fe.getJ().coeff(0, 0)) *
                                                                    _fe.getQuad()[DIM - 1].getW()[j] *
                                                                    (1. / _fe.getJ().coeff(1, 1)) *
                                                                    _sigma.value(_fe.quadrature_point(local_ind, _dof));

            }
            
            
        }

    return;
};


void parallelSolverColoring::_localRHS(FETools::SpectralFE<DIM>& _fe)
{
    // Compute the local RHS vector and compress it onto the global RHS vector using the ID array (_dof.getMap())
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

    return;
    
}


void parallelSolverColoring::_apply_boundary(const Element<DIM>& elem,  const std::map<unsigned int, const Function<DIM> *>& boundary_func)
{  
    //inside this method we apply the direchelet b.c., known the value of the functions to be applied at the border
    
    // get the dof for each element for the DoFHandler
    // and for each dof, see if it is on the boundary (and which boundary,
    // given the element considered is on the boundary too)
    for(unsigned int p = 0; p < _dof.dof_per_cell() && elem.getBound(); p++)
    {
        int gindex = _dof.getMap()[p][elem.getId()-1]-1;
        unsigned short bound = _dof.getPoints()[gindex].getBound();

        _rhs(gindex) = (bound == 0) * _rhs(gindex) +
                              (bound == 1 || bound == 30 || bound == 29) * (boundary_func.at(0)->value(_dof.getPoints()[gindex])) +
                              (bound == 2) * (boundary_func.at(1)->value(_dof.getPoints()[gindex]));
        if constexpr (DIM == 2)
        _rhs(gindex) += (bound == 3 || bound == 28 || bound == 27) * (boundary_func.at(2)->value(_dof.getPoints()[gindex])) +
                                (bound == 4) * (boundary_func.at(3)->value(_dof.getPoints()[gindex]));
        
    }
     
    
    return;
};


void parallelSolverColoring::_export(const std::string& file_name, const bool& mesh_option) const
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







