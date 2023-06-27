#include "psolver2.hpp"

void parallelSolver2::setup(const std::string &file_name, const bool &option)
{
    std::cout << "\nSetting up the mesh, initializing the DoF handler and algebraic structure...\n"
              << std::endl;
    // As a first step, initialize the mesh, the ID array and the DoFhandler object
    if (option == true)
    {
        // 1D test case
        if constexpr (DIM == 1)
        {
            _mesh.setMesh_csv(file_name, _dof.getDeg()[0] + 1);
            _dof.genPoints(_mesh);
        }

        // 2D test case
        else
        {
            _mesh.setMesh_csv(file_name, _dof.getDeg()[0] + 1, _dof.getDeg()[DIM - 1] + 1);
            _dof.genPoints(_mesh);
        }
    }

    else if (option == false)
    {

        // 1D test case
        if constexpr (DIM == 1)
        {
            _mesh.genMesh({0., 1.}, 10, _dof.getDeg()[0] + 1);
            _dof.genPoints(_mesh);
        }

        // 2D test case
        else
        {
            _mesh.genMesh({0., 1.}, {0., 1.}, {10, 10}, {_dof.getDeg()[0] + 1, _dof.getDeg()[DIM - 1] + 1});
            _dof.genPoints(_mesh);
        }
    }
    else
    {
        throw std::runtime_error("Something went wrong, the option chosen for the setup is not supported\n");
    }
    // Then initialise the algebraic structures
    unsigned int size = _dof.getMap()[_dof.dof_per_cell() - 1][_mesh.get_nElems() - 1];
    _system_mat.resize(size, size);
    _system_mat.reserve(size * std::log(size));

    // Insert serially the non zero elements in the matrix: this passage is necessary for the last two versions
    // of the parallel algorithm in order to avoid race conditions while allocating the space required by the
    // non zero coefficients of the matrix
    for (const Element<DIM> &elem : _mesh.getElems())
    {
        // Loop over the X DoFs of the element
        for (unsigned int i = 0; i < _dof.dof_per_cell(); i++)
        {
            // Loop over the Y DoFs of the element
            for (unsigned int j = 0; j < _dof.dof_per_cell(); j++)
            {
                unsigned int row = _dof.getMap()[i][elem.getId() - 1] - 1;
                unsigned int col = _dof.getMap()[j][elem.getId() - 1] - 1;
                _system_mat.coeffRef(row, col) = 0.;
            }
        }
    }

    // Initialize the RHS and solution vectors
    _rhs = VectorXd::Zero(size);
    _sol = VectorXd::Zero(size);
    std::cout << "Solver ready..." << std::endl;
    return;
};

void parallelSolver2::assemble()
{

    // We build a map that, for each boundary tag, stores the
    // corresponding boundary function.
    std::map<unsigned int, const Function<DIM> *> boundary_func;
    boundary_func[0] = &_g;
    boundary_func[1] = &_g;
    if constexpr (DIM == 2)
    {
        // These values of the map are used only in the 2D case
        boundary_func[2] = &_g;
        boundary_func[3] = &_g;
    }

    // Define the array of locks to be used for the parallel write of the global matrix
    std::array<omp_lock_t, NBLOCKS> mat_locks;
    std::array<omp_lock_t, NBLOCKS> rhs_locks;

    // Initialize the locks
    for (omp_lock_t &blocklock : mat_locks)
    {
        omp_init_lock(&blocklock);
    }

    for (omp_lock_t &blocklock : rhs_locks)
    {
        omp_init_lock(&blocklock);
    }

    std::cout << "=================================================================" << std::endl;
    std::cout << "Starting parallel computation of system matrix and RHS vector...\n " << std::endl;
    std::cout << "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ " << std::endl;

    std::shared_ptr<SpectralFE<DIM>> fe_ptr;
    // Define the blocksize, using the NBLOCKS macro
    unsigned int blocksize = (_system_mat.rows() + NBLOCKS - 1) / NBLOCKS;

    // Generate threads
#pragma omp parallel for num_threads(THREADS) private(fe_ptr)
    for (const Element<DIM> &elem : _mesh.getElems())
    {

        if (!fe_ptr)
        {
            fe_ptr = std::make_shared<SpectralFE<DIM>>(_dof.getDeg()[0]);
        }
        fe_ptr->update_current(elem);

        this->_computeStiff(fe_ptr, mat_locks, blocksize);
        this->_computeMass(fe_ptr, mat_locks, blocksize);
        this->_computeRHS(fe_ptr, rhs_locks, blocksize);
        this->_apply_boundary(elem, boundary_func, rhs_locks, blocksize);
    }

    std::cout << "System matrix and RHS vector computed." << std::endl;

    // Destroy the locks
    for (omp_lock_t &blocklock : mat_locks)
    {
        omp_destroy_lock(&blocklock);
    }

    for (omp_lock_t &blocklock : rhs_locks)
    {
        omp_destroy_lock(&blocklock);
    }

    return;
}

void parallelSolver2::solve(const bool &print, std::ostream &out)
{
    // BiCGSTAB method:
    Eigen::setNbThreads(THREADS);
    BiCGSTAB<SparseMatrix<double>> solver;
    solver.compute(_system_mat);

    _sol = solver.solve(_rhs); // solving the system
    if (print)
    {
        if (solver.info() != Success)
        {
            out << "The solver did not converge" << std::endl;
            return;
        }
        out << "\nThe BiCGSTAB has reached convergence in:" << std::endl;
        out << "#iterations: " << solver.iterations() << "\n\n"
            << std::endl;
        out << "The solution is:\n"
            << std::endl;
        out << _sol << std::endl;
    }

    std::cout << "\nSolution computed in #iterations " << solver.iterations() << ".\n"
              << std::endl;

    std::cout << "=================================================================" << std::endl;

    return;
}

void parallelSolver2::process(const std::string &file_name, const bool &mesh_option)
{
    this->_export(file_name, mesh_option);
    return;
}

const SparseMatrix<double> &parallelSolver2::getMat() const
{
    return _system_mat;
}

const VectorXd &parallelSolver2::getRHS() const
{
    return _rhs;
}

const VectorXd &parallelSolver2::getSol() const
{
    return _sol;
}

//  =========================================================
//          PROTECTED METHODS OF THE SOLVER CLASS
//  =========================================================

void parallelSolver2::_computeStiff(std::shared_ptr<FETools::SpectralFE<DIM>> fe_ptr, std::array<omp_lock_t, NBLOCKS> &mat_locks, const unsigned int &blocksize)
{

    // Compute the local stiffeness matrix

    MatrixXd LocStiff = ((fe_ptr->J_cell_invT() * fe_ptr->B_cell()).transpose() * (fe_ptr->D_cell() * fe_ptr->detJ()) * (fe_ptr->J_cell_invT() * fe_ptr->B_cell()));
    for (unsigned int i = 0; i < LocStiff.rows(); i++)
    {
        unsigned int global_index_x = _dof.getMap()[i][fe_ptr->getCurrent().getId() - 1] - 1;
        unsigned short bound = _dof.getPoints()[global_index_x].getBound();
        unsigned short l_i = global_index_x / blocksize;

        // Here evaluate the boundary condtion
        if (bound)
        {
            // Allow the write only if the corresponding block lock isn't owned by another thread: Set the lock
            omp_set_lock(&mat_locks[l_i]);
            // Update the global matrix writing the 1 on the diagonal
            _system_mat.coeffRef(global_index_x, global_index_x) = 1;
            // Unset the lock
            omp_unset_lock(&mat_locks[l_i]);
        }
        else
        {

            for (unsigned int j = 0; j < LocStiff.cols(); j++)
            {

                unsigned int global_index_y = _dof.getMap()[j][fe_ptr->getCurrent().getId() - 1] - 1;

                // Write on _system_mat only if the value is different form 0: limit the writes on the global matrix
                if (LocStiff.coeff(i, j) && _mu.value(fe_ptr->quadrature_point(j, _dof)))
                {
                    // Allow the write only if the corresponding block lock isn't owned by another thread: Set the lock
                    omp_set_lock(&mat_locks[l_i]);
                    // Update the global matrix with non-zero contributions
                    _system_mat.coeffRef(global_index_x, global_index_y) += LocStiff.coeff(i, j) * _mu.value(fe_ptr->quadrature_point(j, _dof));
                    // Unset the lock
                    omp_unset_lock(&mat_locks[l_i]);
                }
            }
        }
    }

    return;
};

void parallelSolver2::_computeMass(std::shared_ptr<FETools::SpectralFE<DIM>> fe_ptr, std::array<omp_lock_t, NBLOCKS> &mat_locks, const unsigned int &blocksize)
{

    // Loop over the quadrature points of the current element
    if constexpr (DIM == 1)
        for (unsigned int i = 0; i < fe_ptr->getNQ()[0]; i++)
        {
            // And compute the local contributuions by means of the quadrature  weights
            unsigned int global_index = _dof.getMap()[i][fe_ptr->getCurrent().getId() - 1] - 1;
            // Here compute the boundary flag
            unsigned short bound = _dof.getPoints()[_dof.getMap()[i][fe_ptr->getCurrent().getId() - 1] - 1].getBound();
            // Compute the block index for the lock to be queried
            unsigned short l_i = global_index / blocksize;

            if (bound == 0) // WRITE THE CONTRIBUTION ONLY IF THE DOF IS NOT ON THE BOUNDARY
            {
                omp_set_lock(&mat_locks[l_i]); // Allow for the update only if the lock is not already owned
                _system_mat.coeffRef(global_index, global_index) +=
                    fe_ptr->getQuad()[0].getW()[i] *
                    (1. / fe_ptr->getJ().coeff(0, 0)) *
                    _sigma.value(fe_ptr->quadrature_point(i, _dof));
                omp_unset_lock(&mat_locks[l_i]);
            }
        }
    else
        for (unsigned int j = 0; j < fe_ptr->getNQ()[DIM - 1]; j++)
        {
            for (unsigned int i = 0; i < fe_ptr->getNQ()[0]; i++)
            {
                // And compute the local contributuions by means of the quadrature  weights
                unsigned int local_ind = (j)*fe_ptr->getNQ()[0] + i;
                unsigned int global_index = _dof.getMap()[local_ind][fe_ptr->getCurrent().getId() - 1] - 1;
                // Here compute the boundary flag
                unsigned short bound = _dof.getPoints()[_dof.getMap()[local_ind][fe_ptr->getCurrent().getId() - 1] - 1].getBound();
                // Compute the block index for the lock to be queried
                unsigned short l_i = global_index / blocksize;

                if (bound == 0) // WRITE THE CONTRIBUTION ONLY IF THE DOF IS NOT ON THE BOUNDARY
                {
                    omp_set_lock(&mat_locks[l_i]); // Allow for the update only if the lock is not already owned
                    _system_mat.coeffRef(global_index, global_index) +=
                        fe_ptr->getQuad()[0].getW()[i] *
                        (1. / fe_ptr->getJ().coeff(0, 0)) *
                        fe_ptr->getQuad()[DIM - 1].getW()[j] *
                        (1. / fe_ptr->getJ().coeff(1, 1)) *
                        _sigma.value(fe_ptr->quadrature_point(local_ind, _dof));
                    omp_unset_lock(&mat_locks[l_i]);
                }
            }
        }

    return;
};

void parallelSolver2::_computeRHS(std::shared_ptr<FETools::SpectralFE<DIM>> fe_ptr, std::array<omp_lock_t, NBLOCKS> &rhs_locks, const unsigned int &blocksize)
{

    // Loop over the quadrature points of the current element

    if constexpr (DIM == 1)
        for (unsigned int i = 0; i < fe_ptr->getNQ()[0]; i++)
        {
            // And compute the local contribution by means of the quadrature weights
            unsigned int global_index = _dof.getMap()[i][fe_ptr->getCurrent().getId() - 1] - 1;
            // Compute the block index for the lock to be queried
            unsigned short l_i = global_index / blocksize;

            // Allow for the update only if the lock is not already owned
            omp_set_lock(&rhs_locks[l_i]);
            _rhs[global_index] += fe_ptr->getQuad()[0].getW()[i] *
                                  (1. / fe_ptr->getJ().coeff(0, 0)) *
                                  _f.value(fe_ptr->quadrature_point(i, _dof));
            omp_unset_lock(&rhs_locks[l_i]);
        }
    else
        for (unsigned int j = 0; j < fe_ptr->getNQ()[DIM - 1]; j++)
        {
            for (unsigned int i = 0; i < fe_ptr->getNQ()[0]; i++)
            {
                // And compute the local contribution by means of the quadrature weights
                unsigned int local_ind = (j)*fe_ptr->getNQ()[0] + i;
                unsigned int global_index = _dof.getMap()[local_ind][fe_ptr->getCurrent().getId() - 1] - 1;
                // Compute the block index for the lock to be queried
                unsigned short l_i = global_index / blocksize;

                // Allow for the update only if the lock is not already owend
                omp_set_lock(&rhs_locks[l_i]);
                // Updte the RHS vecotr
                _rhs[global_index] += fe_ptr->getQuad()[0].getW()[i] *
                                      (1. / fe_ptr->getJ().coeff(0, 0)) *
                                      fe_ptr->getQuad()[DIM - 1].getW()[j] *
                                      (1. / fe_ptr->getJ().coeff(1, 1)) *
                                      _f.value(fe_ptr->quadrature_point(local_ind, _dof));
                // Unset lock after write
                omp_unset_lock(&rhs_locks[l_i]);
            }
        }

    return;
}

void parallelSolver2::_apply_boundary(const Element<DIM> &elem, const std::map<unsigned int, const Function<DIM> *> &boundary_func, std::array<omp_lock_t, NBLOCKS> &rhs_locks, const unsigned int &blocksize)
{

    for (unsigned int p = 0; p < _dof.dof_per_cell() && elem.getBound(); p++)
    {
        int gindex = _dof.getMap()[p][elem.getId() - 1] - 1;
        unsigned short bound = _dof.getPoints()[gindex].getBound();
        // Compute the block index for the lock to be queried
        unsigned short l_i = gindex / blocksize;

        // Allow for the write only if the lock is free
        omp_set_lock(&rhs_locks[l_i]);
        // For each node check that it is on the boundary (and which boundary...) and eventually substitute the value of the right hand side
        _rhs[gindex] = (bound == 0) * _rhs[gindex] +
                       (bound == 1 || bound == 30 || bound == 29) * (boundary_func.at(0)->value(_dof.getPoints()[gindex])) +
                       (bound == 2) * (boundary_func.at(1)->value(_dof.getPoints()[gindex]));
        if constexpr (DIM == 2)
            _rhs[gindex] += (bound == 3 || bound == 28 || bound == 27) * (boundary_func.at(2)->value(_dof.getPoints()[gindex])) +
                            (bound == 4) * (boundary_func.at(3)->value(_dof.getPoints()[gindex]));

        omp_unset_lock(&rhs_locks[l_i]);
        // finally set to 0 the rows corresponging to those nodes on the boundary, apart from elements
        // that represent the solution on that specific node....

        // Differetly form the serial case, this is implicitly done when assembling the global matrix, considering whether
        // a certain dof is on the border or not
    }

    return;
};

void parallelSolver2::_export(const std::string &file_name, const bool &mesh_option) const
{

    // If the solution dielemory is not already present, create it

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

    if constexpr (DIM == 1)
    {
        // Define the points of the mesh
        for (const Point<DIM> &node : _dof.getPoints())
        {
            points->InsertNextPoint(node.getX(), 0., 0.);
        }

        // Set the points of the mesh
        mesh->SetPoints(points);

        // Now store the orignal mesh with the connectivity informations

        vtkSmartPointer<vtkLine> cell = vtkSmartPointer<vtkLine>::New();

        // Store the original mesh, taking into consideration only the points read on the input csv file
        if (mesh_option == false)
            for (const Element<DIM> &elem : _mesh.getElems())
            {
                cell->GetPointIds()->SetId(0, _dof.getMap()[0][elem.getId() - 1] - 1);                // left node
                cell->GetPointIds()->SetId(1, _dof.getMap()[_dof.getDeg()[0]][elem.getId() - 1] - 1); // right node
                cells->InsertNextCell(cell);
            }
        // Otherwise store the refined mesh, taking into consideration even the additional points due to a degree of the
        // Finite Element method used, for r > 1
        else
            for (const Element<DIM> &elem : _mesh.getElems())
            {
                for (unsigned int i = 0; i < (_dof.getDeg()[0]); i++)
                {
                    cell->GetPointIds()->SetId(0, _dof.getMap()[i][elem.getId() - 1] - 1);     // left node
                    cell->GetPointIds()->SetId(1, _dof.getMap()[i + 1][elem.getId() - 1] - 1); // right node
                    cells->InsertNextCell(cell);
                }
            }

        mesh->SetLines(cells);

        // Write the mesh to a VTK file
        writer->SetFileName(("./solution/solution" + file_name + ".vtp").c_str());
        writer->SetInputData(mesh);
        writer->Write();

        // Now store the solution
        solutionArray->SetName("Solution");
        solutionArray->SetNumberOfComponents(1);
        for (unsigned int i = 0; i < _dof.getPoints().size(); i++)
        {
            solutionArray->InsertNextValue(_sol.coeff(i));
        }
        mesh->GetPointData()->AddArray(solutionArray);

        // Now store the exact solution
        exactArray->SetName("Exact Solution");
        exactArray->SetNumberOfComponents(1);
        for (unsigned int i = 0; i < _dof.getPoints().size(); i++)
        {
            exactArray->InsertNextValue(_e.value(_dof.getPoints()[i]));
        }
        mesh->GetPointData()->AddArray(exactArray);

        // Write the solution to the VTK file
        writer->SetFileName(("./solution/solution" + file_name + ".vtp").c_str());
        writer->SetInputData(mesh);
        writer->Write();
    }
    else
    {

        // Loop over the points of the global
        for (const Point<DIM> &node : _dof.getPoints())
        {
            points->InsertNextPoint(node.getX(), node.getY(), 0.);
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

        if (mesh_option == false)
            for (const Element<DIM> &elem : _mesh.getElems())
            {
                cell->GetPointIds()->SetId(0, _dof.getMap()[0][elem.getId() - 1] - 1);                                                // bottom-left corner
                cell->GetPointIds()->SetId(1, _dof.getMap()[_dof.getDeg()[0]][elem.getId() - 1] - 1);                                 // bottom-right
                cell->GetPointIds()->SetId(2, _dof.getMap()[_dof.dof_per_cell() - 1][elem.getId() - 1] - 1);                          // top-right
                cell->GetPointIds()->SetId(3, _dof.getMap()[_dof.dof_per_cell() - _dof.getDeg()[DIM - 1] - 1][elem.getId() - 1] - 1); // top_left
                cells->InsertNextCell(cell);
            }
        // otherwise store the refined mesh, taking into consideration even the additional points due to a degree of the
        // Finite Element method used, for r > 1
        else
            for (const Element<DIM> &elem : _mesh.getElems())
            {

                for (unsigned int i = 0; i < (_dof.getDeg()[0] + 1) * _dof.getDeg()[DIM - 1]; i++)
                {
                    if (i == 0 || ((i + 1) % (_dof.getDeg()[0] + 1) != 0))
                    {
                        cell->GetPointIds()->SetId(0, _dof.getMap()[i][elem.getId() - 1] - 1);                        // bottom-left corner
                        cell->GetPointIds()->SetId(1, _dof.getMap()[i + 1][elem.getId() - 1] - 1);                    // bottom-right
                        cell->GetPointIds()->SetId(2, _dof.getMap()[i + _dof.getDeg()[0] + 2][elem.getId() - 1] - 1); // top-right
                        cell->GetPointIds()->SetId(3, _dof.getMap()[i + _dof.getDeg()[0] + 1][elem.getId() - 1] - 1); // top_left
                        cells->InsertNextCell(cell);
                    }
                }
            }

        mesh->SetPolys(cells);

        // Write the mesh to a VTK file
        writer->SetFileName(("./solution/solution" + file_name + ".vtp").c_str());
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

        // Add the exact solution values as point data to the mesh
        exactArray->SetName("Exact Solution");
        exactArray->SetNumberOfComponents(1);
        for (unsigned int i = 0; i < _dof.getPoints().size(); i++)
        {
            exactArray->InsertNextValue(_e.value(_dof.getPoints()[i]));
        }
        mesh->GetPointData()->AddArray(exactArray);

        // Write the solution to the VTK file
        writer->SetFileName(("./solution/solution" + file_name + ".vtp").c_str());
        writer->SetInputData(mesh);
        writer->Write();
    }

    return;
}
