#include "psolver2.hpp"


void parallelSolver2::setup(const std::string& file_name, const bool &option)
{
    std::cout << "\nSetting up the mesh, initializing the DoF handler and algebraic structure...\n" <<std::endl;
    // As a first step, initialize the mesh, the ID array and the DoFhandler object
    if(option == true)
    {
        _mesh.setMesh_csv(file_name,_dof.getDeg()[0]+1, _dof.getDeg()[DIM-1]+1);
        //_mesh.printMesh();
        _dof.genPoints(_mesh);
        //_dof.printPoints();
    }
        
        
    else if(option == false)
    {
        _mesh.genMesh({0,1},{0,1},{10,10},{_dof.getDeg()[0]+1, _dof.getDeg()[DIM-1]+1});
        //_mesh.printMesh();
        _dof.genPoints(_mesh);
        //_dof.printPoints();
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


void parallelSolver2::assemble()
{

    // Define the array of locks to be used for the parallel write of the global matrix
    std::array<omp_lock_t,NBLOCKS> mat_locks;
    std::array<omp_lock_t, NBLOCKS> rhs_locks;

    // Initialize the locks
    for(omp_lock_t& blocklock : mat_locks)
    {
        omp_init_lock(&blocklock);
    }

    for(omp_lock_t& blocklock : rhs_locks)
    {
        omp_init_lock(&blocklock);
    }

    // Generate threads
    #pragma omp parallel num_threads(4) 
    {
        #pragma omp barrier
        
        SpectralFE<DIM> _fe(_dof.getDeg()[0]);


        #pragma omp single
        {
            std::cout << "================================================================="<< std::endl;
            std::cout << " Starting parallel computation of system matrix and RHS vector...\n - "<< std::endl;
        }

        // Compute the local contribution to the global Stiffeness matrix
        this->_computeStiff(_fe, mat_locks);

        // Compute the local contribution to the global Mass matrix
        //this->_computeMass(_fe,mat_locks);

        // Compute the local contribution to the global RHS vector
        //this->_computeRHS(_fe, rhs_locks);

        // Apply the boundary Dirichelet condition
        //this->_apply_boundary(rhs_locks);

        #pragma omp barrier

        #pragma omp single
        {
            std::cout << " System matrix and RHS vector computed." << std::endl;

        }
  
        
    }
    std::cout << " HERE " << std::endl;
    // Destroy the locks
    for(omp_lock_t& blocklock : mat_locks)
    {
        omp_destroy_lock(&blocklock);
    }

    for(omp_lock_t& blocklock : rhs_locks)
    {
        omp_destroy_lock(&blocklock);
    }

    // Compress the global matrix to then compute the solution
    _system_mat.finalize();

    std::cout << "THAT" << std::endl;

    std::cout << "   \n" << std::endl;


    return;
}

void parallelSolver2::solve(const bool& print, std::ostream& out)
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

void parallelSolver2::process(const std::string & file_name)
{
    //this->_export(file_name);
    //this->_errorL2(true);
    //this->_errorH1(true);
    return;
}


const SparseMatrix<double>& parallelSolver2::getMat() const
{
    return _system_mat;
}

const VectorXd& parallelSolver2::getRHS() const
{
    return _rhs;
}

const VectorXd& parallelSolver2::getSol() const
{
    return _sol;
}


//  =========================================================
//          PROTECTED METHODS OF THE SOLVER CLASS
//  =========================================================


void  parallelSolver2::_computeStiff(FETools::SpectralFE<DIM>& _fe, std::array<omp_lock_t,NBLOCKS>& mat_locks)
{

    // Define the blocksize, using the NBLOCKS macro
    unsigned int blocksize = (_system_mat.rows()+ NBLOCKS -1)/NBLOCKS;
    // Loop over the elements of the mesh
    #pragma omp for
    for (const Element<DIM>& rect : _mesh.getElems())
    {

        // Update the element currently considered by the _fe object
        _fe.update_current(rect);


        // And compute the local stiffeness matrix
        
        MatrixXd LocStiff = ((_fe.J_cell_invT() * _fe.B_cell()).transpose()*(_fe.D_cell()*_fe.detJ())*(_fe.J_cell_invT() * _fe.B_cell()));
        
        for(unsigned int i = 0; i < LocStiff.rows(); i++)
        {
            unsigned int global_index_x = _dof.getMap()[i][_fe.getCurrent().getId()-1]-1;
            // Here compute the boundary flag
            unsigned short bound = _dof.getPoints()[global_index_x].getBound();
            // Compute the block index for the lock to be queried
            unsigned short l_i = global_index_x/blocksize; 
            //std::cout << "l_i: " << l_i<< std::endl;

            // Here evaluate the boundary condtion
            if(bound)
            {
                #pragma omp critical
                    std::cout << "Accessing element: " << global_index_x <<", "<< global_index_x<< std::endl;
                // Allow the write only if the corresponding block lock isn't owned by another thread: Set the lock
                omp_set_lock(&mat_locks[l_i]);
                    // Update the global matrix writing the 1 on the diagonal
                    _system_mat.coeffRef(global_index_x,global_index_x) = 1;
                // Unset the lock
                omp_unset_lock(&mat_locks[l_i]);
                
            }
            else
            {
                
                for(unsigned int j = 0; j <LocStiff.cols(); j++)
                {

                    unsigned int global_index_y = _dof.getMap()[j][_fe.getCurrent().getId()-1]-1;

                    // Write on _system_mat only if the value is different form 0: limit the writes on the global matrix
                    if(LocStiff.coeff(i,j) *_mu.value(_fe.quadrature_point(j,_dof)))
                    {

                        #pragma omp critical
                            std::cout << "Accessing element: " << global_index_x <<", "<< global_index_y<< std::endl;
                        // Allow the write only if the corresponding block lock isn't owned by another thread: Set the lock
                        omp_set_lock(&mat_locks[l_i]);
                            //Update the global matrix with non-zero contributions
                            _system_mat.coeffRef(global_index_x, global_index_y) += LocStiff.coeff(i,j) *_mu.value(_fe.quadrature_point(j,_dof));
                        // Unset the lock
                        omp_unset_lock(&mat_locks[l_i]);
                    }
                        
                }
                
            }
            
        }
        
    }



    return;
    
    
};


void  parallelSolver2::_computeMass(FETools::SpectralFE<DIM>& _fe, std::array<omp_lock_t,NBLOCKS>& mat_locks)
{

    // Define the blocksize, using the NBLOCKS macro
    unsigned int blocksize = (_system_mat.rows() + NBLOCKS - 1)/NBLOCKS;

    // Loop over the elements of the mesh
    #pragma omp for
    for (const Element<DIM>& rect : _mesh.getElems())
    {

        // Update the element currently considered by the _fe object
        _fe.update_current(rect);

        // As for the local mass, this will only bring about its contribution onto the digonal of the system matrix.
        // We compute it in the following way...

        // Loop over the quadrature points of the current element
        for(unsigned int j = 0; j < _fe.getNQ()[DIM-1]; j++)
        {
            for (unsigned int i = 0; i < _fe.getNQ()[0]; i++)
            {
                // And compute the local contributuions by means of the quadrature  weights
                unsigned int local_ind = (j)*_fe.getNQ()[0] + i;
                unsigned int global_index = _dof.getMap()[local_ind][_fe.getCurrent().getId()-1]-1;
                // Here compute the boundary flag
                unsigned short bound = _dof.getPoints()[_dof.getMap()[local_ind][_fe.getCurrent().getId()-1]-1].getBound();
                // Compute the block index for the lock to be queried
                unsigned short l_i = global_index/blocksize;

                // Allow for the update only if the lock is not already owned
                omp_set_lock(&mat_locks[l_i]);
                    _system_mat.coeffRef(global_index, global_index) += (bound == 0) * // WRITE THE CONTRIBUTION ONLY IF THE DOF IS NOT ON THE BOUNDARY
                                                            _fe.getQuad()[0].getW()[i] *
                                                            (1. / _fe.getJ().coeff(0, 0)) *
                                                            _fe.getQuad()[DIM - 1].getW()[j] *
                                                            (1. / _fe.getJ().coeff(1, 1)) *
                                                            _sigma.value(_fe.quadrature_point(local_ind, _dof));
                omp_unset_lock(&mat_locks[l_i]);
            } 
            
        }

    }
    
    return;
};


void parallelSolver2::_computeRHS(FETools::SpectralFE<DIM>& _fe, std::array<omp_lock_t,NBLOCKS>& rhs_locks)
{

    // Define the blocksize, using the NBLOCKS macro
    unsigned int blocksize = (_system_mat.rows() + NBLOCKS - 1)/NBLOCKS;

    // Loop over the elements of the mesh
    #pragma omp for
    for (const Element<DIM>& rect : _mesh.getElems())
    {

        // Update the element currently considered by the _fe object
        _fe.update_current(rect);

       
        // Loop over the quadrature points of the current element
        for(unsigned int j = 0; j < _fe.getNQ()[DIM-1]; j++)
        {
            for (unsigned int i = 0; i < _fe.getNQ()[0]; i++)
            {
                // And compute the local contribution by means of the quadrature weights
                unsigned int local_ind = (j)*_fe.getNQ()[0] + i;
                unsigned int global_index = _dof.getMap()[local_ind][_fe.getCurrent().getId() - 1] - 1;
                // Compute the block index for the lock to be queried
                unsigned short l_i = global_index/blocksize;

                // Allow for the update only if the lock is not already owend
                omp_set_lock(&rhs_locks[l_i]);
                // Updte the RHS vecotr
                _rhs[global_index] += _fe.getQuad()[0].getW()[i] *
                                            (1. / _fe.getJ().coeff(0, 0)) *
                                            _fe.getQuad()[DIM - 1].getW()[j] *
                                            (1. / _fe.getJ().coeff(1, 1)) *
                                            _f.value(_fe.quadrature_point(local_ind, _dof));
                // Unset lock after write
                omp_unset_lock(&rhs_locks[l_i]);
            }
            
        }
    
    }

    return;
    
}


void parallelSolver2::_apply_boundary(std::array<omp_lock_t,NBLOCKS>& rhs_locks)
{  

    // Define the blocksize, using the NBLOCKS macro
    unsigned int blocksize = (_system_mat.rows() + NBLOCKS - 1)/NBLOCKS;

    // We build a map that, for each boundary tag, stores the
    // corresponding boundary function.
    std::map<unsigned int, const Function<DIM> *> boundary_func;
    boundary_func[0]= &_g;
    boundary_func[1]= &_g;
    boundary_func[2]= &_g;
    boundary_func[3]= &_g;

    // Loop over the elements of the mesh
    #pragma omp for
    for (const Element<DIM>& rect : _mesh.getElems())
    {

        //inside this method we apply the direchelet b.c., known the value of the functions to be applied at the border
        //loop over all the nodes of the mesh given
        
        // get the dof for each element for the DoFHandler
        // and for each dof, see if it is on the boundary (and which boundary,
        // given the element considered is on the boundary too)
        for(unsigned int p = 0; p < _dof.dof_per_cell() && rect.getBound(); p++)
        {
            int gindex = _dof.getMap()[p][rect.getId()-1]-1;
            unsigned short bound = _dof.getPoints()[gindex].getBound();
            // Compute the block index for the lock to be queried
            unsigned short l_i = gindex/blocksize;

            // Allow for the write only if the lock is free
            omp_set_lock(&rhs_locks[l_i]);
            // For each node check that it is on the boundary (and which boundary...) and eventually substitute the value of the right hand side
            _rhs[gindex] = (bound==0)*_rhs[gindex] + 
                                (bound==1 || bound==30 ||bound == 29)*(boundary_func.at(0)->value(_dof.getPoints()[gindex]))+ 
                                (bound==2 )*(boundary_func.at(1)->value(_dof.getPoints()[gindex])) + 
                                (bound==3 || bound== 28 ||bound == 27)*(boundary_func.at(2)->value(_dof.getPoints()[gindex])) +
                                (bound==4)*(boundary_func.at(3)->value(_dof.getPoints()[gindex]));
            omp_unset_lock(&rhs_locks[l_i]);
            //finally set to 0 the rows corresponging to those nodes on the boundary, apart from elements
            //that represent the solution on that specific node....


            // Differetly form the serial case, this is implicitly done when assembling the global matrix, considering whether
            // a certain dof is on the border or not



        }

    }
     
    
    return;
};





