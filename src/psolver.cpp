#include "psolver.hpp"


void parallelSolver::setup(const unsigned short &option, const std::string& file_name)
{
    std::cout << "\nSetting up the mesh, initializing the DoF handler and algebraic structure...\n" <<std::endl;
    // As a first step, initialize the mesh, the ID array and the DoFhandler object
    if(option == 1)
    {
        _mesh.setMesh_csv(file_name,_dof.getDeg()[0]+1, _dof.getDeg()[DIM-1]+1);
        //_mesh.printMesh();
        _dof.genPoints(_mesh);
        //_dof.printPoints();
    }
        
        
    else if(option == 2)
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
    _rhs = VectorXd::Zero(size);
    _sol = VectorXd::Zero(size);
    std::cout << "Solver ready..."<<std::endl;
    return;
};


void parallelSolver::assemble()
{
    // generate threads of parallel threads
    #pragma omp parallel num_threads(4) 
    {

        // For every thread define maps that will store the elements 
        // to be inserted in the system matrix and rhs by each thread
        // For performance related aspects, the idea is to store in the map 
        // only the non-zero elements.

    

        std::unordered_map<std::pair<unsigned int, unsigned int>, double> mat_entries;
        std::unordered_map<unsigned int, double> rhs_entries;
        

        
        // Alternative implementation with vectors: in this case, the idea is
        // to add a new element to the vector for every non-zero update of an element of the matrix
        // We then loop over the vectors, reading the correspoing element of the global matrix to be uploaded.
        // Hopefully, since the vectors are data structures optimized for looping, the contiguous accesses in memory
        // will make up fo the presence of more elements

        
        SpectralFE<DIM> _fe(_dof.getDeg()[0]);

        // We build a map that, for each boundary tag, stores the
        // corresponding boundary function.
        std::map<unsigned int, const Function<DIM> *> boundary_func;
        boundary_func[0]= &_g;
        boundary_func[1]= &_g;
        boundary_func[2]= &_g;
        boundary_func[3]= &_g;

        // Store the global index of the elements considered by the thread in a private vector
        std::vector<unsigned int> local_elems(0);
        #pragma omp single
        {
            std::cout << "================================================================="<< std::endl;
            std::cout << " Starting parallel computation of system matrix and RHS vector...\n - "<< std::endl;
        }


        // Loop over the elements of the mesh
        #pragma omp for
        for (const Element<DIM>& rect : _mesh.getElems())
        {
            
            // Push the index of the element locally considered at the end of the local_elems vector
            local_elems.emplace_back(rect.getId());
            // Update the current element considered by the spectral _fe local solver
            _fe.update_current(rect);

            // Compute the local contribution to the global Stiffeness matrix
            this->_localStiff(_fe, mat_entries);
            
            // Compute the local contribution to the global Mass matrix
            this->_localMass(_fe,mat_entries);

            // Compute the local contribution to the global RHS vector
            this->_localRHS(_fe, rhs_entries);

            // Apply the boundary Dirichelet condition
            this->_apply_boundary(rect, rhs_entries, boundary_func);

            // After the parallel computational part, we insert the elements stored in the
            // unordered map serially into the global system and rhs vector
            
            // To do so, we use a pragma critical directive
        }

        #pragma omp single
        {
            std::cout << " System matrix and RHS vector computed.\n          -          \nStarting serial writing... \n - " << std::endl;

        }

        


        // Once we are done with all the computation, serially load the values of the global matrix stored inside
        // the local maps inside the global Eigen matrix

        #pragma omp critical
        {


            auto end = mat_entries.end();

            // To do so, loop over the locally evaluated elements of the mesh and access to the corresponding elements of the matrix
            for(const unsigned int& elem_id : local_elems)
            {
                // For each element loop over its quadrature points
                const Element<DIM>& rect = _mesh.getElems()[elem_id-1];
                //Loop over each quadrature point for each element
                for(unsigned int q = 0; q < rect.getNQ()[0]*rect.getNQ()[1]; ++q)
                {
                    unsigned int global_index_x = _dof.getMap()[q][rect.getId()-1]-1;
                    unsigned short bound = _dof.getPoints()[global_index_x].getBound();
                    

                    if(bound) // Check for boundary condition to avoid accessing 0 elements
                    {
                        // If element is on the boundary, just accesss the diagonal of the matrix
                        _system_mat.coeffRef(global_index_x,global_index_x) = mat_entries.at(std::make_pair(global_index_x, global_index_x));
                        
                    }
                    else
                    {

                        for(unsigned int p = 0; p < rect.getNQ()[0]*rect.getNQ()[1]; ++p)
                        {

                            unsigned int global_index_y = _dof.getMap()[p][rect.getId()-1]-1;
                            // Otherwise, access all the elements on the row
                            // To do so, first check if the element is actually present in the list
                            if(mat_entries.find(std::make_pair(global_index_x, global_index_y)) != end)
                            {
                                //std::cout << " Key: " << global_index_x << ", " << global_index_y << std::endl;
                                _system_mat.coeffRef(global_index_x,global_index_y) = mat_entries.at(std::make_pair(global_index_x, global_index_y));
                            
                            
                            }
                            
                        }

                    }

                    // Finally fill the rhs vector
                    
                    _rhs[global_index_x] = rhs_entries.at(global_index_x);

                }
            }
        }
        #pragma omp barrier

        #pragma omp single
        {
            std::cout << " Serial writing of the system matrix and RHS vector completed."<< std::endl;
            std::cout << "================================================================="<< std::endl;
            std::cout << "System matix: \n" << std::endl;
            //std::cout << _system_mat << std::endl;
        }

        
        
    }

    return;
}

void parallelSolver::solve(const bool& print, std::ostream& out)
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

void parallelSolver::process(const std::string & file_name)
{
    //this->_export(file_name);
    //this->_errorL2(true);
    //this->_errorH1(true);
    return;
}



//  =========================================================
//          PROTECTED METHODS OF THE SOLVER CLASS
//  =========================================================


void  parallelSolver::_localStiff(FETools::SpectralFE<DIM>& _fe, std::unordered_map<std::pair<unsigned int, unsigned int>, double>& mat_entries)
{
    // And compute the local stiffeness matrix
    
    // Compute the local matrix and compress it onto the global matrix using the ID array (_dof.getMap())
    MatrixXd LocStiff = ((_fe.J_cell_invT() * _fe.B_cell()).transpose()*(_fe.D_cell()*_fe.detJ())*(_fe.J_cell_invT() * _fe.B_cell()));
    
    for(unsigned int i = 0; i < LocStiff.rows(); i++)
    {

        // Here compute the boundary flag
        unsigned short bound = _dof.getPoints()[_dof.getMap()[i][_fe.getCurrent().getId()-1]-1].getBound();

        // Here evaluate the boundary condtion
        if(bound)
        {
            // Update the unordered map just by inserting the diagonal element equal to 1
            mat_entries[std::make_pair(_dof.getMap()[i][_fe.getCurrent().getId()-1]-1, _dof.getMap()[i][_fe.getCurrent().getId()-1]-1)]= 1;
            
        }
        else
        {
            
            for(unsigned int j = 0; j <LocStiff.cols(); j++)
            {
                
                // Each thread writes onto its locally defined unordered_map for the entries of the matrix but...
                // diffently form the serial code, instead of later erasing the elments of the matrix coresponding to 
                // the dof on the boundary in the boundary private method, here we simply avoid the write on those elements
                // by checking if is on the border before-hand.
                std::pair<unsigned int, unsigned int> key = std::make_pair(_dof.getMap()[i][_fe.getCurrent().getId()-1]-1, _dof.getMap()[j][_fe.getCurrent().getId()-1]-1);
                // Write on map only if the value is different form 0
                
                if(LocStiff.coeff(i,j) *_mu.value(_fe.quadrature_point(j,_dof)))
                {
                    mat_entries[key] += LocStiff.coeff(i,j) *_mu.value(_fe.quadrature_point(j,_dof));
                }
                    
            }

            
            
        }

        
        
    }



    return;
    
    
};

void  parallelSolver::_localMass(FETools::SpectralFE<DIM>& _fe, std::unordered_map<std::pair<unsigned int, unsigned int>, double>& mat_entries)
{

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

            // Each thread writes onto its locally defined unordered_map for the entries of the matrix
            std::pair<unsigned int, unsigned int> key = std::make_pair(global_index,global_index);
            mat_entries[key] += (bound==0) * // WRITE THE CONTRIBUTION ONLY IF THE DOF IS NOT ON THE BOUNDARY
                                _fe.getQuad()[0].getW()[i] *
                                (1. / _fe.getJ().coeff(0, 0)) *
                                _fe.getQuad()[DIM - 1].getW()[j] *
                                (1. / _fe.getJ().coeff(1, 1)) *
                                _sigma.value(_fe.quadrature_point(local_ind, _dof));

        }
        
        
    }
    
    return;
};


void parallelSolver::_localRHS(FETools::SpectralFE<DIM>& _fe, std::unordered_map<unsigned int, double>& rhs_entries)
{

    
       
    // Loop over the quadrature points of the current element
    for(unsigned int j = 0; j < _fe.getNQ()[DIM-1]; j++)
    {
        for (unsigned int i = 0; i < _fe.getNQ()[0]; i++)
        {
            // And compute the local contribution by means of the quadrature weights
            unsigned int local_ind = (j)*_fe.getNQ()[0] + i;
            unsigned int global_index = _dof.getMap()[local_ind][_fe.getCurrent().getId() - 1] - 1;
            // Each thread writes onto its locally defined unordered_map for the entries of the rhs vector
            rhs_entries[global_index] += _fe.getQuad()[0].getW()[i] *
                                        (1. / _fe.getJ().coeff(0, 0)) *
                                        _fe.getQuad()[DIM - 1].getW()[j] *
                                        (1. / _fe.getJ().coeff(1, 1)) *
                                        _f.value(_fe.quadrature_point(local_ind, _dof));
        }
        
        
    }

    return;
    
}


void parallelSolver::_apply_boundary(const Element<DIM>& rect, std::unordered_map<unsigned int, double>& rhs_entries, const std::map<unsigned int, const Function<DIM> *>& boundary_func)
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

        //for each node check that it is on the boundary (and which boundary...) and eventually substitute the value of the right hand side

        // In the parallel version of the code, the changes committed by the threads are first updated on the locally owned
        // maps and will then reduced on to the _rhs vector consequently
        rhs_entries[gindex] = (bound==0)*_rhs[gindex] + 
                            (bound==1 || bound==30 ||bound == 29)*(boundary_func.at(0)->value(_dof.getPoints()[gindex]))+ 
                            (bound==2 )*(boundary_func.at(1)->value(_dof.getPoints()[gindex])) + 
                            (bound==3 || bound== 28 ||bound == 27)*(boundary_func.at(2)->value(_dof.getPoints()[gindex])) +
                            (bound==4)*(boundary_func.at(3)->value(_dof.getPoints()[gindex]));
        //finally set to 0 the rows corresponging to those nodes on the boundary, apart from elements
        //that represent the solution on that specific node....

        

    }
     
    
    return;
};



// double parallelSolver::_errorL2(const bool& print, std::ostream& out) const
// {
//     VectorXd error(_sol.rows());
//     double sum(0);

//     // loop over all the elements
//     for(const Element<DIM>& e: _mesh.getElems())
//     {
//         //loop over all the quadrature points of the element
//         for(unsigned int j = 0; j < _fe.getNQ()[DIM-1]; j++)
//         {
//             for (unsigned int i = 0; i < _fe.getNQ()[0]; i++)
//             {
//                 // determine the index of local quadrature point currently considered
//                 unsigned int q = (j)*_fe.getNQ()[0]+i;
//                 // and the global index of the node currently considered
//                 unsigned int gindex = _dof.getMap()[q][e.getId()-1]-1;

//                 // compute the local weight
//                 double weight = (_fe.getQuad()[0].getW()[i] *
//                                  (1. / _fe.getJ().coeff(0, 0)) *
//                                  _fe.getQuad()[DIM - 1].getW()[j] *
//                                  (1. / _fe.getJ().coeff(1, 1)));

//                 sum += (_e.value(_dof.getPoints()[gindex]) - _sol.coeff(gindex))*(_e.value(_dof.getPoints()[gindex]) - _sol.coeff(gindex))* weight;
//             }
            
            
//         }

//     }

     
//     if(print)
//     {
//         out << "Norm L2 of error computed: " << std::scientific << std::sqrt(sum) << std::endl;
//     }
        
//     return std::sqrt(sum);
// }



// double parallelSolver::_errorH1(const bool& print,std::ostream& out)
// {
//     // first define the temp variable to store the sum of 

//     double sum(0);
    

//    //  Loop over each element 
//    for(const Element<DIM>& elem : _mesh.getElems())
//    {
//         // For each element get the evaluation of the exact solution at the quadrature nodes
//         // and compute the numerical derivative using the derivative matrix of the element computed by the fe object
//         _fe.update_current(elem);
//         //  An eigen matrix to store the values of the numerical solution on the local quadrature nodes
//         // necessary to then compute the evaluation of the derivatives
//         MatrixXd un(_fe.getNQ()[0],_fe.getNQ()[DIM-1]);
        
//         // Loop over the quadrature points of the element
//         for(unsigned int j = 0; j < _fe.getNQ()[DIM-1]; j++)
//         {
//             for (unsigned int i = 0; i < _fe.getNQ()[0]; i++)
//             {
//                 // determine the index of local quadrature point currently considered
//                 unsigned int q = (j)*_fe.getNQ()[0]+i;
//                 // and the global index of the node currently considered
//                 unsigned int gindex = _dof.getMap()[q][elem.getId()-1]-1;
                
//                 // get the exact solution on the quadrature nodes and place its values in a matrix
//                 un(i,j) = _sol.coeff(gindex);
                
//             }
//         }

//         // std::cout << "Un"<<std::endl;
//         // std::cout <<un <<"\n"<< std::endl;

//         //Compute the numerical derivative
//         MatrixXd ux = (_fe.getJ().coeff(0, 0))*(_fe.getB()[0]* un);
//         MatrixXd uy = (_fe.getJ().coeff(1, 1))*(_fe.getB()[0]* un.transpose()).transpose();

//         // std::cout << "Ux"<<std::endl;
//         // std::cout <<ux <<"\n"<< std::endl;
//         // std::cout << "Uy"<<std::endl;
//         // std::cout <<uy <<"\n"<< std::endl;

//         // Loop over the quadrature points of the element to sum all the differences over the quadrature points
//         for(unsigned int j = 0; j < _fe.getNQ()[DIM-1]; j++)
//         {
//             for (unsigned int i = 0; i < _fe.getNQ()[0]; i++)
//             {
//                 // determine the index of local quadrature point currently considered
//                 unsigned int q = (j)*_fe.getNQ()[0]+i;
//                 // and the global index of the node currently considered
//                 unsigned int gindex = _dof.getMap()[q][elem.getId()-1]-1;
                
//                 // sum all the contributions for the element currenlty considered
//                 sum += ((un.coeff(i, j) - _e.value(_fe.quadrature_point(q, _dof))) *
//                             (un.coeff(i, j) - _e.value(_fe.quadrature_point(q, _dof))) + // (U - U_n)^2
//                         (ux.coeff(i, j) - _e.grad(_fe.quadrature_point(q, _dof))[0]) *
//                             (ux.coeff(i, j) - _e.grad(_fe.quadrature_point(q, _dof))[0]) + // (dUx -dU_nx)^2
//                         (uy.coeff(i, j) - _e.grad(_fe.quadrature_point(q, _dof))[1]) *
//                             (uy.coeff(i, j) - _e.grad(_fe.quadrature_point(q, _dof))[1])) *
//                        (_fe.getQuad()[0].getW()[i] *
//                         (1. / _fe.getJ().coeff(0, 0)) *
//                         _fe.getQuad()[DIM - 1].getW()[j] *
//                         (1. / _fe.getJ().coeff(1, 1))); // (dUy -dU_ny)^2
//             }
//         }
//    }

//    if(print)
//         out << "Norm H1 of error computed: " << std::scientific <<std::sqrt(sum) << std::endl;

//     return  std::sqrt(sum);
    

// }





