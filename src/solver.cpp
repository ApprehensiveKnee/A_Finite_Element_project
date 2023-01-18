#include "solver.hpp"


void serialSolver::setup(unsigned int const &option)
{
    switch (option)
    {
    case 1:
        _mesh.setMesh_csv();
        _map =_mesh.indexMapping(_fe.getDeg(), _mesh.get_nEx(), _fe.getDeg(), _mesh.get_nEy());
        
        // ------------ PRINTING SECTION ----------
        // for (unsigned int j = 0; j < _map[0].size(); j++)
        // {
        //     for (unsigned int i = 0; i < _map.size(); i++)
        //     {
        //         std::cout << _map[i][j] << " " << std::endl;
        //     }
        //     std::cout<< "-------------------" << std::endl;
        // }
        // -----------------------------------------
        
        break;
    case 2:
        if constexpr (DIM == 1)
        {
            //some random data to generate the mesh
            //change to fit the mesh to specific needs
            _mesh.genMesh(0,100,50);
            _map = _mesh.indexMapping(_fe.getDeg(), _mesh.get_nEx());
        }
        else if constexpr (DIM == 2)
        {
            _mesh.genMesh(0,100,0,100,50,50);
            _map=_mesh.indexMapping(_fe.getDeg(), _mesh.get_nEx(), _fe.getDeg(), _mesh.get_nEy());
        }
            
        
        break;
    default:
        break;
    }
    unsigned int size = (_mesh.get_nEx()*_fe.getDeg()+1)*(_mesh.get_nEy()*_fe.getDeg()+1);
    _Stiffness.resize(size, size);
    _Mass.resize(size, size);
    _rhs.resize(_Stiffness.rows());
    _sol.resize(_Stiffness.rows());
    return;
};


void serialSolver::assemble()
{
    this->_calculateStiff();
    this->_calculateMass();
    this->_calculateRHS();
    //this->_apply_boundary();
    // moved to solve method as to facilitate the application of boundary conditions
    return;
}

void serialSolver::solve()
{
    // to solve the system we use the sovers provided by Eigen;
    // as reference here we use the BiCGSTAB method:

    //VectorXd xe = VectorXd::Constant(mat.rows(), 1);                // define exact solution
    //VectorXd b = mat*xe;
    //cout << "Left hand side"<< endl;                                          // compute right-hand side
    //cout << b << endl;

    SparseMatrix<double> mat = _Stiffness + _Mass;
    // -----> update .OSS: ideally we could use a single matrix as member of the class solver
    this->_apply_boundary(mat);
    std::cout << "The rhs is:\n" << _rhs << std::endl;
    
    BiCGSTAB<SparseMatrix<double>> solver;
    solver.compute(mat);

    _sol = solver.solve(_rhs);

    std::cout << "#iterations: " << solver.iterations() << std::endl;
    std::cout << "The solution is:"  << std::endl;                                 // solving
    std::cout << _sol << std::endl;   
    
    // double relative_error = (x - xe).norm() / xe.norm(); // norm() is L2 norm
    // std::cout << "The relative error is:\n" << relative_error << std::endl;                                           // display solution
    return;
}

void serialSolver::process()
{
    this->_export();
    //this->_plot();
     this->_error();
    return;
}


//  --------------------------------------------------

MatrixXd serialSolver::_LocStiff() const
{
    MatrixXd LocStiff(_fe.dof_per_cell(),_fe.dof_per_cell());
    //the dimension of the local matrix depends on the
    //degree of the _fe space used
    unsigned int n_q = _fe.getQPoints().size();
    for(unsigned int i = 0; i<_fe.dof_per_cell(); ++i )
    {
        
        for(unsigned int j = 0; j < _fe.dof_per_cell() ; ++j)
        {
            for(unsigned int q = 0; q < n_q; ++q)
            {
                //std::cout << "Shape grad:"<<<<_fe.shape_grad(i,q);
                LocStiff(i,j) += _mu.value(_fe.quadrature_point(q)) *
                                 _fe.shape_grad(j,q).transpose() *
                                 _fe.shape_grad(i,q) *
                                 _fe.JxW(q);
                                
            }
        }
    }

    return LocStiff;
};

MatrixXd serialSolver::_LocMass() const
{
    MatrixXd LocMass(_fe.dof_per_cell(),_fe.dof_per_cell());
    //the dimension of the local matrix depends on the
    //degree of the _fe space used
    unsigned int n_q = _fe.getQPoints().size();
    for(unsigned int i = 0; i<_fe.dof_per_cell() ; ++i )
    {
        for(unsigned int j = 0; j < _fe.dof_per_cell() ; ++j)
        {
            for(unsigned int q = 0; q < n_q; ++q)
            {
                LocMass(i,j) += _sigma.value(_fe.quadrature_point(q)) *
                                _fe.shape_value(j,q) *
                                _fe.shape_value(i,q) *
                                _fe.JxW(q);
                                
            }
        }
    }

    return LocMass;
};


VectorXd serialSolver::_LocRHS() const
{
    VectorXd LocRHS(_fe.dof_per_cell());
    //the dimension of the local rhs depends on the
    //degree of the _fe space used
    unsigned int n_q = _fe.getQPoints().size();
    for(unsigned int i = 0; i<_fe.dof_per_cell(); ++i )
    {
        for(unsigned int q = 0; q < n_q; ++q)
        {
            LocRHS[i] +=    _f.value(_fe.quadrature_point(q)) *
                            _fe.shape_value(i,q) *
                            _fe.JxW(q);
                            
        }
}

    return LocRHS;
};

void  serialSolver::_calculateStiff()
{
    //first set all elements of the matrix to 0
    _Stiffness.setZero();
    MatrixXd temp;
    //loop over the elements of the mesh
    for (Element_2D rect : _mesh.getElems())
    {
        //update the current element considered by the spectral _fe local solver
        _fe.update_current(rect);
        //and compute the local stiffeness matrix
        temp = this->_LocStiff();
        for (unsigned int i = 0; i< temp.rows(); ++i)
        {
            for(unsigned int j = 0; j< temp.cols(); ++j)
            {
                if(!std::isnan(temp(i,j)))
                _Stiffness.coeffRef(_map[i][rect.getId()-1]-1,_map[j][rect.getId()-1]-1) +=temp(i,j);
            }
        }
        

    }

    std::cout << "Stiffness : COMPUTED" << std::endl;

    return;
    
    
};

void  serialSolver::_calculateMass()
{
    //first set all elements of the matrix to 0
    _Stiffness.setZero();
    MatrixXd temp;

    //loop over the elements of the mesh
    for (Element_2D rect : _mesh.getElems())
    {
        //update the current element considered by the spectral _fe local solver
        _fe.update_current(rect);
        //and compute the local mass matrix
        temp = this->_LocMass();
        //finally update the elements of the global mass matrix
        for (unsigned int i = 0; i< temp.rows(); ++i)
        {
            for(unsigned int j = 0; j< temp.cols(); ++j)
            {
                if(!std::isnan(temp(i,j)))
                _Mass.coeffRef(_map[i][rect.getId()-1]-1,_map[j][rect.getId()-1]-1) +=temp(i,j);
            }
        }
        

    }
    std::cout <<"Mass matrix is:" << _Mass << std::endl;
    
    std::cout << "Mass : COMPUTED" << std::endl;
    return;
};


void serialSolver::_calculateRHS()
{
    //first set all the elements of the vector to 0
    _rhs.setZero();
    VectorXd temp;

    //loop over all the elements of the mesh
    for (Element_2D rect : _mesh.getElems())
    {
        //update the current element considered by the spectral _fe local solver
        _fe.update_current(rect);
        //and compute the local RHS vector
        temp = this->_LocRHS();
        //update the elements of the global mesh
        for (unsigned int i = 0; i< temp.cols(); ++i)
        {
            if(!std::isnan(temp[i]))
            _rhs[_map[i][rect.getId()-1]-1] +=temp[i];
        }

    }
    std::cout << "rhs : COMPUTED" << std::endl;
    return;
    
}


void serialSolver::_apply_boundary( SparseMatrix<double> &mat)
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
        if(_mesh.onBorder(rect))
        {
            _fe.update_current(rect);
            for(unsigned int i = 0; i<_fe.dof_per_cell(); ++i)
            {
                int index = _map[i][rect.getId()-1]-1;
                unsigned short temp = _mesh.onBorder(_fe.getINodes()[i]);
                //for each node check that it is on the boundary (and which boundary...) and eventually substitute the value of the right hand side
                _rhs[index] = (temp==0)*_rhs[_map[i][rect.getId()-1]-1] + 
                                    (temp==1)*(boundary_func[0]->value(_fe.getINodes()[i])) + 
                                    (temp==2)*(boundary_func[1]->value(_fe.getINodes()[i])) + 
                                    (temp==3)*(boundary_func[2]->value(_fe.getINodes()[i])) +
                                    (temp==4)*(boundary_func[3]->value(_fe.getINodes()[i]));
                //finally set to 0 the rows corresponging to those nodes on the boundary, apart from elements
                //that represent the solution on that specific node....
                mat.prune([&index](int i, int j, float) { return i!=index; });
                mat.coeffRef(index,index) = 1;
                
            }
        }
        

    } 
    std::cout << "boundary : APPLIED" << std::endl;
    return;
};

/*
void serialSolver::_plot()
{
    //use gnuplot_iostream
    Gnuplot pipe;
    double x_min = 0, x_max = 1, y_min = 0, y_max = 1;

    // Set the plot style to a surface plot
    pipe << "set pm3d\n";
    pipe << "set view map\n";
    pipe << "set xrange [" << x_min << ":" << x_max << "]\n";
    pipe << "set yrange [" << y_min << ":" << y_max << "]\n";
    pipe << "set xlabel 'x axis' offset -3,-2\n";
    pipe << "set ylabel 'y axis' offset 3,-2\n";
    pipe << "set palette rgbformulae 33,13,10\n";

    std::vector<std::tuple<double, double, double>> points;
    for(unsigned int i = 0; i <_mesh.get_nNodes(); ++i)
    {
        points.emplace_back({_mesh.getNode(i).getX(), _mesh.getNode(i).getY(), _sol[i]});
    }
    
    // plot the function and store the visual representation in a png file
    pipe << "set term png\n";
    pipe << "set output 'solution.png'\n";
    pipe << "splot '-' with points title 'Visualization of the numerical solution'\n";
    pipe.send(xy_pts);
    pipe << "e\n";
    pipe << "pause mouse close\n";
    return;
}

*/


void serialSolver::_export()
{
    std::string matrixFileOut("./solution.mtx");
    Eigen::saveMarket(_sol, matrixFileOut);
    return;
}

void serialSolver::_error(std::ostream& out)
{
    VectorXd error(_sol.rows());
    for(Element_2D elem : _mesh.getElems())
    {
        _fe.update_current(elem);
        for (unsigned int node = 0; node < _fe.dof_per_cell(); node++)
        {
            //compute the difference as the norm of the difference vector obtained between the solution vector numerically obtained with the FEM method
            //and the vector obtained valuating the exact solution on the nodes of the mesh (on which we computed the numerical solution)
            error[_map[node][elem.getId()-1]-1] = _e.value(_fe.quadrature_point(node)) - _sol[_map[node][elem.getId()-1]-1]; 
        }
        
        
    }

   out << "Norm of the error computed " << error.norm() << std::endl;
   return;
}


// void serialSolver::dummy(const Element_2D & rect)
// {
//     std::cout <<_fe.getPun() << std::endl;
//     std::cout << _fe.getDeg() << std::endl;
//     std::cout << _fe.getCurrent().getId() << std::endl;
//     _fe.update_current(rect);
// }