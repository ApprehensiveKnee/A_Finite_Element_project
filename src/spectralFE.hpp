

#ifndef FE
#define FE

#include <numeric>
#include <algorithm>
#include <functional>
#include <tuple>
#include <type_traits>
#include "elements.hpp"
#include "mesh.hpp"
#include "mat_utilities.hpp"

namespace FETools
{


    template <class ElementType>
    class SpectralFE
    {
    public:

        //  ==========================================================================
        // CUSTOM CLASSSES TO IMPLEMENT MATRIX FREE MATRIX PRODUCT ON THE CELL MATRIXES
        //  ==========================================================================
        // since Bref_cell,J_cell and D_cell are made up of copies of the same values repeated in certain positions,
        // this lets us compute the matrix products in a more efficent way simply by accessing the Jacobian of the element
        // D_ref, e Bx_ref/ By_ref
        

        //  Class for Bref_cell extention form Bx_ref and By_ref
        template<typename Derived>
        class ExtentB: public Eigen::EigenBase<Derived>
        {
        public:
            ExtentB() = default;
            ExtentB(const std::array<SparseMatrix<double>,DIM>& B, const std::array<unsigned int, DIM>& nq, const unsigned int& re = r)
            : _B(B),
            _nq(nq),
            _r(re)
            {}

            int rows() const
            {
                return DIM * _nq[0] * (DIM==2?_nq[DIM-1]:1)/*DIM * n_q*/;
            }

            int cols() const
            {
                // REMARK:
                // dofs per cell and number of quadrature points are supposeed to be the same, still
                // this implementation better explains how the B matrix is built (even though it may take slightly more time and resources)
                return DoFHandler(_r).dof_per_cell() /* dof_per_cell */;
            }

            // Return sub-matrix in each iteration
            double operator()(const int& i/* nqx*nqy*DIM*/, const int& j/*nqx*nqy*/) const
            {
                // REMARK
                // Here again, we could have just used a single variable, yet for better clarity, the choice was
                // to "keep the dimensions separated"
                const unsigned int dofx = DoFHandler(_r).getDeg()[0] +1;
                const unsigned int dofy = DoFHandler(_r).getDeg()[1] +1 ;

                if constexpr(DIM == 2)
                {
                    // lines of Bx_ref
                    if(i%2 == 0 /*even rows (x dedicated)*/)
                    {
                        if((i>>1)/_nq[0] == j/dofx) // i and j indexing on an element on the "diagonal" block
                            return _B[0].coeff((i/2)%_nq[0],j%dofx);
                        else 
                            return 0.;
                    }
                    // lines of By_ref
                    else /*odd rows*/
                    {
                        if((i>>1)%_nq[1] == j%dofy) // i and j indexing on an elemente on the diagonal block
                            return _B[1].coeff((i>>1)/_nq[1], j/dofy);
                        else 
                            return 0.;
                    }

                }
                else
                {
                    if(i/_nq[0] == j/dofx) // i and j indexing on an element on the "diagonal" block
                            return _B[0].coeff(i%_nq[0],j%dofx);
                        else 
                            return 0.;
                } 
                    
            }

            double transpose(const int& i/* nqx*nqy*DIM*/, const int& j/*nqx*nqy*/) const
            {
                return this->operator()(j,i);
                    
            }


        private:
            const std::array<SparseMatrix<double>,DIM>& _B;
            const std::array<unsigned int,DIM>& _nq;
            const unsigned int& _r;
        };

        //  Class for D_cell extension form D_ref matrix
        
        template<typename Derived>
        class ExtentD: public Eigen::EigenBase<Derived>
        {
        public:
            ExtentD() = default;
            ExtentD(const SparseMatrix<double> &A)
            : _A(A)
            {}


            int rows() const
            {
                return _A.rows() * DIM;
            }

            int cols() const
            {
                return _A.cols() * DIM;
            }

            // Return sub-matrix in each iteration
            double operator()(const int& i, const int& j) const
            {
                if( i == j)
                    return _A.coeff(i/DIM, j/DIM);
                else
                    return 0.;
            }

            // Operator * overloading to support multiplication with custom expression class (coeff * D_ref)
            SparseMatrix<double> operator*(const double& coeff) const
            {
                SparseMatrix<double> result(this->rows(), this->rows());

                for (int i = 0; i < this->rows(); ++i)
                        result.coeffRef(i, i) += coeff * this->operator()(i, i);

                return result;
            }

        private:
            const SparseMatrix<double>& _A;
        };

        //  Class for J_cell extention from Jacobian of the element
        template<typename Derived>
        class ExtentJ: public Eigen::EigenBase<Derived>
        {
        public:
            ExtentJ() = default;
            ExtentJ(const MatrixXd &A, const unsigned int& nq)
            : _A(A),
            _nq(nq)
            {}

            int rows() const
            {
                return _A.rows()*_nq;
            }

            int cols() const
            {
                return _A.cols()*_nq;
            }

            // Return sub-matrix in each iteration
            double operator()(const int& i, const int& j) const
            {

                if constexpr(DIM == 2)
                {

                    if(i%2 == 0 /* row is even*/)
                    {
                        if(j%2 == 0 /*col is even*/ && i == j /*on the diagonal*/ )
                        {
                            return _A(0, 0);
                        }
                        else if(j%2 == 1 /*col is */ && j == i+1 /*on the upper diagonal*/)
                        {
                            return _A(0,1);
                        }
                        else
                            return 0.;
                    }
                    else /*row is odd*/
                    {
                        if(j%2 == 0 /*col is even*/ && j == i-1 /*on the lower diagonal*/)
                        {
                            return _A(1, 0);
                        }
                        else if(j%2 == 1/*col is odd*/ && j == i/*on the diaglona*/)
                        {
                            return _A(1,1);
                        }
                        else
                            return 0.;
                    }

                }
                else
                {
                    if(i == j /*on the diagonal*/)
                        return _A(0,0);
                }
                    
            }

            // Operator * overloading to support multiplication with custom expression class (J_cell^invT * B_cell)
            MatrixXd operator*(const ExtentB<SparseMatrix<double>> &B) const
            {

                SparseMatrix<double> result(this->rows(), B.cols());
                for (int i = 0; i <this->rows(); ++i)
                {
                    for (int j = 0; j < B.cols(); ++j)
                    {
                        for (int k = 0; k < this->cols(); ++k)
                            result.coeffRef(i, j) += this->operator()(i, k) * B(k, j);
                    }

                }
                return result;
            }


            double transpose(const int& i, const int& j) const
            {
                return this->operator()(j,i);
            }


        private:
            const MatrixXd& _A;
            const unsigned int _nq;
        };

        

        //  ==========================================================================




        //constructor
        SpectralFE(const unsigned int &re = r,
                const unsigned int &nqx = r + 1,
                const unsigned int &nqy = r + 1)
            : _current_elem(),
            _qr(),
            _r(re),
            _nq(classInit(nqx, nqy)),
            _quad_points(),
            _D_ref(nqx * DIM==2?nqy:1 /* n_q */, nqx * DIM==2?nqy:1 /* n_q */),
            _B_ref(classInit(nqx, nqy, r)), // dof_per_cell = n_q
            _J_invT(DIM, DIM),
            _Dcell(_D_ref),
            _Jcell_invT(_J_invT, _nq[0]*((DIM==2)?_nq[DIM-1]:1)),
            _Bcell(_B_ref, _nq, _r),
            pun(0){};

        SpectralFE(const ElementType &element,
                const unsigned int &re = r,
                const unsigned int &nqx = r + 1,
                const unsigned int &nqy = r + 1)
            : _current_elem(element),
            _qr(),
            _r(r),
            _nq(classInit(nqx, nqy, r)),
            _quad_points(),
            _D_ref(nqx * DIM == 2 ? nqy : 1 /* n_q */, nqx * DIM == 2 ? nqy : 1 /* n_q */),
            _B_ref(classInit(nqx, nqy)), // dof_per_cell = n_q
            _J_invT(DIM, DIM),
            _Dcell(_D_ref),
            _Jcell_invT(_J_invT, _nq[0]*(DIM==2)?_nq[DIM-1]:1),
            _Bcell(_B_ref, _nq, _r),
            pun(0){};
        //standar getters
        const ElementType& getCurrent() const
        {
            return _current_elem;
        };

        const std::array<Quadrature, DIM>& getQuad() const
        {
            return _qr;
        };

        const unsigned int& getDeg() const
        {
            return _r;
        };

        const std::array<unsigned int, DIM>& getNQ() const
        {
            return _nq;
        }

        void setNQ(const std::array<unsigned int,DIM>& exactness)
        {
            _nq = exactness;
            return;
        }

        const std::vector<Point>& getQPoints() const
        {
            return _quad_points;
        };

        const SparseMatrix<double>& getD() const
        {
            return _D_ref;
        }

        const std::array<SparseMatrix<double>, DIM>& getB( ) const
        {
            return _B_ref;
        }

        const MatrixXd& getJ() const
        {
            return _J_invT;
        } 

        const ExtentD<SparseMatrix<double>>& D_cell() const
        {
            return _Dcell;
        }

        const ExtentJ<MatrixXd>& J_cell_invT() const
        {
            return _Jcell_invT;
        }

        const ExtentB<SparseMatrix<double>>& B_cell() const
        {
            return _Bcell;
        }

        double detJ() const
        {
            if constexpr(DIM == 2)
                return _current_elem.jacobian().determinant();
            else
                return _current_elem.jacobian();
            
        }

        const unsigned int& getPun() const
        {
            return pun;
        };

        // a method to move up on the examination and computation of local matrixes
        // updating the element of the mesh to be computed
        void move_up(const Mesh<ElementType> &mesh)
        {
            if(this->getPun() < (mesh.get_nElems()))
            {
                if constexpr(DIM == 2)
                {
                    if(std::is_same<ElementType, Element_2D>::value)
                        this->update_current(mesh.getElem(pun));
                    else
                    {
                        throw std::runtime_error("Something went wrong, the dimensonal paramether doesn't correspond with the type of mesh element\n");
                    }
                }
                else
                {
                    if(std::is_same<ElementType, Element_1D>::value)
                        this->update_current(mesh.getElem(pun));
                    else
                    {
                        throw std::runtime_error("Something went wrong, the dimensonal paramether doesn't correspond with the type of mesh element\n");
                    }

                }

                pun++;
                    
            }
            else 
            {
                std::cout << "The SpectralFE object has reached the end of the mesh"<<std::endl;
                pun = -1;
            }
            return;
        };

        // a method to update the element currently considered (used in move_up method)
        void update_current( const ElementType &geoele )
        {
            _current_elem = geoele;
            this->_update_quad();
            // update the Jacobian
            this->_update_J();
            if(_current_elem.getNQ() == this->getNQ())
            {
                return;
            }
            else
            {
                //update the degree of exactenss for the quadrature on the current element
                this->setNQ(_current_elem.getNQ());
                // update the evaluation of quadrature points
                this->_update_D();
                if constexpr (DIM ==2)
                {
                    
                    
                    // In this case, recompute Bx
                    if(_current_elem.getNQ()[0] != this->getNQ()[0])
                    {
                        this->_update_B(0);
                    }
                    // In this case, recompute By
                    if(_current_elem.getNQ()[1] != this->getNQ()[1])
                    {
                        this->_update_B(1);
                    }
                    //update quadrature nodes and weights, as well as the degree

                }

            }
            
            return;
        };

        // a method to update first set the reference matrixes
        void set()
        {
            // First initialization of the reference matrixes
            this->_update_D();
            this->_update_B(0);
            this->_update_B(1);
            return;
        }
        
        // member function to map the quadrature points on to the specific element
        Point quadrature_point(const unsigned int &q,const DoFHandler& dof) const
        {

            //  FIRST APPROACH:
            // compute the quadrature point usinge the direct map function with the quadrature pointson the reference element

            // if constexpr (DIM == 2)
            // {
            //     double x, y;
            //     std::tie(x,y) = _current_elem.direct_map((this->getQPoints()[q]).getX(),(this->getQPoints()[q]).getY());
            //     Point my_point(x,y);
            //     return my_point;
            // }
            // else
            // {
            //     double x;
            //     x = _current_elem.direct_map((this->getQPoints()[q]).getX());
            //     Point my_point(x);
            //     return my_point;
            // }

            //  SECOND APPROACH
            // get the quarature point for the global mesh computed by a DoFHandler object ALREADY INITIALISED
            
            // using the index of the element currently considered, get the coordinates of the quadrature points
            
            return dof.getPoints()[dof.getMap()[q][this->getCurrent().getId()-1]-1];

            // this second approach requires the elements to each have the same number of quadarture points, but,
            // on the other hand, prevents us from performing refìdundant computation
            
            
        };

        
        //default destructor
        ~SpectralFE() = default;

    private:
        
        //firstly, a method to compute the coordinates and  weithts of the quadrature points along a certain dimension
        //used to approximate integrals on the current element
        std::vector<double> _comp_quad_c(const unsigned short& dir)
        {
            _qr[dir].LGL_quadratures(this->getNQ()[dir]);
            return _qr[dir].getN();
        };

        std::vector<double> _comp_quad_w(const unsigned short& dir)
        {
            
            _qr[dir].LGL_quadratures(this->getNQ()[dir]);
            return _qr[dir].getW();
        };

        // a method to update the vector of quadrature points for the current element
        void _update_quad()
        {
            // delete the quadrature coordinates of the past element
            // OSS. in this case, the degree of the FE spaces along the two dimensions
            // are supposed to be the same and identical to the _r member. Furthermore, 
            //for each element, the number of quadrature points is supposed to be the same:
            // the quadrature points are computed over the refence element [-1,1]*[-1,1],
            // as such the updated quadrature points and weights will always be the same for every element--> REDUNDANT WORK.
            // Yet, in the case we would like to change the degree of the quadrature formula over a
            // specific element, this function would make this process possible and less difficult to implement.
            // For those reasons, the _update_quad() method is called only in the case the _n member of the element
            // is different from the one of the previous element. Such check is implemented in the update_current() method
            if constexpr (DIM == 2)
            {
                //update the quadrature points over the reference element
                {
                    _quad_points.clear();

                    //loop over y coordinate
                    for(auto j : this->_comp_quad_c(1))
                    {
                        //loop over x coordinate
                        for(auto i : this->_comp_quad_c(0))
                        {
                            Point my_point(i,j);
                            _quad_points.emplace_back(my_point);
                        }
                    }

                }
                
                

            }
            else
            {
                // Update the quadrature points over the reference element

                {
                    _quad_points.clear();

                    //loop over x coordinate
                    for(auto i : this->_comp_quad_c(0))
                    {
                        Point my_point(i);
                        _quad_points.emplace_back(my_point);
                    }

                }
                
                

            }
            

            return;
        };
        
        // a method to update the D_ref matrix
        void _update_D()
        {
            // Clear _D_ref and resize
            _D_ref.resize(_nq[0]*((DIM ==2)?_nq[DIM-1]:1),_nq[0]*((DIM==2)?_nq[DIM-1]:1));
            // Compute the quadrature nodes on the quadrature points and insert them in D
            unsigned int index = 0;

            if constexpr(DIM == 2)
            {
                //loop over the y coordinate
                for(auto j : this->_comp_quad_w(1))
                {
                    //loop over x coordinate
                    for(auto i : this->_comp_quad_w(0))
                    {
                        _D_ref.coeffRef(index,index)=(i*j);
                        index++;
                    }
                }

            }
            else
            {
                //loop over the x coordinate
                for(auto i : this->comp_quad_w(0))
                {
                    _D_ref.coeffRef(index,index)=(i);
                    index++;
                }

            }
            

            return;
        }

        // a method to update the Bx_ref and By_ref matrix:
        // i.e. to compute the spectral Legendre Gauss Lobatto derivative matrix d at the np LGL nodes x (on [-1,1]),
        // please note: this function needs the correct values of the  quadrature point on the 1D reference interval to work correctly,
        // as such it is necessary to update the nodes of the Quadrature object before computing the matrix. This is, if the number of quadrature points
        // changed between the dimensions inside a certain element or between elements
        void _update_B( const unsigned int& dir)
        {
            std::vector<double> x = this->_comp_quad_c(dir);
            unsigned int np = x.size();
            
            // First of all, as a precautionary measure, we erase the content of the _der_matrix member
            _B_ref[dir].resize(this->getNQ()[0],this->getNQ()[0]);
            const unsigned int n = np -1;
            const unsigned int dof = DoFHandler(this->getDeg()).getDeg()[dir] + 1;
            std::vector<double> lnx(np);
            //compute the legendre polynomials over the LGL nodes
            FETools::Quadrature::legendre_pol(x, n, lnx);
            for(unsigned int j = 0; j < dof ; j++)
            {
                for(unsigned int i = 0; i < np; ++i)
                {
                    if(i != j)
                    {
                        _B_ref[dir].coeffRef(i,j) = lnx[i]/((x[i]-x[j])*lnx[j]);
                    }
                }
            }
            _B_ref[dir].coeffRef(0,0) = -0.25*n*np;
            _B_ref[dir].coeffRef(np-1,np-1) = 0.25*n*np;
        
            

            return;

        }
        
        /// a method to update the jacobian of the current element
        void _update_J()
        {
            if constexpr(DIM == 2)
            {
                _J_invT.resize(DIM, DIM);
                _J_invT = _current_elem.jacobian().inverse().transpose();

            }
            else
            {
                _J_invT.resize(DIM,DIM);
                _J_invT << 1/_current_elem.jacobian();
            }
                
            return;
        };

        // a method for simple class initialization
        std::array<unsigned int, DIM> classInit(const unsigned int& nx,const unsigned int& ny)
        {
            if constexpr (DIM == 2)
            {
                return {nx,ny};
            }
            else
            {
                return {nx};
            }
        }

        std::array<SparseMatrix<double>, DIM> classInit(const unsigned int& nqx, const unsigned int& nqy, const unsigned int& re )
        {
            // Instantiate a dof handler object to get the dof number
            unsigned int dof(DoFHandler(re).dof_per_cell());
            if constexpr (DIM == 2)
            {
                return {SparseMatrix<double>(nqx,dof),SparseMatrix<double>(nqy,dof)};
            }
            else
            {
                return {SparseMatrix<double>(nqx,dof)};
            }
        }


    private:

        
        // First, the element currently considered
        ElementType _current_elem;
        // To compute the quadrature coordinates and weights necessary
        std::array<Quadrature,DIM> _qr;
        // The degree of the FE space considered  --> number of quadrature nodes in each direction of the element
        // (for this implementation, _r is supposed to be the same for both directions)
        const unsigned int _r;
        // The number of quadratue points used for the current element (for this implementation, supposed to be the same for the two dimensions);
        std::array<unsigned int, DIM> _nq;
        //  A vector containig the coordinates of the quadrature points of the current elements
        std::vector<Point> _quad_points;
        // An Eigen Diagonal Matrix storing the values of the quadrature points for the single quadrature points stored inn _quad_points
        // ----> computed just once if the number of quadrature points are the same for all the elements. Otherwhise there is a need to recompute
        // it in the update_current() method...
        SparseMatrix<double> _D_ref;
        // Eigen Sparse Matrixes storing the values of the gradient of the basis functions over the x and y directions, given the quadrature points of the elements
        // ----> computed just one if the number of quadrature points is the same for all the elements. Otherwhise there is a need to recompute 
        // it in the update_current() method
        std::array<SparseMatrix<double>, DIM> _B_ref;
        // Jacobian of the element
        MatrixXd _J_invT;
        // Extensions of the previous local matrixes
        ExtentD<SparseMatrix<double>> _Dcell;
        ExtentJ<MatrixXd> _Jcell_invT;
        ExtentB<SparseMatrix<double>> _Bcell;
        // A counter to specify the element of the mesh currently examined
        unsigned int pun;
    };
}


#endif