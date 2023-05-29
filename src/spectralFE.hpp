

//===========================================================
//         HEADER FILE FOR THE FE SPECTRAL CLASS 
//===========================================================

#ifndef FE
#define FE

#include <numeric>
#include <algorithm>
#include <functional>
#include <tuple>
#include "elements.hpp"
#include "mesh.hpp"
#include "mat_utilities.hpp"
#include "quadrature.hpp"
#include "functions.hpp"

namespace FETools
{


    template <unsigned short DIM>
    class SpectralFE
    {
    public:


        // Constructor
        SpectralFE(const unsigned int &re = r)
            : _current_elem(),
              _qr(),
              _r(re),
              _quad_points(),
              _D_ref((re+1) * DIM == 2 ? (re+1) : 1 /* n_q */, (re+1) * DIM == 2 ? (re+1) : 1 /* n_q */),
              _B_ref(initMatArray(re)), // dof_per_cell = n_q
              _J_invT(DIM, DIM),
              _Dcell(_D_ref),
              _Jcell_invT(_J_invT, _current_elem),
              _Bcell(_B_ref, _current_elem, _r){};

        // Default getters
        const Element<DIM>& getCurrent() const
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
            return _current_elem.getNQ();
        }

        const std::vector<Point<DIM>>& getQPoints() const
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

        const ExtendMat::ExtendedD<SparseMatrix<double>,DIM>& D_cell() const
        {
            return _Dcell;
        }

        const ExtendMat::ExtendedJ<MatrixXd,DIM>& J_cell_invT() const
        {
            return _Jcell_invT;
        }

        const ExtendMat::ExtendedB<SparseMatrix<double>,DIM>& B_cell() const
        {
            return _Bcell;
        }

        double detJ() const
        {
            if constexpr(DIM == 2)
                return _current_elem.template jacobian<Matrix2d>().determinant();
            else
                return _current_elem.template jacobian<double>();
            
        }


        // a method to update the element currently considered
        void update_current(const Element<DIM> &geoele)
        {
        
            
            if(geoele.getNQ() == this->getNQ() && this->getQuad()[0].getW().size() != 0)
            {
                // In this first case, we just need to update the jacobian of the element after updating the 
                // element itself:
                _current_elem = geoele;
                // update the Jacobian
                this->_update_J();
                return;
            }
            else
            {
                // If the number of quadrature points along the two directions are different for the
                // new element with respect to the previous one, we have to update the quadrature points, the D and the
                // B matrixes
                unsigned int nqx = _current_elem.getNQ()[0];
                unsigned int nqy;
                if constexpr (DIM == 2)
                {
                    nqy = _current_elem.getNQ()[1];
                }
                // copute empty flag
                bool emptyx = this->getQuad()[0].getW().size() == 0;
                bool emptyy = this->getQuad()[1].getW().size() == 0;
                // update the current element  
                _current_elem = geoele;
                // update the Jacobian
                this->_update_J();
                // update the quadrature points
                this->_update_quad();
                // update the evaluation of quadrature points
                this->_update_D();
                // In the 1D case just recompute the Bx matrix
                if(nqx != this->getNQ()[0] || emptyx)
                {
                    this->_update_B(0);
                }
                if constexpr (DIM ==2)
                {
                    
                    // In this case, recompute By
                    if(nqy != this->getNQ()[1] || emptyy)
                    {
                        this->_update_B(1);
                    }
                    

                }

            }
            
            return;
        };
        
        // member function to map the quadrature points on to the specific element
        Point<DIM> quadrature_point(const unsigned int &q,const DoFHandler<DIM>& dof) const
        {

            //  FIRST APPROACH:
            // compute the quadrature point usinge the direct map function with the quadrature pointson the reference element

            // if constexpr (DIM == 2)
            // {
            //     double x, y;
            //     std::tie(x,y) = _current_elem.directMap<double>((this->getQPoints()[q]).getX(),(this->getQPoints()[q]).getY());
            //     Point<DIM> my_point(x,y);
            //     return my_point;
            // }
            // else
            // {
            //     double x;
            //     x = _current_elem.directMap<Matrix2d>((this->getQPoints()[q]).getX());
            //     Point<DIM> my_point(x);
            //     return my_point;
            // }

            //  SECOND APPROACH
            // Get the quarature point for the global mesh computed by a DoFHandler object ALREADY INITIALISED
            
            // Using the index of the element currently considered, get the coordinates of the quadrature points
            
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
                            Point<DIM> my_point(i,j, 0);
                            
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
                        Point<DIM> my_point(i, 0, 0);
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
            _D_ref.resize(this->getNQ()[0]*((DIM ==2)?this->getNQ()[DIM-1]:1),this->getNQ()[0]*((DIM==2)?this->getNQ()[DIM-1]:1));
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
                for(auto i : this->_comp_quad_w(0))
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
            const unsigned int dof = DoFHandler<DIM>(this->getDeg(), this->getDeg()).getDeg()[dir] + 1;
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
                _J_invT = _current_elem.template jacobian<Matrix2d>().inverse().transpose();

            }
            else
            {
                _J_invT.resize(DIM,DIM);
                _J_invT << 1/_current_elem.template jacobian<double>();
            }
                
            return;
        };


        // a static method to allow for correct initializiaion of array of Eigen::SparseMatrix
        // given the dimension of the problem (DIM)
        static std::array<SparseMatrix<double>,DIM> initMatArray(const unsigned int &re)
        {
            if constexpr (DIM == 1)
                return {SparseMatrix<double>((re+1), (re+1))};
            else
                return {SparseMatrix<double>((re+1), (re+1)*(re+1)), SparseMatrix<double>((re+1), (re+1)*(re+1))};
        }


    private:

        
        // First, the element currently considered (taken as reference as to avoid making a copy)
        Element<DIM> _current_elem;
        // To compute the quadrature coordinates and weights necessary
        std::array<Quadrature,DIM> _qr;
        // The degree of the FE space considered  --> number of quadrature nodes in each direction of the element
        // (Here we assume the deg of the FE spaces is the same for all directions, to ease off the implementation of the code)
        const unsigned int _r;
        //  A vector containig the coordinates of the quadrature points of the current elements
        std::vector<Point<DIM>> _quad_points;
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
        ExtendMat::ExtendedD<SparseMatrix<double>,DIM> _Dcell;
        ExtendMat::ExtendedJ<MatrixXd, DIM> _Jcell_invT;
        ExtendMat::ExtendedB<SparseMatrix<double>, DIM> _Bcell;
    };
}


#endif