

//  ==========================================================================
// CUSTOM CLASSSES TO IMPLEMENT MATRIX FREE MATRIX PRODUCT ON THE CELL MATRIXES
//  ==========================================================================


#ifndef UTIL
#define UTIL

#include <iostream>
#include <array>
#include <vector>
#include <string>
#include "mesh.hpp"
#include <Eigen/Dense>
#include <Eigen/SparseCore>

// since Bref_cell,J_cell and D_cell are made up of copies of the same values repeated in certain positions,
// this lets us compute the matrix products in a more efficent way simply by accessing the Jacobian of the element
// D_ref, and Bx_ref/ By_ref.
// The following classes are views(proxies) to implement matrix multiplications between 
// the local matrixes without explicitily defining them.
        
namespace ExtendMat
{

    //  Class for Bref_cell extention form Bx_ref and By_ref
    template<typename Derived, unsigned int DIM>
    class ExtendedB: public Eigen::EigenBase<Derived>
    {
    public:
        ExtendedB() = default;
        ExtendedB(const std::array<SparseMatrix<double>,DIM>& B, const Element<DIM>& geoele, const unsigned int& re = r)
        : _B(B),
        _current(geoele),
        _r(re)
        {}

        int rows() const
        {
            return DIM * _current.getNQ()[0] * (DIM==2?_current.getNQ()[DIM-1]:1)/*DIM * n_q*/;
        }

        int cols() const
        {
            // REMARK:
            // dofs per cell and number of quadrature points are supposeed to be the same, still
            // this implementation better explains how the B matrix is built 
            // Ideally, we could just get the number of the cols form nq, since the
            // number of quadarture points and dofs is supposed to be the same
            return DoFHandler<DIM>(_r, _r).dof_per_cell() /* dof_per_cell */;
        }

        // Return sub-matrix in each iteration
        double operator()(const unsigned int& i/* nqx*nqy*DIM*/, const unsigned int& j/*nqx*nqy*/) const
        {
            // REMARK
            // Here again, we could have just used a single variable, yet for better clarity, the choice was
            // to "keep the dimensions separated"
            const unsigned int dofx = DoFHandler<DIM>(_r, _r).getDeg()[0] +1;
            const unsigned int dofy = DoFHandler<DIM>(_r, _r).getDeg()[1] +1 ;

            if constexpr(DIM == 2)
            {
                // lines of Bx_ref
                if(i%2 == 0 /*even rows (x dedicated)*/)
                {
                    if((i>>1)/_current.getNQ()[0] == j/dofx) // i and j indexing on an element on the "diagonal" block
                        return _B[0].coeff((i/2)%_current.getNQ()[0],j%dofx);
                    else 
                        return 0.;
                }
                // lines of By_ref
                else /*odd rows*/
                {
                    if((i>>1)%_current.getNQ()[1] == j%dofy) // i and j indexing on an elemente on the diagonal block
                        return _B[1].coeff((i>>1)/_current.getNQ()[1], j/dofy);
                    else 
                        return 0.;
                }

            }
            else
            {
                if(i/_current.getNQ()[0] == j/dofx) // i and j indexing on an element on the "diagonal" block
                        return _B[0].coeff(i%_current.getNQ()[0],j%dofx);
                    else 
                        return 0.;
            } 
                
        }

        double transpose(const unsigned int& i/* nqx*nqy*DIM */, const unsigned int& j/* nqx*nqy */) const
        {
            return this->operator()(j,i);
                
        }

        

    private:
        const std::array<SparseMatrix<double>,DIM>& _B;
        const Element<DIM>& _current;
        const unsigned int& _r;
    };

    //  Class for D_cell extension form D_ref matrix
    
    template<typename Derived, unsigned int DIM>
    class ExtendedD: public Eigen::EigenBase<Derived>
    {
    public:
        ExtendedD() = default;
        ExtendedD(const SparseMatrix<double> &A)
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
        double operator()(const unsigned int& i, const unsigned int& j) const
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
    template<typename Derived, unsigned int DIM>
    class ExtendedJ: public Eigen::EigenBase<Derived>
    {
    public:
        ExtendedJ() = default;
        ExtendedJ(const MatrixXd &A, const Element<DIM>& geoele)
        : _A(A),
        _current(geoele)
        {}

        int rows() const
        {
            return _A.rows()*(DIM==1?_current.getNQ()[0]:_current.getNQ()[0]*_current.getNQ()[DIM-1]);
        }

        int cols() const
        {
            return _A.cols()*(DIM==1?_current.getNQ()[0]:_current.getNQ()[0]*_current.getNQ()[DIM-1]);
        }

        // Return sub-matrix in each iteration
        double operator()(const unsigned int& i, const unsigned int& j) const
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
        MatrixXd operator*(const ExtendedB<SparseMatrix<double>,DIM> &B) const
        {

            SparseMatrix<double> result(this->rows(), B.cols());
            for (int i = 0; i <this->rows(); ++i)
            {
                for (int j = 0; j < B.cols(); ++j)
                {
                    for (int k = 0; k < this->cols(); ++k)
                    {
                        result.coeffRef(i, j) += this->operator()(i, k) * B(k, j);
                    }
                        
                }


            }
            return result;
        }


        double transpose(const unsigned int& i, const unsigned int& j) const
        {
            return this->operator()(j,i);
        }

        

    private:
        const MatrixXd& _A;
        const Element<DIM>& _current;
    };

}




//  ==========================================================================






#endif





