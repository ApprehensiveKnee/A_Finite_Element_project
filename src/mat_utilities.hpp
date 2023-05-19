

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
// The following classes are views (PROXIES) to implement matrix multiplications between 
// the local matrixes without explicitily defining them.

//PLEASE NOTE: the defined custom classes are not compatible for 
// operations such as matrix product with Eigen classes. Still I have defined some methods
// to allow for the creation of eigen matrixes which are the same as the classes defined,
// may their usage be required

// Definition of concept later used for GeneralProduct template restriction
// It should accept all objects with the following requirements (ideally it would also include eigen matrixes)
template<typename M>
concept IsGeneralMatrix = requires(M mat)
{
    //{ typename M::Scalar{} } -> std::same_as<double>;
    { mat.rows() } -> std::same_as<int>;
    { mat.cols() } -> std::same_as<int>;
    { mat(0,0) } -> std::same_as<double>;
};

namespace ExtendMat
{

    //==========================================================================================================
    //  Class for Bref_cell extension form Bx_ref and By_ref
    template<typename Derived, unsigned short DIM>
    class ExtendedB: public Eigen::EigenBase<Derived>
    {
    public:
        ExtendedB() = default;
        ExtendedB(const std::array<Derived,DIM>& B, const Element<DIM>& geoele, const unsigned int& re = r, const bool& transpose = false)
        : _B(B),
        _current(geoele),
        _r(re),
        _transpose(transpose)
        {}

        int rows() const
        {
            if(_transpose)
                return DoFHandler<DIM>(_r, _r).dof_per_cell(); /* dof_per_cell */
            else
                return DIM * _current.getNQ()[0] * (DIM==2?_current.getNQ()[DIM-1]:1);/*DIM * n_q*/
        }

        int cols() const
        {
            // REMARK:
            // dofs per cell and number of quadrature points are supposeed to be the same, still
            // this implementation better explains how the B matrix is built 
            // Ideally, we could just get the number of the cols form nq, since the
            // number of quadarture points and dofs is supposed to be the same
            if(_transpose)
                return DIM * _current.getNQ()[0] * (DIM==2?_current.getNQ()[DIM-1]:1);/*DIM * n_q*/
            else
                return DoFHandler<DIM>(_r, _r).dof_per_cell(); /* dof_per_cell */
        }

        // Access operator to implement matrix free computation: no need to explicitly generate the matrix
        double operator()(const unsigned int& i/* nqx*nqy*DIM - no transpose */, const unsigned int& j/*nqx*nqy - no transpose*/) const
        {
            
            // REMARK
            // Here again, we could have just used a single variable, yet for better clarity, the choice was
            // to "keep the dimensions separated"
            if(_transpose)
            {
                const unsigned int dofx = DoFHandler<DIM>(_r, _r).getDeg()[0] +1;
                const unsigned int dofy = DoFHandler<DIM>(_r, _r).getDeg()[1] +1;

                if constexpr(DIM == 2)
                {
                    // B
                    // lines of Bx_ref
                    if(j%2 == 0 /*even cols (x dedicated)*/)
                    {
                        if((j>>1)/_current.getNQ()[0] == i/dofx) // i and j indexing on an element on the "diagonal" block
                            return _B[0].coeff((j/2)%_current.getNQ()[0],i%dofx);
                        else 
                            return 0.;
                    }
                    // lines of By_ref
                    else /*odd colss*/
                    {
                        if((j>>1)%_current.getNQ()[1] == i%dofy) // i and j indexing on an elemente on the diagonal block
                            return _B[1].coeff((j>>1)/_current.getNQ()[1], i/dofy);
                        else 
                            return 0.;
                    }

                }
                else
                {
                    if(j/_current.getNQ()[0] == i/dofx) // i and j indexing on an element on the "diagonal" block
                            return _B[0].coeff(j%_current.getNQ()[0],i%dofx);
                        else 
                            return 0.;
                }

            }
            else
            {

                const unsigned int dofx = DoFHandler<DIM>(_r, _r).getDeg()[0] +1;
                const unsigned int dofy = DoFHandler<DIM>(_r, _r).getDeg()[1] +1;

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
                
                
        }

        ExtendedB<Derived,DIM> transpose() const
        {
            return ExtendedB<Derived,DIM>(_B, _current, _r, !_transpose);
        }

        //_______________________________ EIGEN SUPPORT ______________________________
        // Some methods to return the Eigen corresponging matrixes (not really going to be used
        // they are just for completeness)
        MatrixXd mat() const
        {
            // Define the Eigen matrix to be returned
            MatrixXd matrix(this->rows(), this->cols());
            for (unsigned int i = 0; i < this->rows(); i++)
            {
                for (unsigned int j = 0; j < this->cols(); j++)
                {
                    matrix(i,j) = this->operator()(i,j);
                }
                
            }
            return matrix;
            
        }

        MatrixXd matT() const
        {
            MatrixXd matrix(this->cols(), this->rows());
            for (unsigned int i = 0; i < this->cols(); i++)
            {
                for (unsigned int j = 0; j < this->rows(); j++)
                {
                    matrix(i,j) = this->operator()(j,i);
                }
                
            }
            return matrix;
            
        }
        //_______________________________ EIGEN SUPPORT ______________________________

        // Destructor

        ~ExtendedB() = default;

        

    private:
        const std::array<Derived,DIM>& _B;
        const Element<DIM>& _current;
        const unsigned int& _r;
        const bool _transpose;
    };

    //==========================================================================================================
    //  Class for D_cell extension form D_ref matrix
    // The traspose of this matrix would be the matrix itself since it is diagonal

    template<typename Derived, unsigned short DIM>
    class ExtendedD: public Eigen::EigenBase<Derived>
    {
    public:
        ExtendedD() = default;
        ExtendedD(const Derived &A)
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

        // Access operator to allow for matrix free computation: no need to explicitly generate the matrix
        double operator()(const unsigned int& i, const unsigned int& j) const
        {
            if( i == j)
                return _A.coeff(i/DIM, j/DIM);
            else
                return 0.;
        }


        //_______________________________ EIGEN SUPPORT ______________________________
        // Some methods to return the Eigen corresponging matrixes (not really going to be used
        // they are just for completeness)

        SparseMatrix<double> mat() const
        {
            // Define the Eigen matrix to be returned
            SparseMatrix<double> matrix(this->rows(), this->cols());
            for (unsigned int i = 0; i < this->rows(); i++)
            {
                for (unsigned int j = 0; j < this->cols(); j++)
                {
                    matrix.coeffRef(i,j) = this->operator()(i,j);
                }
                
            }
            return matrix;
            
        }

        SparseMatrix<double> matT() const
        {
            SparseMatrix<double> matrix(this->cols(), this->rows());
            for (unsigned int i = 0; i < this->cols(); i++)
            {
                for (unsigned int j = 0; j < this->rows(); j++)
                {
                    matrix.insert(i,j) = this->operator()(j,i);
                }
                
            }
            return matrix;
            
        }

        //_______________________________ EIGEN SUPPORT ______________________________

        // Destructor

        ~ExtendedD() = default;

    private:
        const Derived& _A;
    };





    //==========================================================================================================

    //  Class for J_cell extention from Jacobian of the element
    template<typename Derived, unsigned short DIM>
    class ExtendedJ: public Eigen::EigenBase<Derived>
    {
    public:
        ExtendedJ() = default;
        ExtendedJ(const Derived &A, const Element<DIM>& geoele, const bool& transpose = false)
        : _A(A),
        _current(geoele),
        _transpose(transpose)
        {}

        int rows() const
        {
            if(_transpose)
                return _A.cols()*(DIM==1?_current.getNQ()[0]:_current.getNQ()[0]*_current.getNQ()[DIM-1]);
            else
                return _A.rows()*(DIM==1?_current.getNQ()[0]:_current.getNQ()[0]*_current.getNQ()[DIM-1]);
        }

        int cols() const
        {
            if(_transpose)
                return _A.rows()*(DIM==1?_current.getNQ()[0]:_current.getNQ()[0]*_current.getNQ()[DIM-1]);
            else
                return _A.cols()*(DIM==1?_current.getNQ()[0]:_current.getNQ()[0]*_current.getNQ()[DIM-1]);
        }

        // Access operator to allow for matrix free computation
        double operator()(const unsigned int& i, const unsigned int& j) const
        {

            if(_transpose)
            {

                if constexpr(DIM == 2)
                {

                    if(j%2 == 0 /* col is even*/)
                    {
                        if(i%2 == 0 /*row is even*/ && j == i /*on the diagonal*/ )
                        {
                            return _A.coeff(0, 0);
                        }
                        else if(i%2 == 1 /*row is */ && i == j+1 /*on the lower diagonal*/)
                        {
                            return _A.coeff(0,1);
                        }
                        else
                            return 0.;
                    }
                    else /*col is odd*/
                    {
                        if(i%2 == 0 /*row is even*/ && i == j-1 /*on the upper diagonal*/)
                        {
                            return _A.coeff(1, 0);
                        }
                        else if(i%2 == 1/*row is odd*/ && j == i/*on the diaglona*/)
                        {
                            return _A.coeff(1,1);
                        }
                        else
                            return 0.;
                    }

                }
                else
                {
                    if(i == j /*on the diagonal*/)
                        return _A.coeff(0,0);
                    else
                        return 0.;
                }

            }
            else
            {
                if constexpr(DIM == 2)
                {

                    if(i%2 == 0 /* row is even*/)
                    {
                        if(j%2 == 0 /*col is even*/ && i == j /*on the diagonal*/ )
                        {
                            return _A.coeff(0, 0);
                        }
                        else if(j%2 == 1 /*col is */ && j == i+1 /*on the upper diagonal*/)
                        {
                            return _A.coeff(0,1);
                        }
                        else
                            return 0.;
                    }
                    else /*row is odd*/
                    {
                        if(j%2 == 0 /*col is even*/ && j == i-1 /*on the lower diagonal*/)
                        {
                            return _A.coeff(1, 0);
                        }
                        else if(j%2 == 1/*col is odd*/ && j == i/*on the diaglona*/)
                        {
                            return _A.coeff(1,1);
                        }
                        else
                            return 0.;
                    }

                }
                else
                {
                    if(i == j /*on the diagonal*/)
                        return _A.coeff(0,0);
                    else
                        return 0.;
                }

            }
                
        }

        

        //_______________________________ EIGEN SUPPORT ______________________________
        // Some methods to return the Eigen corresponging matrixes (not really going to be used
        // they are just for completeness)
        MatrixXd mat() const
        {
            // Define the Eigen matrix to be returned
            MatrixXd matrix(this->rows(), this->cols());
            for (unsigned int i = 0; i < this->rows(); i++)
            {
                for (unsigned int j = 0; j < this->cols(); j++)
                {
                    matrix(i,j) = this->operator()(i,j);
                }
                
            }
            return matrix;
            
        }

        MatrixXd matT() const
        {
            MatrixXd matrix(this->cols(),this->rows());
            for (unsigned int i = 0; i < this->cols(); i++)
            {
                for (unsigned int j = 0; j < this->rows(); j++)
                {
                    matrix(i,j) = this->operator()(j,i);
                }
                
            }
            return matrix;
            
        }
        //_______________________________ EIGEN SUPPORT ______________________________

        
        // Destructor
        ~ExtendedJ() = default;

    private:
        const Derived& _A;
        const Element<DIM>& _current;
        const bool _transpose;
    };

    //==========================================================================================================

    
    // Finally define a class that will serve the general purpose of representing the 
    // products between the previously defined classes:
    // To better define this class, use c++ concepts to specify some necessary 
    // proprieties of the matrix involved in matrix product

    // PLEASE NOTE: I thought that implementing such proxy to allow for matrix free computation would speed up the code
    //  as all the accesses would have been cached, but it turns out that simply managing it all with
    // Eigen matrix is acutually a far better solution in terms of performace.


    template <class T, class U>
    requires ((IsGeneralMatrix<T> || IsDouble<T>) && IsGeneralMatrix<U>)
    class GeneralProduct
    {
    public:
        // Constructor
        GeneralProduct(const T& A,const U& B, const bool& transpose = false):
            _A(A),
            _B(B),
            _transpose(transpose)
            {if constexpr(IsGeneralMatrix<T>)
                    // check form correspondence of number of rows and cols
                    if(A.cols() != B.rows())
                    {
                        throw std::runtime_error("The number of cols of first operand must be equal to the number of cols of the second operand\n");
                    }
                        
            };

        // Definition of rows and cols methods
        int rows() const
        {
            if(_transpose)
            {
                return _B.cols();
            }
            else
            {
                if constexpr (IsDouble<T>)
                    return _B.rows();
                else
                    return _A.rows();
            }
            
        }

        int cols() const
        {
            if(_transpose)
            {
                if constexpr (IsDouble<T>)
                {
                    return _B.cols();
                }
                else
                    return _A.rows();
            }
            else
            {
                
                return _B.cols();
            }
            
        }

        // Access operator
        double operator()(const unsigned int& i, const unsigned int& j) const
        {
            if(_transpose)
            {
                // Here define the different cases
                if constexpr (IsDouble<T>) // first operand is a double
                {
                    return _A * _B(j,i);
                }
                else // both operands are matrixes
                {
                    double sum(0);
                    for (unsigned int k = 0; k < _A.cols(); k++)
                    {
                        sum += _A(j,k)*_B(k,i);
                        // sum += _A.transpose(i,k)*_B.transpose(k,j);
                    }
                    
                    return sum;
                }
            }
            else
            {
                // Here define the different cases
                if constexpr (IsDouble<T>) // first operand is a double
                {
                    return _A * _B(i,j);
                }
                else // both operands are matrixes
                {
                    double sum(0);
                    for (unsigned int k = 0; k < _A.cols(); k++)
                    {
                        sum += _A(i,k)*_B(k,j);
                    }
                    return sum;
                }

            }
            
                
        }


        ExtendMat::GeneralProduct<T,U> transpose() const
        {
            return ExtendMat::GeneralProduct<T,U>(_A,_B, !_transpose);
        }

        // Getters (for debugging)

        const T& getA() const 
        {
            return _A;
        }

        const U& getB() const
        {
            return _B;
        }

        // Destructor
        ~GeneralProduct() = default;

    private:
        const T _A; // const reference to the first operand
        const U _B; // const reference to the second operand
        const bool _transpose;
    };


    //  ==========================================================================

    //  ========================    OPERATOR * OVERLOADING  ========================

    template<typename Derived>
    concept EigenMatrix = std::is_base_of_v<Eigen::EigenBase<Derived>, Derived>;

    // ========================= OVERLOADING OF OPERATOR * USING EIGEN CLASSES ===============================

    template<EigenMatrix Derived, unsigned short DIM>
    const SparseMatrix<double> operator*(const ExtendMat::ExtendedD<Derived, DIM>& D, const double& coeff)
    {
        SparseMatrix<double> result(D.rows(), D.cols());

        for (int i = 0; i < D.rows(); ++i)
                result.coeffRef(i, i) = coeff * D(i, i);

        return result;
    }

    template<EigenMatrix Derived1, EigenMatrix Derived2, unsigned short DIM>
    const MatrixXd operator*(const ExtendMat::ExtendedJ<Derived1,DIM> &J, const ExtendMat::ExtendedB<Derived2,DIM> &B)
    {
        MatrixXd result = MatrixXd::Zero(J.rows(), B.cols());

        for (int i = 0; i <J.rows(); ++i)
        {
            for (int j = 0; j < B.cols(); ++j)
            {
                for (int k = 0; k < J.cols(); ++k)
                {
                    result(i,j) += J(i, k) * B(k, j);
                }
                    
            }


        }
        return result;
    }


    



   //  ========================= OVERLOADING OF OPERATOR * USING EIGEN CLASSES ===============================

   //  ========================= OVERLOADING OF OPERATOR * WITH CUSTOM CLASSES ===============================

    // Overloading of operator * to allow for product between "virtual matrix classes" ExtendedB and ExtendedJ
    /*

    template<EigenMatrix Derived1, EigenMatrix Derived2, unsigned short DIM>
    GeneralProduct<ExtendedB<Derived1,DIM>, ExtendedJ<Derived2, DIM>> operator *( const ExtendedB<Derived1,DIM>& B,const ExtendedJ<Derived2, DIM>& J)
    {
        return GeneralProduct<ExtendedB<Derived1,DIM>, ExtendedJ<Derived2, DIM>>(B, J);
    }

    
    template<EigenMatrix Derived1,EigenMatrix Derived2, unsigned short DIM>
    GeneralProduct<ExtendedJ<Derived1, DIM>, ExtendedB<Derived2,DIM>> operator *(const ExtendedJ<Derived1, DIM>& J,const ExtendedB<Derived2,DIM>& B)
    {
        return GeneralProduct<ExtendedJ<Derived1,DIM>, ExtendedB<Derived2, DIM>>(J, B);
    }




    // Operator * overloading to support multiplication with custom expression class (coeff * D_ref)


    template<EigenMatrix Derived, unsigned short DIM>
    GeneralProduct<double, ExtendedD<Derived,DIM>> operator*(const ExtendedD<Derived,DIM>& D, const double& coeff)
    {
        return GeneralProduct<double, ExtendedD<Derived,DIM>>(coeff,D);
    }

    // Operator * overloading to support multiplication with custom expression class (D_ref*coeff)
    template<EigenMatrix Derived, unsigned short DIM>
    GeneralProduct<double, ExtendedD<Derived,DIM>> operator*(const double& coeff, const ExtendedD<Derived,DIM>& D)
    {
        return GeneralProduct<double, ExtendedD<Derived,DIM>>(coeff, D);
    }

    */




    // overloading of operator * to allow for product between general product objects
    template <class T, class U,class X, class Y>
    ExtendMat::GeneralProduct<ExtendMat::GeneralProduct<T,U>, ExtendMat::GeneralProduct<X,Y>> operator*(const ExtendMat::GeneralProduct<T,U>& operand1,const ExtendMat::GeneralProduct<X,Y>& operand2 )
    {
        return ExtendMat::GeneralProduct<ExtendMat::GeneralProduct<T,U>, ExtendMat::GeneralProduct<X,Y>>(operand1,operand2, false);
    }

    // overloading of operator * to allow for product between general product objects and doubles
    template <class T, class U>
    ExtendMat::GeneralProduct< double, ExtendMat::GeneralProduct<T,U>> operator*(const double& operand1,const ExtendMat::GeneralProduct<T,U>& operand2 )
    {
        return ExtendMat::GeneralProduct<double, ExtendMat::GeneralProduct<T,U>>(operand1,operand2, false);
    }

    template <class T, class U>
    ExtendMat::GeneralProduct< double, ExtendMat::GeneralProduct<T,U>> operator*(const ExtendMat::GeneralProduct<T,U>& operand1, const double& operand2 )
    {
        return ExtendMat::GeneralProduct<double, ExtendMat::GeneralProduct<T,U>>(operand2,operand1, false);
    }

    //  ========================= OVERLOADING OF OPERATOR * WITH CUSTOM CLASSES ===============================

}


//  ========================    OPERATOR * OVERLOADING  ========================
#endif





