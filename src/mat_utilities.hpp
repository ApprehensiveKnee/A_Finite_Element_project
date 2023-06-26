

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

// Definition of a concept later used for GeneralProduct template restriction
// It should accept all objects with the following requirements (ideally it would also include eigen matrixes)
template <typename M>
concept IsGeneralMatrix = requires(M mat) {
    //{ typename M::Scalar{} } -> std::same_as<double>;
    {
        mat.rows()
    } -> std::same_as<int>;
    {
        mat.cols()
    } -> std::same_as<int>;
    {
        mat(0, 0)
    } -> std::same_as<double>;
};

namespace ExtendMat
{

    //==========================================================================================================
    //  Class for Bref_cell extension form Bx_ref and By_ref
    template <typename Derived, unsigned short DIM>
    class ExtendedB : public Eigen::EigenBase<Derived>
    {
    public:
        ExtendedB() = default;
        ExtendedB(const std::array<Derived, DIM> &B, const Element<DIM> &geoele, const bool &transpose = false)
            : _B(B),
              _current(geoele),
              _transpose(transpose)
        {
        }

        int rows() const
        {
            if (_transpose)
                return _current.getNQ()[0] * (DIM == 2 ? _current.getNQ()[DIM - 1] : 1); /* dof_per_cell */
            else
                return DIM * _current.getNQ()[0] * (DIM == 2 ? _current.getNQ()[DIM - 1] : 1); /*DIM * n_q*/
        }

        int cols() const
        {
            // REMARK: just get the number of the cols form nq, since the
            // number of quadarture points and dofs is supposed to be the same
            if (_transpose)
                return DIM * _current.getNQ()[0] * (DIM == 2 ? _current.getNQ()[DIM - 1] : 1); /*DIM * n_q*/
            else
                return _current.getNQ()[0] * (DIM == 2 ? _current.getNQ()[DIM - 1] : 1); /* dof_per_cell */
        }

        // Access operator to implement matrix free computation: no need to explicitly generate the matrix
        double operator()(const unsigned int &i /* nqx*nqy*DIM - no transpose */, const unsigned int &j /*nqx*nqy - no transpose*/) const
        {

            if (_transpose)
            {
                const unsigned int dofx = _current.getNQ()[0];
                const unsigned int dofy = _current.getNQ()[1];

                if constexpr (DIM == 2)
                {
                    // B
                    // lines of Bx_ref
                    if (j % 2 == 0 /*even cols (x dedicated)*/)
                    {
                        if ((j >> 1) / _current.getNQ()[0] == i / dofx) // i and j indexing on an element on the "diagonal" block
                            return _B[0].coeff((j / 2) % _current.getNQ()[0], i % dofx);
                        else
                            return zero_value;
                    }
                    // lines of By_ref
                    else /*odd cols*/
                    {
                        if ((j >> 1) % _current.getNQ()[1] == i % dofy) // i and j indexing on an elemente on the diagonal block
                            return _B[1].coeff((j >> 1) / _current.getNQ()[1], i / dofy);
                        else
                            return zero_value;
                    }
                }
                else
                {
                    if (j / _current.getNQ()[0] == i / dofx) // i and j indexing on an element on the "diagonal" block
                        return _B[0].coeff(j % _current.getNQ()[0], i % dofx);
                    else
                        return zero_value;
                }
            }
            else
            {

                const unsigned int dofx = _current.getNQ()[0];
                const unsigned int dofy = _current.getNQ()[1];

                if constexpr (DIM == 2)
                {
                    // lines of Bx_ref
                    if (i % 2 == 0 /*even rows (x dedicated)*/)
                    {
                        if ((i >> 1) / _current.getNQ()[0] == j / dofx) // i and j indexing on an element on the "diagonal" block
                            return _B[0].coeff((i / 2) % _current.getNQ()[0], j % dofx);
                        else
                            return zero_value;
                    }
                    // lines of By_ref
                    else /*odd rows*/
                    {
                        if ((i >> 1) % _current.getNQ()[1] == j % dofy) // i and j indexing on an elemente on the diagonal block
                            return _B[1].coeff((i >> 1) / _current.getNQ()[1], j / dofy);
                        else
                            return zero_value;
                    }
                }
                else
                {
                    if (i / _current.getNQ()[0] == j / dofx) // i and j indexing on an element on the "diagonal" block
                        return _B[0].coeff(i % _current.getNQ()[0], j % dofx);
                    else
                        return zero_value;
                }
            }
        }

        ExtendedB<Derived, DIM> transpose() const
        {
            return ExtendedB<Derived, DIM>(_B, _current, !_transpose);
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
                    matrix(i, j) = this->operator()(i, j);
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
                    matrix(i, j) = this->operator()(j, i);
                }
            }
            return matrix;
        }
        //_______________________________ EIGEN SUPPORT ______________________________

        // Updater for the reference matrix
        void update(const Element<DIM> &elem, const std::vector<double> &x, const unsigned int &r, const unsigned short dir)
        {
            _current = elem;
            unsigned int np = x.size();
            // First of all, as a precautionary measure, we erase the content of the _der_matrix member
            _B[dir].resize(elem.getNQ()[0], elem.getNQ()[0]);
            const unsigned int n = np - 1;
            const unsigned int dof = DoFHandler<DIM>(r, r).getDeg()[dir] + 1;
            std::vector<double> lnx(np);
            // compute the legendre polynomials over the LGL nodes
            FETools::Quadrature::legendre_pol(x, n, lnx);
            for (unsigned int j = 0; j < dof; j++)
            {
                for (unsigned int i = 0; i < np; ++i)
                {
                    if (i != j)
                    {
                        _B[dir].coeffRef(i, j) = lnx[i] / ((x[i] - x[j]) * lnx[j]);
                    }
                }
            }
            _B[dir].coeffRef(0, 0) = -0.25 * n * np;
            _B[dir].coeffRef(np - 1, np - 1) = 0.25 * n * np;
        }

        const std::array<Derived, DIM> &getRef() const
        {
            return _B;
        }

        // Destructor

        ~ExtendedB() = default;

    private:
        // Eigen Sparse Matrixes storing the values of the gradient of the basis functions over the x and y directions, given the quadrature points of the elements
        // ----> computed just one if the number of quadrature points is the same for all the elements. Otherwhise there is a need to recompute
        // it in the update_current() method
        std::array<Derived, DIM> _B;
        Element<DIM> _current;
        const bool _transpose;
        static constexpr double zero_value = 0.0;
    };

    //==========================================================================================================
    //  Class for D_cell extension form D_ref matrix
    // The traspose of this matrix would be the matrix itself since it is diagonal

    template <typename Derived, unsigned short DIM>
    class ExtendedD : public Eigen::EigenBase<Derived>
    {
    public:
        ExtendedD() = default;
        ExtendedD(const Derived &A)
            : _A(A)
        {
        }

        int rows() const
        {
            return _A.rows() * DIM;
        }

        int cols() const
        {
            return _A.cols() * DIM;
        }

        // Access operator to allow for matrix free computation: no need to explicitly generate the matrix
        double operator()(const unsigned int &i, const unsigned int &j) const
        {
            if (i == j)
                return _A.coeff(i / DIM, j / DIM);
            else
                return zero_value;
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
                    matrix.coeffRef(i, j) = this->operator()(i, j);
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
                    matrix.insert(i, j) = this->operator()(j, i);
                }
            }
            return matrix;
        }

        //_______________________________ EIGEN SUPPORT ______________________________

        // Updater for the reference matrix

        void update(const Element<DIM> &elem, std::array<std::vector<double>, DIM> &weights)
        {
            // Clear _A and resize
            _A.resize(elem.getNQ()[0] * ((DIM == 2) ? elem.getNQ()[DIM - 1] : 1), elem.getNQ()[0] * ((DIM == 2) ? elem.getNQ()[DIM - 1] : 1));
            // Compute the quadrature nodes on the quadrature points and insert them in D
            unsigned int index = 0;

            if constexpr (DIM == 2)
            {
                // loop over the y coordinate
                for (auto j : weights[1])
                {
                    // loop over x coordinate
                    for (auto i : weights[0])
                    {
                        _A.coeffRef(index, index) = (i * j);
                        index++;
                    }
                }
            }
            else
            {
                // loop over the x coordinate
                for (auto i : weights[0])
                {
                    _A.coeffRef(index, index) = (i);
                    index++;
                }
            }
        }

        const Derived &getRef() const
        {
            return _A;
        }

        // Destructor

        ~ExtendedD() = default;

    private:
        // An Eigen Diagonal Matrix storing the values of the quadrature points for the single quadrature points stored inn _quad_points
        // ----> computed just once if the number of quadrature points are the same for all the elements. Otherwhise there is a need to recompute
        // it in the update_current() method...
        Derived _A;
        static constexpr double zero_value = 0.0;
    };

    //==========================================================================================================

    //  Class for J_cell extention from Jacobian of the element
    template <typename Derived, unsigned short DIM>
    class ExtendedJ : public Eigen::EigenBase<Derived>
    {
    public:
        ExtendedJ() = default;
        ExtendedJ(const Derived &A, const Element<DIM> &geoele, const bool &transpose = false)
            : _A(A),
              _current(geoele),
              _transpose(transpose)
        {
        }

        int rows() const
        {
            if (_transpose)
                return _A.cols() * (DIM == 1 ? _current.getNQ()[0] : _current.getNQ()[0] * _current.getNQ()[DIM - 1]);
            else
                return _A.rows() * (DIM == 1 ? _current.getNQ()[0] : _current.getNQ()[0] * _current.getNQ()[DIM - 1]);
        }

        int cols() const
        {
            if (_transpose)
                return _A.rows() * (DIM == 1 ? _current.getNQ()[0] : _current.getNQ()[0] * _current.getNQ()[DIM - 1]);
            else
                return _A.cols() * (DIM == 1 ? _current.getNQ()[0] : _current.getNQ()[0] * _current.getNQ()[DIM - 1]);
        }

        // Access operator to allow for matrix free computation
        double operator()(const unsigned int &i, const unsigned int &j) const
        {

            if (_transpose)
            {

                if constexpr (DIM == 2)
                {

                    if (j % 2 == 0 /* col is even*/)
                    {
                        if (i % 2 == 0 /*row is even*/ && j == i /*on the diagonal*/)
                        {
                            return _A.coeff(0, 0);
                        }
                        else if (i % 2 == 1 /*row is */ && i == j + 1 /*on the lower diagonal*/)
                        {
                            return _A.coeff(0, 1);
                        }
                        else
                            return 0.;
                    }
                    else /*col is odd*/
                    {
                        if (i % 2 == 0 /*row is even*/ && i == j - 1 /*on the upper diagonal*/)
                        {
                            return _A.coeff(1, 0);
                        }
                        else if (i % 2 == 1 /*row is odd*/ && j == i /*on the diaglona*/)
                        {
                            return _A.coeff(1, 1);
                        }
                        else
                            return 0.;
                    }
                }
                else
                {
                    if (i == j /*on the diagonal*/)
                        return _A.coeff(0, 0);
                    else
                        return 0.;
                }
            }
            else
            {
                if constexpr (DIM == 2)
                {

                    if (i % 2 == 0 /* row is even*/)
                    {
                        if (j % 2 == 0 /*col is even*/ && i == j /*on the diagonal*/)
                        {
                            return _A.coeff(0, 0);
                        }
                        else if (j % 2 == 1 /*col is */ && j == i + 1 /*on the upper diagonal*/)
                        {
                            return _A.coeff(0, 1);
                        }
                        else
                            return 0.;
                    }
                    else /*row is odd*/
                    {
                        if (j % 2 == 0 /*col is even*/ && j == i - 1 /*on the lower diagonal*/)
                        {
                            return _A.coeff(1, 0);
                        }
                        else if (j % 2 == 1 /*col is odd*/ && j == i /*on the diaglona*/)
                        {
                            return _A.coeff(1, 1);
                        }
                        else
                            return 0.;
                    }
                }
                else
                {
                    if (i == j /*on the diagonal*/)
                        return _A.coeff(0, 0);
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
                    matrix(i, j) = this->operator()(i, j);
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
                    matrix(i, j) = this->operator()(j, i);
                }
            }
            return matrix;
        }
        //_______________________________ EIGEN SUPPORT ______________________________

        // Updater of the reference matrix
        void update(const Element<DIM> &elem)
        {

            _current = elem;
            if constexpr (DIM == 2)
            {
                _A.resize(DIM, DIM);
                _A = _current.template jacobian<Matrix2d>().inverse().transpose();
            }
            else
            {
                _A.resize(DIM, DIM);
                _A << 1 / _current.template jacobian<double>();
            }

            return;
        }

        const Derived &getRef() const
        {
            return _A;
        }

        // Destructor
        ~ExtendedJ() = default;

    private:
        Derived _A;
        Element<DIM> _current;
        const bool _transpose;
        static constexpr double zero_value = 0.0;
    };

    //==========================================================================================================

    // Finally define a class that will serve the general purpose of representing the
    // products between the previously defined classes:
    // To better define this class, use c++ concepts to specify some necessary
    // proprieties of the matrix involved in matrix product

    // PLEASE NOTE: The implementation of such proxy to allow for matrix free computation was thought to speed up the code
    //  as all the accesses to the elements of the product matrix would have been cached, but it turns out that simply managing the
    //  product between the previously defined custom classes by overloading the * operator by generating an intermidiate Eigen Matrix
    //  works far better in terms of performance.
    //  Still, I decided to leave the implementation of the GeneralProduct as a

    template <class T, class U>
        requires((IsGeneralMatrix<T> || IsDouble<T>) && IsGeneralMatrix<U>)
    class GeneralProduct
    {
    public:
        // Constructor
        GeneralProduct(const T &A, const U &B, const bool &transpose = false) : _A(A),
                                                                                _B(B),
                                                                                _transpose(transpose)
        {
            if constexpr (IsGeneralMatrix<T>)
                // check form correspondence of number of rows and cols
                if (A.cols() != B.rows())
                {
                    throw std::runtime_error("The number of cols of first operand must be equal to the number of cols of the second operand\n");
                }
        };

        // Definition of rows and cols methods
        int rows() const
        {
            if (_transpose)
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
            if (_transpose)
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
        double operator()(const unsigned int &i, const unsigned int &j) const
        {
            if (_transpose)
            {
                // Here define the different cases
                if constexpr (IsDouble<T>) // first operand is a double
                {
                    return _A * _B(j, i);
                }
                else // both operands are matrixes
                {
                    double sum(0);
                    for (unsigned int k = 0; k < _A.cols(); k++)
                    {
                        sum += _A(j, k) * _B(k, i);
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
                    return _A * _B(i, j);
                }
                else // both operands are matrixes
                {
                    double sum(0);
                    for (unsigned int k = 0; k < _A.cols(); k++)
                    {
                        sum += _A(i, k) * _B(k, j);
                    }
                    return sum;
                }
            }
        }

        ExtendMat::GeneralProduct<T, U> transpose() const
        {
            return ExtendMat::GeneralProduct<T, U>(_A, _B, !_transpose);
        }

        // Getters (for debugging)

        const T &getA() const
        {
            return _A;
        }

        const U &getB() const
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

    template <typename Derived>
    concept EigenMatrix = std::is_base_of_v<Eigen::EigenBase<Derived>, Derived>;

    // ========================= OVERLOADING OF OPERATOR * USING EIGEN CLASSES ===============================

    template <EigenMatrix Derived, unsigned short DIM>
    const SparseMatrix<double> operator*(const ExtendMat::ExtendedD<Derived, DIM> &D, const double &coeff)
    {
        SparseMatrix<double> result(D.rows(), D.cols());

        for (int i = 0; i < D.rows(); ++i)
            result.coeffRef(i, i) = coeff * D(i, i);

        return result;
    }

    template <EigenMatrix Derived1, EigenMatrix Derived2, unsigned short DIM>
    const MatrixXd operator*(const ExtendMat::ExtendedJ<Derived1, DIM> &J, const ExtendMat::ExtendedB<Derived2, DIM> &B)
    {
        MatrixXd result = MatrixXd::Zero(J.rows(), B.cols());

        for (int i = 0; i < J.rows(); ++i)
        {
            for (int j = 0; j < B.cols(); ++j)
            {
                for (int k = 0; k < J.cols(); ++k)
                {
                    result(i, j) += J(i, k) * B(k, j);
                }
            }
        }
        return result;
    }

    //  ========================= OVERLOADING OF OPERATOR * USING EIGEN CLASSES ===============================

    //  ========================    OPERATOR * OVERLOADING  ========================

}

#endif
