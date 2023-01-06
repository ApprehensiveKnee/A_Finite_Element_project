#ifndef UTIL
#define UTIL

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <array>
#include <iterator>
#include <functional>
#include <tuple>
#include "elements.hpp"


//A simple function to compute the difference between two diffent points in any diemensio
double Difference(const std::vector<double> &p1, const std::vector<double> &p2, std::ostream &out = std::cout);



//______________________________________________________________


//a class to store the coordinates of the relevant points (quadrature coordinates)
class Point
{
private:
    std::array<double, DIM> _coord;
public:
    Point(const double &x, const double &y = 0):
        _coord(classInit(x,y))
        {};
    //standard getter
    std::tuple<double,double> getCoord() const;
    //a function to print the point
    void printPoint(std::ostream &out = std::cout) const;
    //a static function to initialize the class
    static std::array<double,DIM> classInit(const double &x, const double &y);
    
    ~Point() = default;

};


//______________________________________________________________

//Pure class to implement the evaluation of the function considered for the problem---> interface
class Function
{
public:
    //standard contructor
    Function(){};
    //define a virtual method that will be overridden by the derived classes
    // please observe this method is defined for 2D case only
    virtual double value(const Point p) = 0;
    //define virtual contructor
    virtual ~Function(){};

};


//______________________________________________________________


// a class to implement all the functions necessary to evaluate the Gauss Legendre Lobotto quadrature nodes and members
class Quadrature
{
private:
    std::vector<double> _nodes;
    std::vector<double> _weights;
public:
    //standard contructor
    Quadrature(){};

    //A member function to evaluate the Jacobi polinomial and derivatives at x in [-1,1] given the degree n and the paramethers alpha and beta
    static void jacobi_pol(const std::vector<double> &x, const unsigned int &n, const double &_alpha, const double &_beta, std::array<std::vector<double>,3>& pol, std::array<std::vector<double>,3>& der);


    template <typename Iterator>
    // A member function to evaluate the n roots of the Jacobi polynomial, obtained by using Newton method and deflation process.
    static void jacobi_roots(const unsigned int &n, const double &_alpha, const double &_beta, const Iterator start, const Iterator end)
    {
        if( n < 1)
        {
            std::cout << "The polynomial degree should be greater than 0" << std::endl;
            return;
        }
        else
        {
            // assume the vector given has already the needed space allocated,
            // no extra control is requested
            std::vector<double> x0 (1,std::cos(M_PI/(2*n)));
            double x1(0.);
            double tol=1e-14; unsigned int kmax=15;
            //for(unsigned int j = 0; j < n; ++j)
            unsigned int j = 1;
            for (Iterator it = start; it !=end; ++it, j++)
            {   
                
                double diff = tol + 1;
                for(unsigned int kiter = 0; kiter <= kmax && diff >= tol; kiter++)
                {

                    std::array<std::vector<double>,3> pol, der;

                    jacobi_pol(x0, n, _alpha, _beta, pol, der);
                    // deflation process q(x)=p(x)*(x-x_1)*... (x-x_{j-1})
                    // q(x)/q'(x)=p(x)/[p'(x)-p(x)*\sum_{i<j} 1/(x-x_i)]
                    double ss = std::accumulate(start,
                                                it, 
                                                1./x0[0],
                                                [&](const double &a, const double &b){return a + 1./ (x0[0]-b);});
                    x1 = x0[0]-pol[0][0]/(der[0][0]-ss*pol[0][0]);
                    double diff = std::abs(x1-x0[0]);
                    x0[0]=x1;
                }
                x0[0] = .5 * (x1+std::cos((2*(j+1)-1)*M_PI/(2*n)));
                *it = x1;
        
            }
            if(std::isnan(*(start+(j+1)/2))){ *(start+(j+1)/2) = 0;}
            std::sort(start, end);
        }
        return;
    };


    // A member function to evaluate Legendre polynomials of degree n at x coordinates, using the three term relation
    static void legendre_pol(const std::vector<double> &x, const unsigned int &n, std::vector<double> &pol);

    // A member function to evaluate the nodes and weights of the Legendre Gauss Lobatto formulae on the interval [-1, 1];
    void LGL_quadratures(const unsigned int &n/*number of nodes*/);

    // A member function to evaluate the nodes and weights of the Legendre Gauss Lobatto formulae on the interval [a, b];
    void LGL_quadratures(const unsigned int &n/*number of nodes*/,const double &a, const double &b);

    // standard getters
    std::vector<double> getN() const;
    std::vector<double> getW() const;

    ~Quadrature() = default;
};

//______________________________________________________________

class SpectralFE
{
private:

    //first, the element currently considered
    Rectangle _current_elem;
    //to compute the quadrature coordinates and weights necessary
    Quadrature _qr;
    //the degree of the FE space considered  --> number of quadrature nodes per each direction of the element
    unsigned int _r;
    //a vector containig the quadrature points of the elements
    std::vector<Point> _quad_points;
    //a vector containing the nodes of the element
    std::vector<Point> _int_nodes;
    

public:
    //constructor
    SpectralFE(const unsigned int &r)
        :_r(r)
        {};

    SpectralFE(const Rectangle& element, const Quadrature& qr, const unsigned int &r)
        :_current_elem(element),
        _qr(qr),
        _r(r)
        {};
    //standar getters
    Rectangle getCurrent() const;
    Quadrature getQuad() const;
    unsigned int getDeg() const;
    std::vector<Point> getQPoints() const;
    std::vector<Point> getINodes() const;

    
    // a method to compute the number of degree of freedom per element
    unsigned int _dof_per_cell();

    

    void _update_current( const Rectangle geoele );

    // member functions to compute the basis functions and gradient
    double shape_value( const unsigned int &i /*index of the quadrature node*/,
                        Point q /*the quadrature point of which we compute the evaluation*/);
    Vector2d shape_grad(const unsigned int &i /*index of the quadrature node*/,
                         Point q /*the quadrature point of which we compute the evaluation*/);

protected:
    //firstly, a method to compute the coordinates and  weithts of the quadrature points along a certain dimension
    //used to approximate integrals on the current element
    std::vector<double> _comp_quad_c();
    std::vector<double> _comp_quad_w();
    // a method to update the vector of quadrature points for the current element
    void _update_quad();
    // a method to update the vector of nodes within the current element;
    void _update_nodes();

    
    //some methods to compute the algebraic objects necessary to compute the elements of the local matrix
    // - the jacobian
    // - its determinant
    // - the inverse of the jacobian
    // - the inverse transpose of the jacobian

    Matrix2d J();
    double detJ();
    Matrix2d J_inv();
    Matrix2d J_inv_t();
    
    
    //function for legendre polynomilas and first derivative (iterative)
    static std::tuple<std::function <double(const double &)> ,std::function <double(const double &)>>legendre_pol(const unsigned int &n);
    //function for legendre polynomilas, first and second order derivative (iterative)
    static std::tuple<std::function <double(const double &)> ,std::function <double(const double &)>, std::function <double(const double &)>>legendre_der(const unsigned int &n);
    //function for Lagrange polynomials (1D) on the reference element [-1,1]
    static std::function <double(const double &)> lagrange_pol(const double &xi/*the index of the basis function, i.e. the index of the node on which the basis function isn't zero*/,
                                                                 const unsigned int &n/*the degree of the polynomial function, i.e the degree of the FEm space considered*/);

    //function for Lagrange polynomials first order derivative (1D) on the reference element [-1,1]
    static std::function <double(const double &)> lagrange_der(const double &xi/*the index of the basis function, i.e. the index of the node on which the basis function isn't zero*/,
                                                                 const unsigned int &n/*the degree of the polynomial function, i.e the degree of the FEm space considered*/);

    // function for Lagrange polynomials (2D basis function l_i(x)*l_j(y)) on the reference element [-1,1]*[-1,1]
    static std::function <double(const double &, const double &)> phi(const double &xi, 
                                                                         const double &xj,
                                                                         const unsigned int &n);

    // function for Lagrange polynomials partial derivative (2D basis function l_i(x)*dl_j(y)) on the reference element [-1,1]*[-1,1] 
    // with respect to the y variable
    static std::function <double(const double &, const double &)> dphi(const double &xi, 
                                                                         const double &xj,
                                                                         const unsigned int &n);

    
    //! compute phiDer
    // void _comp_phiDer();
    
};
#endif