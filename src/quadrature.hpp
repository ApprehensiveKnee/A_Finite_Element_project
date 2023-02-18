

//=================================================================
//              HEADER FILE FOR THE QUADRATURE CLASS
//=================================================================


#ifndef QUAD
#define QUAD


#include <numeric>
#include <algorithm>
#include <functional>
#include <tuple>
#include <type_traits>
#include "elements.hpp"



//______________________________________________________________

namespace FETools
{
    // a class to implement all the functions necessary to evaluate the Gauss Legendre Lobotto quadrature nodes and members
    class Quadrature
    {
    private:
        std::vector<double> _nodes;
        std::vector<double> _weights;
    public:
        //standard contructor
        Quadrature()
            :_nodes(),
            _weights()
            {};

        //A member function to evaluate the Jacobi polinomial and derivatives at x in [-1,1] given the degree n and the paramethers alpha and beta
        static void jacobi_pol(const std::vector<double> &x, const unsigned int &n, const double &_alpha, const double &_beta, std::array<std::vector<double>,3>& pol, std::array<std::vector<double>,3>& der);


        template <typename Iterator>
        // A member function to evaluate the n roots of the Jacobi polynomial, obtained by using Newton method and deflation process.
        static void jacobi_roots(const unsigned int &n, const double &_alpha, const double &_beta, const Iterator start, const Iterator end)
        {
            if( n < 1)
            {
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
                    if(std::isnan(*(it)) || std::abs(*(it)) < tol){ *(it) = 0;}
                }
                std::sort(start, end);
            }
            return;
        };


        // A member function to evaluate Legendre polynomials of degree n at x coordinates, using the three term relation
        static void legendre_pol(const std::vector<double> &x, const unsigned int &n, std::vector<double> &pol);


        // A member function to evaluate the nodes and weights of the Legendre Gauss Lobatto formulae on the interval [-1, 1];
        // the solution is then stored inside the _nodes and _weights members
        void LGL_quadratures(const unsigned int &n/*number of nodes*/);

        // A member function to evaluate the nodes and weights of the Legendre Gauss Lobatto formulae on the interval [a, b];
        void LGL_quadratures(const unsigned int &n/*number of nodes*/,const double &a, const double &b);



        // standard getters
        const std::vector<double>& getN() const;
        const std::vector<double>& getW() const;

        ~Quadrature() = default;
    };



}

#endif

