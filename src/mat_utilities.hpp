
//=================================================================
// HEADER FILE FOR THE CLASS QUADRATURE AND FE SPECTRAL DEFINITIONS
//=================================================================

#ifndef UTIL
#define UTIL

#include <numeric>
#include <algorithm>
#include <functional>
#include <tuple>
#include <type_traits>
#include "elements.hpp"
#include "mesh.hpp"


//______________________________________________________________


namespace Functions
{
    
    //Pure class for scalar function taking values in 2D domains
    class Function
    {
    public:
        //standard contructor
        Function(){};
        //define a virtual method that will be overridden by the derived classes
        virtual double value(const Point &) const= 0;
        virtual std::array<double, DIM> grad(const Point&) const{return {0.,0.};};
        //define virtual contructor
        virtual ~Function(){};

    };

    //then some general functions derived from the function base class
    //we hereby consider just scalar functions

    class DiffusionCoefficient : public Function
    {
    public:
        // Constructor.
        DiffusionCoefficient()
        {}

        // Evaluation.
        virtual double
        value(const Point & /*p*/) const override
        {
        return 1.0;
        }

    };

    // Reaction coefficient.
    class ReactionCoefficient : public Function
    {
    public:
        // Constructor.
        ReactionCoefficient()
        {}

        // Evaluation.
        virtual double
        value(const Point& /*p*/) const override
        {
        return 1.0;
        }
    };

    // Forcing term.
    class ForcingTerm : public Function
    {
    public:
        // Constructor.
        ForcingTerm()
        {}

        // Evaluation.
        virtual double
        value(const Point & p) const override
        {
            return (20*M_PI*M_PI + 1)*std::sin(2*M_PI*p.getX())*std::sin(4*M_PI*p.getY());
        }
    };

    // The exact solution (must be known)
    class ExactSolution : public Function
    {
    public:
        // Constructor.
        ExactSolution()
        {}

        // Evaluation.
        virtual double
        value(const Point & p) const override
        {
            return std::sin(2*M_PI*p.getX())*std::sin(4*M_PI*p.getY());
        }

        std::array<double, DIM> 
        grad(const Point &p) const override
        {
            if constexpr(DIM == 2)
                return{2 * M_PI *std::cos(2*M_PI*p.getX())*std::sin(4*M_PI*p.getY()), 4 * M_PI * std::sin(2*M_PI*p.getX())*std::cos(4*M_PI*p.getY())};
            else
                return{2 * M_PI *std::cos(2*M_PI*p.getX())*std::sin(4*M_PI*p.getY())};
        }
    };

    // Some constant functions to use, for instance, with Dirichelet b.c

    class functionZero : public Function
    {
    public:
        // Constructor.
        functionZero()
        {}

        // Evaluation.
        virtual double
        value(const Point & /*p*/) const override
        {
        return 0.0;
        }
    };

    class functionOne : public Function
    {
    public:
        // Constructor.
        functionOne()
        {}

        // Evaluation.
        virtual double
        value(const Point & /*p*/) const override
        {
        return 1.0;
        }
    };
}




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





