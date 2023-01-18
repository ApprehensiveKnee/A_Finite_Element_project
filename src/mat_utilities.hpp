#ifndef UTIL
#define UTIL

#include <numeric>
#include <algorithm>
#include <functional>
#include <tuple>
#include <type_traits>
#include "elements.hpp"


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
        // please observe this method is defined for 2D case only
        virtual double value(const Point &) const= 0;
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
            return (20*M_PI*M_PI + 1)*std::sin(2*M_PI*p.getX())+std::sin(4*M_PI*p.getY());
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
        void LGL_quadratures(const unsigned int &n/*number of nodes*/);

        // A member function to evaluate the nodes and weights of the Legendre Gauss Lobatto formulae on the interval [a, b];
        void LGL_quadratures(const unsigned int &n/*number of nodes*/,const double &a, const double &b);

        // standard getters
        std::vector<double> getN() const;
        std::vector<double> getW() const;

        ~Quadrature() = default;
    };



    template <class ElementType>
    class SpectralFE
    {
    private:

        
        //first, the element currently considered
        ElementType _current_elem;
        //to compute the quadrature coordinates and weights necessary
        Quadrature _qr;
        //the degree of the FE space considered  --> number of quadrature nodes per each direction of the element
        const unsigned int _r;
        //a vector containig the quadrature points of the elements
        std::vector<Point> _quad_points;
        //a vector containing the quadrature wights for the quadrature points
        std::vector<double> _quad_weights;
        //a vector containing the nodes of the element
        std::vector<Point> _int_nodes;
        // a counter to specify the element of the mesh currently examined
        unsigned int pun;


    public:
        //constructor
        SpectralFE(const unsigned int &r)
            :_current_elem(),
            _qr(),
            _r(r),
            _quad_points(),
            _quad_weights(),
            _int_nodes(),
            pun(0)
            {};

        SpectralFE(const ElementType& element, const unsigned int &r)
            :_current_elem(element),
            _qr(),
            _r(r),
            _quad_points(),
            _quad_weights(),
            _int_nodes(),
            pun(0)
            {};
        //standar getters
        const ElementType& getCurrent() const
        {
            return _current_elem;
        };

        const Quadrature& getQuad() const
        {
            return _qr;
        };

        const unsigned int& getDeg() const
        {
            return _r;
        };

        const std::vector<Point>& getQPoints() const
        {
            return _quad_points;
        };


        const std::vector<Point>& getINodes() const
        {
            return _int_nodes;
        };

        const std::vector<double>& getWeights() const
        {
            return _quad_weights;
        };

        const unsigned int& getPun() const
        {
            return pun;
        };

        
        // a method to compute the number of degree of freedom per element
        unsigned int dof_per_cell() const
        {   
            if constexpr (DIM == 2)
                return _current_elem.internalN(this->getDeg())*_current_elem.internalN(this->getDeg());
            else
                return _current_elem.internalN(this->getDeg());
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
            _update_quad();
            _update_nodes();
            return;
        };

        
        // member function to map the quadrature points on to the specific element
        // 2D CASE
        Point quadrature_point(const unsigned int &q) const
        {
            if constexpr (DIM == 2)
            {
                double x, y;
                std::tie(x,y) = _current_elem.direct_map((this->getQPoints()[q]).getX(),(this->getQPoints()[q]).getY());
                Point my_point(x,y);
                return my_point;
            }
            else
            {
                double x;
                x = _current_elem.direct_map((this->getQPoints()[q]).getX());
                Point my_point(x);
                return my_point;
            }
            
        };
        
        // member functions to compute the basis functions and gradient
        
        double shape_value(const unsigned int &i /*index of the quadrature node*/,
                        const unsigned int &q /*index of the quadrature point of which we compute the evaluation*/) const
        {
            if constexpr (DIM == 2)
            {
                //first we get the node of the element represented by the index i, using the internal standard node numbering:
                // we do so by looking into the _int_nodes member, containing the coordiates of all the internal nodes of the current element
                // and we get the funtion to evaluate the basis function in a certain point (quadrature points)
                auto phi_i =FETools::SpectralFE<ElementType>::phi((this->getINodes())[i].getX(),(this->getINodes())[i].getY(),this->getDeg()+1);

                //now we evaluate the basis function over the quadrature point

                return phi_i((this->getQPoints())[q].getX(),(this->getQPoints())[q].getY()) > tol? phi_i((this->getQPoints())[q].getX(),(this->getQPoints())[q].getY()): 0;
            }
            else
            {
                //corresponding 1D case using lagrangiang basis functions
                auto phi_i =FETools::SpectralFE<ElementType>::lagrange_pol((this->getINodes())[i].getX(),this->getDeg()+1);
                return phi_i((this->getQPoints())[q].getX()) > tol? phi_i((this->getQPoints())[q].getX()):0;
            }
            
        }

        VectorXd shape_grad(const unsigned int &i /*index of the quadrature node*/,
                            const unsigned int &q /*index of the quadrature point of which we compute the evaluation*/) const
        {
            if constexpr (DIM == 2)
            {
                //first we get the node of the element represented by the index i, using the internal standard node numbering:
                // we do so by looking into the _int_nodes member, containing the coordiates of all the internal nodes of the current element
                // we then get the funtion to evaluate the derivative of the basis function in a certain point (quadrature points)
                auto dphi_i =FETools::SpectralFE<ElementType>::dphi((this->getINodes())[i].getX(),(this->getINodes())[i].getY(),this->getDeg()+1);

                // now we evaluate the partial derivatives of the basis function over the quadrature point

                VectorXd grad(2);
                grad[0]= dphi_i((this->getQPoints())[q].getY(),(this->getQPoints())[q].getX())>inf?dphi_i((this->getQPoints())[q].getY(),(this->getQPoints())[q].getX()):inf, // d_phi(x_p,y_p)/dx
                grad[1]= dphi_i((this->getQPoints())[q].getX(),(this->getQPoints())[q].getY())>inf?dphi_i((this->getQPoints())[q].getX(),(this->getQPoints())[q].getY()):inf; // d_phi(x_p,y_p)/dy
                //return grad;
                //se ritornassimo direttamente grad, questo corrisponderebbe a una valutazione della derivata della funzione di base
                //sull'elemento di riferimento. Ma bisogna poi ancora considerare la mappatura sull'elemento specifico
                if(grad.norm()>tol)
                    return (this->_J_inv_t())*grad;
                else 
                    return VectorXd::Constant(2, 0);

            }
            else
            {
                //corresponding 1D case
                auto dphi_i = FETools::SpectralFE<ElementType>::lagrange_der((this->getINodes())[i].getX(),this->getDeg()+1);
                VectorXd grad(1);
                grad[0]= dphi_i((this->getQPoints())[q].getX())?dphi_i((this->getQPoints())[q].getX()):inf; // dphi(x_q)/dx
                if(grad.norm() > tol)
                    return (this->_J_inv_t())*grad;
                else 
                    return VectorXd::Constant(1, 0);
            }
            
        };
        //member functions to compute the product of J and quadrature weights
        double JxW(const unsigned int &q /*quadratue point associated to the weight*/) const
        {
            //firstly, let's compute the determinant of the jacobian
            return (this->getWeights()[q]*(this->_det_J()));
            
        };
        //default destructor
        ~SpectralFE() = default;

    //private:
        
        //firstly, a method to compute the coordinates and  weithts of the quadrature points along a certain dimension
        //used to approximate integrals on the current element
        std::vector<double> _comp_quad_c()
        {
            _qr.LGL_quadratures(this->getDeg()+1);
            return _qr.getN();
        };

        std::vector<double> _comp_quad_w()
        {
            _qr.LGL_quadratures(this->getDeg()+1);
            return _qr.getW();
        };

        // a method to update the vector of quadrature points for the current element
        // 2D CASE
        void _update_quad()
        {
            // delete the quadrature coordinates of the past element
            // OSS. in this case, the degree of the FE spaces along the two dimensions
            // are supposed to be the same and identical to _r member. Furthermore, 
            //for each element, the number of quadrature points is supposed to be the same:
            // the quadrature points are computed over the refence element [-1,1]*[-1,1],
            // as such the updated quadrature points and weights will always be the same for every element--> REDUNDANT WORK.
            // Yet, in the case we would like to change the degree of the quadrature formula over a
            // specific element, this function would make this process possible and less difficult to implement
            if constexpr (DIM == 2)
            {

                {
                    _quad_points.clear();
                    auto temp = this->_comp_quad_c();

                    //loop over y coordinate
                    for(auto j : temp)
                    {
                        //loop over x coordinate
                        for(auto i : temp)
                        {
                            Point my_point(i,j);
                            _quad_points.emplace_back(my_point);
                        }
                    }

                }
                
                {
                    _quad_weights.clear();
                    std::vector<double> temp = this->_comp_quad_w();
                    //loop over the y coordinate
                    for(auto j : temp)
                    {
                        //loop over x coordinate
                        for(auto i : temp)
                        {
                            _quad_weights.emplace_back(i*j);
                        }
                    }

                }

            }
            else
            {

                {
                    _quad_points.clear();
                    auto temp = this->_comp_quad_c();

                    //loop over x coordinate
                    for(auto i : temp)
                    {
                        Point my_point(i);
                        _quad_points.emplace_back(my_point);
                    }

                }
                
                {
                    _quad_weights.clear();
                    std::vector<double> temp = this->_comp_quad_w();
                    //loop over x coordinate
                    for(auto i : temp)
                    {
                        _quad_weights.emplace_back(i);
                    }

                }

            }
            

            return;
        };
        // a method to update the vector of nodes within the current element;
        // 2D CASE
        void _update_nodes()
        {
            if constexpr(DIM ==2)
            {

                //delete the nodes of the past element
                _int_nodes.clear();
                //compute the internal nodes

                // ---->  _current_elem.getNodes()[0].printNode();
                auto internal_nodes_x = _current_elem.getNodes()[0].nodes(_r, _current_elem.getNodes()[1]);
                auto internal_nodes_y = _current_elem.getNodes()[0].nodes(_r, _current_elem.getNodes()[2]);
                //loop over y_coordinate
                for( auto j : internal_nodes_y)
                {
                    //loop over x coordinate
                    for( auto i : internal_nodes_x)
                    {
                        // since we need the coordinates of the nodes over the reference element,
                        // for each node we apply the inverse affine transformation
                        double r,s;
                        std::tie(r,std::ignore) = _current_elem.inverse_map(i[0],i[1]);
                        if(std::abs(r) < tol) r= 0;
                        std::tie(std::ignore, s) = _current_elem.inverse_map(j[0],j[1]);
                        if(std::abs(s) < tol) s = 0;
                        Point my_point(r,s);
                        _int_nodes.emplace_back(my_point);
                    }
                }

            }
            else
            {

                //delete the nodes of the past element
                _int_nodes.clear();
                //compute the internal nodes

                // ---> _current_elem.getNodes()[0].printNode();
                auto internal_nodes_x = _current_elem.getNodes()[0].nodes(_r, _current_elem.getNodes()[1]);
                //loop over x_coordinate
                for( auto i : internal_nodes_x)
                {
                    
                    // since we need the coordinates of the nodes over the reference element,
                    // for each node we apply the inverse affine transformation
                    double r = _current_elem.inverse_map(i[0]);
                    Point my_point(r);
                    _int_nodes.emplace_back(my_point);

                }
            }
            
            return;
        };


        
        //some methods to compute the algebraic objects necessary to compute the elements of the local matrix
        // - the jacobian
        // - its determinant
        // - the inverse of the jacobian
        // - the inverse transpose of the jacobian

        MatrixXd _J() const
        {
            if constexpr(DIM == 2)
                return _current_elem.jacobian();
            else
            {
                MatrixXd J(1,1);
                J << _current_elem.jacobian();
                return J;

            }

            
        };
        double _det_J() const
        {
            if constexpr(DIM == 2)
                return _current_elem.jacobian().determinant();
            else
                return _current_elem.jacobian();
            
        }
        MatrixXd _J_inv() const
        {
            if constexpr(DIM == 2)
                return _current_elem.jacobian().inverse();
            else
            {
                MatrixXd J(1,1);
                J << _current_elem.jacobian();
                J.inverse();
                return J;

            }
            
        }
        MatrixXd _J_inv_t() const
        {
            if constexpr(DIM == 2)
                return _current_elem.jacobian().inverse().transpose();
            else
            {
                MatrixXd J(1,1);
                J << _current_elem.jacobian();
                J.inverse();
                // transposition doesn't affect the matrix
                return J;

            }
            
        }

        
        
        //function for legendre polynomilas and first derivative (iterative)
        static std::tuple<std::function <double(const double &)> ,std::function <double(const double &)>>legendre_pol(const unsigned int &n)
        {
            //legendre polynomial of degree 1 and derivative
            std::function <double(const double &)> l1 = [](const double &x){return x;};
            std::function <double(const double &)> l1d = [](const double &x){return 1.;};
            //legendre polynomial of degree 0 and derivative
            std::function <double(const double &)> l0 = [](const double &x){return 1.;};
            std::function <double(const double &)> l0d = [](const double &x){return 0.;};

            std::function <double(const double &)> ln;
            std::function <double(const double &)> lnd;
            if(n == 0){
                return {l0,l0d};
            }
            else if(n == 1)
            {
                return {l1,l1d};
            }
            for(unsigned int k = 1; k < n; ++k){
                //iteratively contrusct the legendre polynomial of order k +1 <= n (Ln(x)) and derivative ((Ln)'(x))
                ln = [k, l1, l0](const double &x){return (((2*k + 1)*x*l1(x) - k*l0(x))/(k+1));};
                lnd = [k, l1, l1d, l0d](const double &x){return (((2*k + 1)*(x*l1d(x) + l1(x)) - k*l0d(x))/(k+1));};

                //update the previous polynomials and correspective derivatives
                l0 = [l1](const double &x){return l1(x);};
                l0d = [l1d](const double &x){return l1d(x);};
                l1= [ln](const double &x){return ln(x);};
                l1d= [lnd](const double &x){return lnd(x);};
                
            }
            return {ln,lnd};
                
        };

        //function for legendre polynomilas, first and second order derivative (iterative)
        static std::tuple<std::function <double(const double &)> ,std::function <double(const double &)>, std::function <double(const double &)>>legendre_der(const unsigned int &n)
        {
            //legendre polynomial of degree 1 and derivative
            std::function <double(const double &)> l1 = [](const double &x){return x;};
            std::function <double(const double &)> l1d = [](const double &x){return 1.;};
            std::function <double(const double &)> l1dd = [](const double &x){return 0.;};
            //legendre polynomial of degree 0 and derivative
            std::function <double(const double &)> l0 = [](const double &x){return 1.;};
            std::function <double(const double &)> l0d = [](const double &x){return 0.;};
            std::function <double(const double &)> l0dd = [](const double &x){return 0.;};

            std::function <double(const double &)> ln;
            std::function <double(const double &)> lnd;
            std::function <double(const double &)> lndd;
            if(n == 0){
                return {l0,l0d,l0dd};
            }
            else if(n == 1)
            {
                return {l1,l1d,l1dd};
            }
            for(unsigned int k = 1; k < n; ++k){
                //iteratively contrusct the legendre polynomial of order k +1 <= n (Ln(x)) and derivative ((Ln)'(x))
                ln = [k, l1, l0](const double &x){return (((2*k + 1)*x*l1(x) - k*l0(x))/(k+1));};
                lnd = [k, l1, l1d, l0d](const double &x){return (((2*k + 1)*(x*l1d(x) + l1(x)) - k*l0d(x))/(k+1));};
                lndd = [k, l1d, l0dd, l1dd](const double &x){return (((2*k + 1)*(x*l1dd(x) + 2*l1d(x)) - k*l0dd(x))/(k+1));};

                //update the previous polynomials and correspective derivatives
                l0 = [l1](const double &x){return l1(x);};
                l0d = [l1d](const double &x){return l1d(x);};
                l0dd = [l1dd](const double &x){return l1dd(x);};
                l1= [ln](const double &x){return ln(x);};
                l1d= [lnd](const double &x){return lnd(x);};
                l1dd = [lndd](const double &x){return lndd(x);};
                
            }
            return {ln,lnd, lndd};
        };

        //function for Lagrange polynomials (1D) on the reference element [-1,1]
        static std::function <double(const double &)> lagrange_pol(const double &xi/*the index of the basis function, i.e. the index of the node on which the basis function isn't zero*/,
                                                                    const unsigned int &n/*the degree of the polynomial function, i.e the degree of the FEm space considered*/)
        {
            //first, we get the legendre polynomials of degree n-1
            auto [ln,lnd] = FETools::SpectralFE<ElementType>::legendre_pol(n-1);
            //then generate the function to evaluate the Lagrangian basis for a given coordinate
            std::function<double(const double &)> l_i = [&xi, &n, lnd, ln](const double &x){return (-1./((n-1)*(n))*(1. - x*x)*lnd(x)/((x - xi)*ln(xi)));};
            return l_i;
        };

        //function for Lagrange polynomials first order derivative (1D) on the reference element [-1,1]
        static std::function <double(const double &)> lagrange_der(const double &xi/*the index of the basis function, i.e. the index of the node on which the basis function isn't zero*/,
                                                                    const unsigned int &n/*the degree of the polynomial function, i.e the degree of the FEm space considered*/)
        {
            //first, we get the legendre polynomials of degree n-1 and derivative of first order and second order
            auto [ln,lnd, lndd] = FETools::SpectralFE<ElementType>::legendre_der(n-1);
            //then generate the function to evaluate the derivate Lagrangian basis for a given coordinate
            // OSS. for x = xi the function is not defined (-> +inf): we will handle this points by returning value 1
            std::function<double(const double &)> dl_i = [&xi, &n, lndd ,lnd, ln](const double &x)
            {
                if(x == -1.)
                {
                    double x_nex = x + std::numeric_limits<double>::epsilon();
                    return (-1./(n*(n-1)*ln(xi)))*((-(1 + x_nex*x_nex - 2*x_nex*xi)*lnd(x) + (1-x_nex*x_nex)*(x_nex-xi)*lndd(x_nex))/((x_nex-xi)*(x_nex-xi)));
                }
                else if(x == 1.)
                {
                    double x_nex = x - std::numeric_limits<double>::epsilon();
                    return (-1./(n*(n-1)*ln(xi)))*((-(1 + x_nex*x_nex - 2*x_nex*xi)*lnd(x) + (1-x_nex*x_nex)*(x_nex-xi)*lndd(x_nex))/((x_nex-xi)*(x_nex-xi)));

                }
                else if(x ==xi)
                {
                    return 0.;
                }
                else
                {
                    return (-1./(n*(n-1)*ln(xi)))*((-(1 + x*x - 2*x*xi)*lnd(x) + (1-x*x)*(x-xi)*lndd(x))/((x-xi)*(x-xi)));
                }
                    
            };
            
            return dl_i;

        };


        static std::tuple<std::function<double(const double &)>,std::function<double(const double &)>> lagrange_interp(const double &xi/*the index of the basis function, i.e. the index of the node on which the basis function isn't zero*/,
                                                                    const unsigned int &n/*the degree of the polynomial function, i.e the degree of the FEm space considered*/)
        {

            //fist of all, determine the cooridinates of the internal nodes on the 1D refernce interval [-1,1];
            std::vector<double> eps(n);
            for(unsigned int  i = 0; i < n; i++)
            {
                eps[i] = -1. + i*(2./(n-1));
            }
            // now iteratively contruct the lagrange interpolating polynomials
            std::function <double(const double &)> ln = [](const double &x){return 1.;};
            std::function <double(const double &)> lnd = [](const double &x){return 0.;};
            for(unsigned int j = 0; j < n; j++)
            {

                std::function <double(const double &)> l_j;
                if(eps[j]!= xi)
                {
                    l_j = [&xi, eps, j](const double &x){return (x-eps[j])/(xi -eps[j]);};
                    std::function <double(const double &)> l_m = [&xi, eps, j](const double &x){return 1./(xi -eps[j]);};
                    ln = [l_j,ln](const double &x){return ln(x)*l_j(x);};
                    lnd = [l_m, lnd](const double &x){return lnd(x)+l_m(x);};
                }
                
            }
            ln = [ln,&xi](const double &x){return x!=xi?ln(x):1.;};
            lnd = [ln,lnd, &xi](const double &x){
                                                if(x!=-1. || x!=1.)
                                                    return x!=xi?lnd(x)*ln(x):0.;
                                                else if (x == -1.)
                                                {
                                                    double x_nex = x + std::numeric_limits<double>::epsilon();
                                                    return lnd(x_nex)*ln(x_nex);
                                                }
                                                else
                                                {
                                                    double x_nex = x - std::numeric_limits<double>::epsilon();
                                                    return lnd(x_nex)*ln(x_nex);
                                                } 
                                            };
            
            return {ln,lnd};
        }
        // function for Lagrange polynomials (2D basis function l_i(x)*l_j(y)) on the reference element [-1,1]*[-1,1]
        static std::function <double(const double &, const double &)> phi(const double &xi, 
                                                                            const double &xj,
                                                                            const unsigned int &n)
        {
            //lagrangian basis function for the two dimension
            std::function <double(const double &)> l_i;
            std::tie(l_i, std::ignore)= FETools::SpectralFE<ElementType>::lagrange_interp(xi,n);
            std::function <double(const double &)> l_j;
            std::tie(l_j,std::ignore)= FETools::SpectralFE<ElementType>::lagrange_interp(xj,n);
            //the basis funcion for the 2D case is given as tensor product of the two lagrangian polynomials
            std::function <double(const double&, const double &)> phi = [l_i, l_j](const double &x, const double &y){ return l_i(x)*l_j(y);};

            return phi;
        }

        // function for Lagrange polynomials partial derivative (2D basis function l_i(x)*dl_j(y)) on the reference element [-1,1]*[-1,1] 
        // with respect to the y variable
        static std::function <double(const double &, const double &)> dphi(const double &xi, 
                                                                            const double &xj,
                                                                            const unsigned int &n)
        {
            //lagrangian basis function for the two dimension
            std::function <double(const double &)> l_i;
            std::tie(l_i, std::ignore )= FETools::SpectralFE<ElementType>::lagrange_interp(xi,n);
            std::function <double(const double &)> dl_j;
            std::tie(std::ignore, dl_j)= FETools::SpectralFE<ElementType>::lagrange_interp(xj,n);
            //the basis funcion for the 2D case is given as tensor product of the two lagrangian polynomials
            std::function <double(const double&, const double &)> phi = [l_i, dl_j](const double &x, const double &y){ return l_i(x)*dl_j(y);};

            return phi;
        }

        
    };
}

#endif