#include "mat_utilities.hpp"


double Difference(const std::vector<double> &p1, const std::vector<double> &p2, std::ostream &out)
{
    if(p1.size() != p2.size()){
        out << "ERROR: the points do not have the same size";
        return 0.;
    }
    double temp(0);
    for (unsigned int i = 0; i < p1.size(); ++i){
        temp += (p2[i]-p1[i])*(p2[i]-p1[i]); // (x1-x2)^2 + (y1-y2)^2 + ...
    }
    return std::sqrt(temp);
};

//_____________________ POINT CLASS MEMEBERS___________________

std::tuple<double,double> Point::getCoord() const
{
    return {_coord[0],_coord[1]};
};


void Point::printPoint(std::ostream &out) const
{
    auto [temp1,temp2] = this->getCoord();
    out << temp1 << ", " << temp2 << std::endl;
    return;
};

std::array<double, DIM> Point::classInit(const double &x, const double &y)
{
    std::array<double, DIM> temp;
    if(DIM == 1)
    {
        temp[0] = x;
        return temp;
    }
    else if(DIM == 2)
    {
        temp[0] = x;
        temp[1] = y;
        return temp;
    }
};

//______________________________________________________________


//_________________________ QUADRATURE CLASS MEMBERS ____________________________


void Quadrature::jacobi_pol(const std::vector<double> &x, const unsigned int &n, const double &_alpha, const double &_beta, std::array<std::vector<double>,3>& pol, std::array<std::vector<double>,3>& der)
{

    

    auto apb = _alpha+_beta;
    auto ab2 = (_alpha*_alpha) - (_beta*_beta);


    //we guess the vectors passed as arguments are newly created 
    //and not already initialized, so we reserve the required space and 
    //the procede to fill them
    pol[0].resize(x.size(), 1.);
    der[0].resize(x.size(), 0.);

    pol[1].resize(x.size(), 0.);
    der[1].resize(x.size(), 0.);

    pol[2].resize(x.size(), 0.);
    der[2].resize(x.size(), 0.);

    //then we fill the vectors with the initial values
    // std::iota(pol[0].begin(), pol[0].end(), 1.);
    // std::iota(der[0].begin(), der[0].end(), 0.);

    // std::iota(pol[1].begin(), pol[1].end(), 0.);
    // std::iota(der[1].begin(), der[1].end(), 0.);

    // std::iota(pol[2].begin(), pol[2].end(), 0.);
    // std::iota(der[2].begin(), der[2].end(), 0.);

    
    if(n == 0)
    {
        return;
    }
    else if(n == 1)
    {
        pol[1] = pol[0];
        pol[2] = pol[1];

        der[1] = der[0];
        der[2] = der[1];

        for(unsigned int i = 0; x.size(); ++i){
            pol[0][i] = (_alpha - _beta+(apb+2)*x[i])*0.5;
            der[0][i] = 0.5*(apb+2);
        }
        return;
    }
    else
    {
        pol[1] = pol[0];
        pol[2] = pol[1];

        der[1] = der[0];
        der[2] = der[1];

        

        for(unsigned int i = 0; i <x.size(); ++i){
            pol[0][i] = (_alpha - _beta +(apb+2)*x[i])*0.5;
            der[0][i] = 0.5*(apb+2);
        }

        
        for(unsigned int k = 1; k < n; ++k ){
            //some ausiliary values
            unsigned int k2 = k*2; 
            unsigned int k2ab=k2+_alpha+_beta;

            pol[2] = pol[1];
            pol[1] = pol[0];

            der[2] = der[1];
            der[1] = der[0];

            unsigned int a1 = 2*(k+1)*(k+1+apb)*k2ab; 
            //Gamma(x+n)/Gamma(x) = (x+n-1) * ... * (x+1) * x 
            unsigned int a21 = (k2ab+1)*ab2;
            unsigned int a22 = (k2ab+2)*(k2ab+1)*k2ab; 
            unsigned int a3 =2*(k+_alpha)*(k+_beta)*(k2ab+2);
            for(unsigned int i = 0; i <x.size(); ++i){
                pol[0][i] = ((a21+a22*x[i])*pol[1][i]-a3*pol[2][i])/a1;
                der[0][i]= (a22*pol[1][i]+(a21+a22*x[i])*der[1][i]-a3*der[2][i])/a1;
            } 

            
        }
        return;
        
    }




};


//template void Quadrature::jacobi_roots<std::iterator>(const unsigned int &n, const double &_alpha, const double &_beta, std::iterator start, std::iterator end);


// void Quadrature::jacobi_roots(const unsigned int &n, const double &_alpha, const double &_beta, std::vector<double> roots)
// {
//     if( n < 1)
//     {
//         std::cout << "The polynomial degree should be greater than 0" << std::endl;
//         return;
//     }
//     else
//     {
//         roots.resize(n);
//         std::vector<double> x0 (1,std::cos(M_PI/(2*n)));
//         double x1(0.);
//         double tol=1.e-14; unsigned int kmax=15;
//         for(unsigned int j = 0; j < n; ++j)
//         {
//             double diff = tol + 1;
//             for(unsigned int kiter = 0; kiter <= kmax || diff >= tol; kiter++)
//             {
//                 std::array<std::vector<double>,3> pol, der;
//                 jacobi_pol(x0, n, _alpha, _beta, pol, der);
//                 // deflation process q(x)=p(x)*(x-x_1)*... (x-x_{j-1})
//                 // q(x)/q'(x)=p(x)/[p'(x)-p(x)*\sum_{i<j} 1/(x-x_i)]
//                 double ss = std::accumulate(roots.begin(),
//                                             roots.begin() + j - 1, 
//                                             0.0,
//                                             [&](double a, double b){return a + 1./ (x0[0]-b);});
//                 x1 = x0[0]-pol[0][0]/(der[0][0]-ss*pol[0][0]);
//                 double diff = std::abs(x1-x0[0]);
//                 x0[0]=x1;
//             }
//             x0[0] = .5 * (x1+std::cos((2*(j+2)-1)*M_PI/(2*n)));
//             roots[j]=x1;
//         }
//         std::sort(roots.begin(), roots.end());
//     }
// };


void Quadrature::legendre_pol(const std::vector<double> &x, const unsigned int &n, std::vector<double> &pol)
    {
        pol.resize(x.size());
        
        if(n == 0)
        {
            std::iota(pol.begin(), pol.end(), 1.);
            return;
        }
        else
        {
            std::vector<double> p1(pol.size(), 1.);
            std::vector<double> p2(x);
            std::vector<double> p3 = p2;

        
            for( unsigned int k = 1; k < n; ++k)
            {
                for(unsigned int i = 0; i <pol.size(); ++i)
                    p3[i] = ((2*k+1)*x[i]*p2[i]-k*p1[i])/(k+1);

                p1 = p2;
                p2 = p3;
            }

            pol = p3;
            
        }

        

    }


void Quadrature::LGL_quadratures(const unsigned int &n/*number of nodes*/)
{

    if(n == 1)
    {
        std::cout <<"The number of quadrature nodes should be greater than 1" << std::endl;
        return;
    }
    else
    {

        _nodes.resize(n);
        _weights.resize(n);

        _nodes[0]= -1.;
        _nodes[n-1]= 1.;
        jacobi_roots(n-2,1.,1., _nodes.begin()+1/*from the second element*/, _nodes.end()-1/*to the last one*/);
        double coeff = 2./(n*(n-1));
        legendre_pol(_nodes, n-1, _weights);

        

        for(unsigned int i = 0; i <_nodes.size(); ++i){
            _weights[i] = coeff/(_weights[i]*_weights[i]);
        }

        //____________________________

        std::cout <<"LGL nodes are:" << std::endl;
        for(auto i : _nodes)
        {
            std::cout << " " << i << " ";
        }
        std::cout << std::endl;
        std::cout <<"LGL weights are:" << std::endl;
        for(auto i : _weights)
        {
            std::cout << " " << i << " ";
        }
        std::cout << std::endl;

        //____________________________

    }
};


void Quadrature::LGL_quadratures(const unsigned int &n/*number of nodes*/,const double &a, const double &b)
{

    if(n == 1)
    {
        std::cout <<"The number of quadrature nodes should be greater than 1" << std::endl;
        return;
    }
    else
    {
        _nodes.resize(n);
        _weights.resize(n);

        double bma = (b - a)/2;
        double bpa = (b + a)/2;

        _nodes[0]= -1.;
        _nodes[n-1]= 1.;
        jacobi_roots(n-2,1,1, _nodes.begin()+1/*from the second element*/, _nodes.end()-1/*to the last one*/);
        double coeff = 2./n*(n-1);
        legendre_pol(_nodes, n-1, _weights);
        for(unsigned int i = 0; i <_nodes.size(); ++i){
            _weights[i] = (coeff/(_weights[i]*_weights[i]))*bma;
            _nodes[i] = bma*_nodes[i]+bpa;
        }

    }

    
};


std::vector<double> Quadrature::getN() const{return _nodes;};
std::vector<double> Quadrature::getW() const{return _weights;};

//______________________________________________________________



//______________________SPECTRAL FE CLASS MEMBERS______________________


Rectangle SpectralFE::getCurrent() const
{
    return _current_elem;
};

Quadrature SpectralFE::getQuad() const
{
    return _qr;
};

unsigned int SpectralFE::getDeg() const
{
    return _r;
};

std::vector<Point> SpectralFE::getQPoints() const
{
    return _quad_points;
};

std::vector<Point> SpectralFE::getINodes() const
{
    return _int_nodes;
};

std::vector<double> SpectralFE::_comp_quad_c()
{
    _qr.LGL_quadratures(this->getDeg()+1);
    return _qr.getN();
};


std::vector<double> SpectralFE::_comp_quad_w()
{
    _qr.LGL_quadratures(this->getDeg()+1);
    return _qr.getW();
};

void SpectralFE::_update_quad()
{
    // delete the quadrature coordinates of the past element
    // OSS. in this case, the degree of the FE spaces along the two dimensions
    // are supposed to be the same and identical to _r member
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
    return;
};

void SpectralFE::_update_nodes()
{
    //delete the nodes of the past element
    _int_nodes.clear();

    _current_elem.getNodes()[0].printNode();
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
            std::tie(std::ignore, s) = _current_elem.inverse_map(j[0],j[1]);
            Point my_point(r,s);
            _int_nodes.emplace_back(my_point);
        }
    }
    return;
};


unsigned int SpectralFE::_dof_per_cell()
{
    return _current_elem.internalN(this->getDeg())*_current_elem.internalN(this->getDeg());

};


Matrix2d SpectralFE::J()
{
    return _current_elem.jacobian();
};

double SpectralFE::detJ()
{
    return _current_elem.jacobian().determinant();
}

Matrix2d SpectralFE::J_inv()
{
    return _current_elem.jacobian().inverse();
}

Matrix2d SpectralFE::J_inv_t()
{
    return _current_elem.jacobian().inverse().transpose();
}


// double shape_value(  /*quadrature point coordinates*/)
// {
//     std::function <double(const double&, const double &)> shape = phi()
// };

// Matrix2d shape_grad( /*quadrature point coordinates*/)
// {

// };


void SpectralFE::_update_current( const Rectangle geoele )
{
    _current_elem = geoele;
    _update_quad();
    _update_nodes();
    return;
};



std::tuple <std::function <double(const double &)>,std::function <double(const double &)>> SpectralFE::legendre_pol(const unsigned int &n)
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


std::tuple<std::function <double(const double &)> ,std::function <double(const double &)>, std::function <double(const double &)>> SpectralFE::legendre_der(const unsigned int &n)
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


std::function <double(const double &)> SpectralFE::lagrange_pol(const double &xi,const unsigned int &n)
{
    //first, we get the legendre polynomials of degree n-1
    auto [ln,lnd] = SpectralFE::legendre_pol(n-1);
    //then generate the function to evaluate the Lagrangian basis for a given coordinate
    std::function<double(const double &)> l_i = [&xi, &n, lnd, ln](const double &x){return (x != xi)?(-1./((n-1)*(n))*(1. - x*x)*lnd(x)/((x - xi)*ln(xi))):1.;};
    return l_i;
};


std::function <double(const double &)> SpectralFE::lagrange_der(const double &xi,const unsigned int &n)
{
    //first, we get the legendre polynomials of degree n-1 and derivative of first order and second order
    auto [ln,lnd, lndd] = SpectralFE::legendre_der(n-1);
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


std::function <double(const double &, const double &)> SpectralFE::phi(const double &xi, const double &xj,const unsigned int &n)
{
    //lagrangian basi function for the two dimension
    std::function <double(const double &)> l_i = SpectralFE::lagrange_pol(xi,n);
    std::function <double(const double &)> l_j = SpectralFE::lagrange_pol(xj,n);
    //the basis funcion for the 2D case is given as tensor product of the two lagrangian polynomials
    std::function <double(const double&, const double &)> phi = [l_i, l_j](const double &x, const double &y){ return l_i(x)*l_j(y);};

    return phi;
}

std::function <double(const double &, const double &)> SpectralFE::dphi(const double &xi, const double &xj,const unsigned int &n)
{
    //lagrangian basi function for the two dimension
    std::function <double(const double &)> l_i = SpectralFE::lagrange_pol(xi,n);
    std::function <double(const double &)> dl_j = SpectralFE::lagrange_der(xj,n);
    //the basis funcion for the 2D case is given as tensor product of the two lagrangian polynomials
    std::function <double(const double&, const double &)> phi = [l_i, dl_j](const double &x, const double &y){ return l_i(x)*dl_j(y);};

    return phi;
}


double SpectralFE::shape_value(const unsigned int &i, Point q)
{

    //first we get the node of the element represented by the index i, using the internal standard node numbering:
    // we do so by looking into the _int_nodes member, containing the coordiates of all the internal nodes of the current element
    double xi,yi;
    std::tie(xi,yi) = (this->getINodes())[i].getCoord();


    // we then get the funtion to evaluate the basis function in a certain point (quadrature points)
    auto phi_i =SpectralFE::phi(xi,yi,this->getDeg()+1);

    //now we evaluate the basis function over the quadrature point

    double x,y;
    
    std::tie(x,y) = q.getCoord();

    return phi_i(x,y);
}

Vector2d SpectralFE::shape_grad(const unsigned int &i, Point q)
{
    //first we get the node of the element represented by the index i, using the internal standard node numbering:
    // we do so by looking into the _int_nodes member, containing the coordiates of all the internal nodes of the current element
    double xi,yi;
    std::tie(xi,yi) = (this->getINodes())[i].getCoord();

    // we then get the funtion to evaluate the derivative of the basis function in a certain point (quadrature points)
    auto dphi_i =SpectralFE::dphi(xi,yi,this->getDeg()+1);

    // now we evaluate the partial derivatives of the basis function over the quadrature point
    double x,y;
    std::tie(x,y) = q.getCoord();
    Vector2d grad(2);
    grad << dphi_i(x,y), dphi_i(y,x);
    return grad;
};






//_______________________________________________________________