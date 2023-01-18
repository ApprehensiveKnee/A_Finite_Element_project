#include "mat_utilities.hpp"




//_________________________ QUADRATURE CLASS MEMBERS ____________________________


void FETools::Quadrature::jacobi_pol(const std::vector<double> &x, const unsigned int &n, const double &_alpha, const double &_beta, std::array<std::vector<double>,3>& pol, std::array<std::vector<double>,3>& der)
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

        for(unsigned int i = 0; i<x.size(); ++i){
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



void FETools::Quadrature::legendre_pol(const std::vector<double> &x, const unsigned int &n, std::vector<double> &pol)
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
                {
                    p3[i] = ((2*k+1)*x[i]*p2[i]-k*p1[i])/(k+1);

                }
                    
                p1 = p2;
                p2 = p3;
            }

            pol = p3;
            
        }

        

    }


void FETools::Quadrature::LGL_quadratures(const unsigned int &n/*number of nodes*/)
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


        //   - PRINTING SECTION -
        //____________________________

        // std::cout <<"LGL nodes are:" << std::endl;
        // for(auto i : _nodes)
        // {
        //     std::cout << " " << i << " ";
        // }
        // std::cout << std::endl;
        // std::cout <<"LGL weights are:" << std::endl;
        // for(auto i : _weights)
        // {
        //     std::cout << " " << i << " ";
        // }
        // std::cout << std::endl;

        //____________________________

    }
};


void FETools::Quadrature::LGL_quadratures(const unsigned int &n/*number of nodes*/,const double &a, const double &b)
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


std::vector<double> FETools::Quadrature::getN() const{return _nodes;};
std::vector<double> FETools::Quadrature::getW() const{return _weights;};

//______________________________________________________________



//______________________SPECTRAL FE CLASS MEMBERS______________________


// as the class has been defined as template class, all the definitions
// of the member functions are grouped inside the corresponding header file


//_______________________________________________________________