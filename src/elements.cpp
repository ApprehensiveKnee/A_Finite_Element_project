#include "elements.hpp"


//_____________________ POINT CLASS MEMEBERS___________________

const std::array<double,DIM>& Point::getCoord() const{return _coord;};

const unsigned short& Point::getBound() const {return _boundary;}

void Point::setBound(const unsigned short& boundary_flag){_boundary = boundary_flag;return;};

const double& Point::getX() const{return _coord[0];};

const double& Point::getY() const{return _coord[1];};

void Point::printPoint(std::ostream &out) const
{
    auto temp = this->getCoord();
    out << temp[0] << ", " <<temp[1] << std::endl;
    out << "boundary: " << this->getBound() << std::endl;
    
    return;
};

std::array<double, DIM> Point::classInit(const double &x, const double &y)
{
    std::array<double, DIM> temp;
    if constexpr (DIM == 1)
    {
        temp[0] = x;
        return temp;
    }
    else if constexpr (DIM == 2)
    {
        temp[0] = x;
        temp[1] = y;
        return temp;
    }
};

//______________________________________________________________


// ______________ DEFINITIONS FOR NODE MEMBER FUNC. ______________


const unsigned int& Node::getId() const{return _id;};

void Node::setX(const double &r){_coord[0]=r; return;};

void Node::setY(const double &s){_coord[1]=s; return;};

void Node::printNode(std::ostream &out) const
{
    out << "==============================" << std::endl;
    out << this->getId() << ": ";
    this->printPoint(out);
    out << "==============================" << std::endl;
    out<<std::endl;
    return;
};


std::vector<std::array<double,DIM>> Node::nodes(const unsigned int &n, const Node &v) const
{
    std::vector<std::array<double,DIM>> temp(Element_1D::internalN(n)); // vector containing the coordinates of the nodes along a certain dimension of the element
    std::array<double,DIM> v1= this->getCoord()/*get the coordinates of vert1*/;
    std::array<double,DIM> v2= v.getCoord()/*get the coordinates of vert2*/;

    if constexpr (DIM == 1)
    {

        double hx = std::abs(v1[0] - v2[0])/n;// step over x
        for(unsigned int i = 0; i< n; i ++) /*in this case, the last element v2 is manually inserted as to avoid round-off error discrepancies*/
        {
        temp[i][0] = v1[0] + hx*i;
        }
        temp[n][0]=v2[0];
        return temp;

    }
    else if constexpr(DIM == 2)
    {
        double hx = std::abs(v1[0] - v2[0])/n;// step over x
        double hy = std::abs(v1[1] - v2[1])/n;// step over y

        for(unsigned int i = 0; i<n ; ++i) /*in this case, the last element v2 is manually inserted as to avoid round-off error discrepancies*/
        {
        temp[i][0]=v1[0] + hx*i;
        temp[i][1]=v1[1] + hy*i;
        }

        temp[n][0]=v2[0];
        temp[n][1]=v2[1];
        return temp;
        
    }
    else
    {
        std::cout << "DIM is out of range"<< std::endl;
        return temp;
    }
}

std::vector<Point> Node::points(const unsigned int &n, const Node &v) const
{
    std::vector<Point> temp(Element_1D::internalN(n)); // vector containing the coordinates of the nodes along a certain dimension of the element
    std::array<double,DIM> v1= this->getCoord()/*get the coordinates of vert1*/;
    unsigned short b1 = this->getBound()/*get thr boundary of vert1*/;
    std::array<double,DIM> v2= v.getCoord()/*get the coordinates of vert2*/;
    unsigned short b2 = v.getBound()/*get the boundary of vert2*/;

    if constexpr (DIM == 1)
    {

        double hx = std::abs(v1[0] - v2[0])/n;// step over x

        for(unsigned int i = 1; i< n; i ++) /*in this case, the last element v2 is manually inserted as to avoid round-off error discrepancies*/
        {
            temp[i] =Point(v1[0] + hx*i, 0,(i == 0)? b1:0) ;
        }
        temp[n]=Point(v2, b2);
        return temp;

    }
    else if constexpr(DIM == 2)
    {
        double hx = std::abs(v1[0] - v2[0])/n;// step over x
        double hy = std::abs(v1[1] - v2[1])/n;// step over y

        for(unsigned int i = 0; i<n ; ++i) /*in this case, the last element v2 is manually inserted as to avoid round-off error discrepancies*/
        {
            unsigned short b(0);
            if(b1 == b2)
            {
                b = b1;
            }
            // ________ SPECIAL CASES FOR FIRST VERTEX ________

            else if(b1==31-1)// start point == vertex 1
            {
                //either used to compute points on side 1 or 4
                if(b2 == 31-2) // last == vertex 2
                    b = 1; //side 1;
                else if(b2 == 31-4) // last == vertex 4
                    b = 4; // side 4;
                else 
                    b = b2;
            }
            else if(b1 == 31-2) // start point == vertex 2
            {
                //used to compute points on side 2
                b = 2;
            }
            else if(b1 == 31 -4) // start point == vertex 4
            {
                //used to compute points on side 3
                b = 3;
            }
            // ________ SPECIAL CASES FOR LAST VERTEX ________

            else if(b2 == 31 - 3) // last point == vertex 3
            {
                //used to compute either points on side 2 or 3
                if(b1 == 31 -4) // first == vertex 4
                    b = 3; //side 3
                else if( b1 == 31 - 3) // first == vertex 2
                    b = 2; //side 2
                else 
                    b = b1;
            }
            else if(b2 == 31 -2) // last point == vertex 2
            {
                //used to compute points on side 1
                b = 1;
            }
            else if(b2 == 31 - 4)// last point == vertex 4
            {
                //used to compute points on side 4
                b = 4;
            }
            else
            {
                //otherwise b is an internal point
                b = 0;
            }
            
            temp[i] = Point(v1[0] + hx*i, v1[1] + hy*i, (i==0)?b1:b);
        }

        temp[n]=Point(v2, b2);
        return temp;
        
    }
    else
    {
        std::cout << "DIM is out of range"<< std::endl;
        return temp;
    }
}


//______________________________________________________________


// ______________ DEFINITIONS FOR MESH ELEMENT MEMBER FUNC. ______________

const unsigned int& MeshElement::getId() const{return _element_id;};

const unsigned short& MeshElement::getBound() const{return _boundary;};

const std::array<unsigned int, DIM>& MeshElement::getNQ() const{return _nq;};

void MeshElement::setBound(const unsigned short& boundary_flag){_boundary = boundary_flag; return;}

unsigned int MeshElement::internalN(const unsigned int &nx){return nx + 1;};

std::array<unsigned int,DIM> MeshElement::classInit(const unsigned int &nx, const unsigned int &ny)
{
    std::array<unsigned int, DIM> temp;
    if constexpr (DIM == 1)
    {
        temp[0] = nx;
        return temp;
    }
    else if constexpr (DIM == 2)
    {
        temp[0] = nx;
        temp[1] = ny;
        return temp;
    }
};

//______________________________________________________________


// ______________ DEFINITIONS FOR 1D ELEMENT MEMBER FUNC. ______________



const std::array<Node,2>& Element_1D::getNodes() const{return _nodes;};

void Element_1D::setNodes(const std::array<Node,2> &nodes){_nodes = nodes; return;};

void Element_1D::printElem(std::ostream &out) const
{
    //print the global index of the element
    out<<"Vertices of Element #"<<this->getId()<<" are:" << std::endl;
    out << "----------"<<std::endl;
    for (auto i : _nodes)
    {
        //print the ids of the nodes of the element
        out <<" "<< i.getId() <<" ";
        
    }
    
    out <<std::endl << "----------"<<std::endl;
};


double Element_1D::jacobian() const
{
    double J = (this->getNodes()[1].getX() - this->getNodes()[0].getX())/2.;
    return J;
};


double Element_1D::inverse_map(const double &x) const
{
    
    double J = (this->jacobian());
    double x_ref = x/J - _nodes[0].getX();
    return x_ref -1;
};

double Element_1D::direct_map(const double &x_ref) const
{
    double J = (this->jacobian());
    double x = J*(x_ref+1.)+ this->getNodes()[0].getX();
    return x;
};



//______________________________________________________________


// ______________ DEFINITIONS FOR 2D ELEMENT MEMBER FUNC. ______________


const std::vector<Node>& Element_2D::getNodes() const{return _nodes;};


void Element_2D::setNodes(const std::vector<Node> &nodes)
{
    _nodes = nodes; 
    return;
};


void Element_2D::printElem(std::ostream &out) const
{
    //print the global index of the element
    out<<"Vertices of Element #"<<this->getId()<<" are:" << std::endl;
    out << "----------"<<std::endl;
    for (auto i : _nodes)
    {
        //print the ids of the nodes of the element
        out <<" "<< i.getId() <<" ";
        
    }
    
    out <<std::endl << "----------"<<std::endl;
};


Matrix2d Element_2D::jacobian() const
{
    //x_k = B*x_ref + c where c = v1_k and B is the jacobian J = [v2_k - v1_x , v3_k - v1_k], x = [ x1_k , x2_k]^T, x_ref = [x1_ref, x2_ref]^T

    Matrix2d J(2,2);
    J << (_nodes[1].getX() - _nodes[0].getX())/2.,
        (_nodes[1].getY() - _nodes[0].getY())/2.,
        (_nodes[2].getX() - _nodes[0].getX())/2.,
        (_nodes[2].getY() - _nodes[0].getY())/2.;
    return J;
};


std::tuple<double,double> Element_2D::inverse_map(const double &x, const double &y) const
{
    Matrix2d J(this->jacobian());
    Vector2d c(2), x_ref(2), x_k(2);
    c << _nodes[0].getX(), _nodes[0].getY();
    x_k << x,y;
    x_ref = J.inverse()*(x_k-c);
    return { x_ref(0) - 1.,x_ref(1) - 1.}; 
    //added -1  and multiplied by 2 ,otherwise the inverse transformation would bring us to
    // the refernce element [0,1]x[0,1]

};

std::tuple<double,double> Element_2D::direct_map(const double &x_r, const double &y_r) const
{
    Matrix2d J(this->jacobian());
    Vector2d c(2), x_ref(2), x_k(2);
    c << _nodes[0].getX(), _nodes[0].getY();
    x_ref << x_r+1.,y_r+1.;
    x_k = J*x_ref + c;
    return { x_k(0),x_k(1)}; 
    //added -1  and multiplied by 2 ,otherwise the inverse transformation would bring us to
    // the refernce element [0,1]x[0,1]

};

//______________________________________________________________


