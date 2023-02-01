//===========================================================
// HEADER FILE FOR THE CLASS NODE AND ELEMENTS DEFINITIONS
//===========================================================

#ifndef ELEM
#define ELEM
 
#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <numeric>
#include <iterator>
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

// Define some constant quantities...

inline constexpr unsigned DIM=2;
inline constexpr double tol= 1.e-9;
inline constexpr unsigned r = 2;

using namespace Eigen;

//______________________________________________________________


//a class to store the coordinates of the relevant points (quadrature coordinates)
class Point
{
protected:
    // a member array to store the coordinates
    std::array<double, DIM> _coord;
    // a member to store boundary information:
    // if considering the mesh read form csv of the one generated inside the code, the convetion used is the following:
    // 
    //                 side 3
    //            V4 __________ V3
    //              |          |
    //              |          |
    //      side 4  |          |  side 2    // Consider the following mesh
    //              |          |
    //              |__________|
    //            V1            V2
    //                 side 1
    //
    // 0 -> internal node
    // 1 -> node on the bottom border
    // 2 -> node on the right border
    // 3 -> node on the top border
    // 4 -> node on the left border
    // The following numbering is assumed to better handle the conditions on the vertexes
    // 31 -1 -> for V1
    // 31 -2 -> for V2
    // 31 -3 -> for V3
    // 31- 4 -> for V4
    // In 1D the first vertex is simply indicated with 1, the last with 2
    unsigned short _boundary;
public:
    Point(){};
    Point(const double &x, const double &y = 0., const unsigned short& boundary = 0):
        _coord(classInit(x,y)),
        _boundary(boundary)
        {};
    Point(const std::array<double, DIM> coord, const unsigned short& boundary = 0):
        _coord(coord),
        _boundary(boundary)
        {}
    //standard getters
    const std::array<double,DIM>& getCoord() const;
    const unsigned short& getBound() const;
    void setBound(const unsigned short& boundary_flag);
    const double& getX() const;
    const double& getY() const;
    //a function to print the point
    void printPoint(std::ostream &out = std::cout) const;
    //a static function to initialize the class
    static std::array<double,DIM> classInit(const double &x, const double &y);
    
    ~Point() = default;

};

//______________________________________________________________

class Node : public Point
{
private:
    //global id of the node
    unsigned int _id;
public:
    Node(){};
    //constructor given the coordinates and id (2D case)
    Node(const unsigned int &id, const double &x, const double &y = 0., const unsigned short& boundary = 0)
        :Point(x,y,boundary),
        _id(id)
        {};
    //standard id getter
    const unsigned int& getId() const;
    //standar setters
    void setX(const double &r);
    void setY(const double &s);
    //standard printer
    void printNode(std::ostream &out = std::cout) const;
    //a method to return a vector containing the coordinates of all the internal nodes 
    //along a certain dimension of the element based off the degree of the polinomial used for such dimension
    std::vector<std::array<double,DIM>> nodes(const unsigned int &n, const Node &v) const;

    //a correponding method returning a vector of points
    std::vector<Point> points(const unsigned int &n, const Node &v) const;

    //standard destructor
    ~Node() = default;
};


class MeshElement
{
protected:
    //id of the element
    unsigned int _element_id;
    // a member to store the boundary flag
    //  the convetion adopted is the same as the Point boundary flag
    unsigned short _boundary;
    // members to store the number of quadrature points to compute the integrals on the element,both for the x and y dimensions
    std::array<unsigned int, DIM> _nq;
public:
    MeshElement(){};
    MeshElement(const unsigned int &id, const unsigned short& boundary, const unsigned int& nx, const unsigned int& ny = 0)
        :_element_id(id),
        _boundary(boundary),
        _nq(classInit(nx,ny))
        {}
    //getter:gets the id of the element
    const unsigned int& getId() const;
    //getter: gets the boundary flag, to implement a tree search of the boundary nodes
    const unsigned short& getBound() const;
    //getter: gets the degree of exactness
    const std::array<unsigned int, DIM>& getNQ() const;
    //setter: sets the boundary flag after the initialization of the element
    void setBound(const unsigned short& boundary_flag);
    // printer: prints the id and vertexes of the element
    //virtual void printElem(std::ostream &out = std::cout) const;
    // a method to compute the number of internal nodes for each element, given the degree of the polunamilas used
    // in each dimension for the spectral elements
    static unsigned int internalN(const unsigned int &nx);
    // a method for class initializiation
    static std::array<unsigned int,DIM> classInit(const unsigned int &x, const unsigned int &y);

    virtual ~MeshElement() = default;
};

//______________________________________________________________


//Mesh element, 1D case

class Element_1D : public MeshElement
{
private:
    
    // vector of Node elements corresponding to vert indexes
    std::array<Node,2> _nodes;
public:
    //default constructor
    Element_1D()
        :MeshElement()
        {}
    //constructor give the nodes of the element
        Element_1D(const unsigned int &id,
                   const Node &node1 = Node(1, -1., 0), const Node &node2 = Node(2, 1., 0),
                   const unsigned short &boundary = 0,
                   const unsigned int &nqx = r+1)
            : MeshElement(id, boundary, nqx),
              _nodes({node1, node2})
            {}
    //constructor given a vector of Node elements and index element
        Element_1D(const unsigned int &id,
                   const std::array<Node, 2> &nodes = {Node(1, -1., 0), Node(2, 1., 0)},
                   const unsigned short &boundary = 0,
                   const unsigned int &nqx = r+1)
            : MeshElement(id, boundary, nqx),
              _nodes(nodes)
            {}
    //getter:gets the nodes of the element
    const std::array<Node,2>& getNodes() const;
    //setter for the nodes array (to get the coordinates)
    void setNodes(const std::array<Node,2> &nodes);
    // printer: prints the id and vertexes of the element
    void printElem(std::ostream &out = std::cout) const;
    // a method to compute the jacobian (double in 1D)
    double jacobian() const;
    // a member to compute the inverse transformation onto the reference element 
    double inverse_map(const double &x) const;
    double direct_map(const double &x_ref) const;
    
    //default destructor
    ~Element_1D() = default;

};


//______________________________________________________________

// Base class, 2D element

class Element_2D : public MeshElement
{
protected:
    // vector of Node elements corresponding to vert indexes
    std::vector<Node> _nodes;

    // the following is a viual representation of the ordering of vertices adopted for the construction
    // of the rectanglular element to be used for the SEM mesh
    //
    //                
    //            V3 __________ V4
    //              |          |
    //              |          |
    //              |          |  
    //              |          |
    //              |__________|
    //            V1            V2
    //                 

public:
    //default constructor
    Element_2D()
        :MeshElement()
        {}
    //constructor taking the nodes of the element
        Element_2D(const unsigned int &id,
                   const Node &node1 = Node(1, -1., -1.),
                   const Node &node2 = Node(2, 1., -1.),
                   const Node &node3 = Node(3, -1., 1.),
                   const Node &node4 = Node(4, 1., 1.),
                   const unsigned short& boundary = 0,
                   const unsigned int& nqx = r+1,
                   const unsigned int& nqy = r+1)
            : MeshElement(id, boundary,nqx,nqy),
            _nodes({node1, node2, node3, node4})
        {}
    //constructor by Node vector-passing
        Element_2D(const unsigned int &id,
                   const std::vector<Node> &nodes = {
                       Node(1, -1., -1.),
                       Node(2, 1., -1.),
                       Node(3, -1., 1.),
                       Node(4, 1., 1.),
                   },
                   const unsigned short &boundary = 0,
                   const unsigned int& nqx = r+1,
                   const unsigned int& nqy = r+1) 
            : MeshElement(id, boundary,nqx,nqy),
            _nodes(nodes)
        {}
    //getter:gets the nodes of the element
    const std::vector<Node>& getNodes() const;
    //setter for the nodes array (to get the coordinates)
    void setNodes(const std::vector<Node> &nodes);
    //printer: prints the ids and coordinates of the vertexes
    void printElem(std::ostream &out = std::cout) const;
    // a method to compute the jacobian of a specific element, describing the affine tranformation, using eigen classes
    Matrix2d jacobian() const;
    // a method to compute the inverse transformation onto the refernce element [-1,1]x[-1,1]
    std::tuple<double,double> inverse_map(const double &x, const double &y) const;
    // a member to compute the direct transformation from the reference element
    std::tuple<double,double> direct_map(const double &x_ref, const double &y_ref) const;
    
    //default destructor
    ~Element_2D(){};
};





// class Triangle :public Element_2D
// {
// public:
//     Triangle()
//         :Element_2D()
//         {}
//     //constructor given the global indexes of the vertexes
//     Triangle(const unsigned int &elem, const unsigned int &vert1 = 1, const unsigned int &vert2 = 2, const unsigned int &vert3 = 3)
//         :Element_2D(elem, vert1, vert2, vert3)
//         {}
    
//     Triangle(const unsigned int &elem, std::vector<unsigned int> &vertexes)
//         :Element_2D(elem, vertexes)
//         {}

//     Triangle(const unsigned int &elem, std::vector<Node> &nodes)
//         :Element_2D(elem, nodes)
//         {}
//     // virtual destructor overridden
//     ~Triangle() = default;
// };





// class Rectangle : public Element_2D
// {


// public:
//     Rectangle()
//         :Element_2D()
//         {}
//     //constructor given the global indexes of the vertexes
//     Rectangle(const unsigned int &elem,const unsigned int &vert1 = 1, const unsigned int &vert2 = 2, const unsigned int &vert3 = 3, const unsigned int &vert4 = 4)
//         :Element_2D(elem, vert1, vert2, vert3, vert4 )
//         {}
    
//     Rectangle(const unsigned int &elem, std::vector<unsigned int> &vertexes)
//         :Element_2D(elem, vertexes)
//         {}
    
//     Rectangle(const unsigned int &elem, std::vector<Node> &nodes)
//         :Element_2D(elem,  nodes)
//         {}
//     //default destructor overridden
//     ~Rectangle() = default;
// };


//______________________________________________________________



#endif