#ifndef ELEM
#define ELEM
 
#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <iterator>
#include <fstream>
#include <sstream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

// Consider xD problems
inline constexpr unsigned DIM=2;
using namespace Eigen;


//______________________________________________________________

class Node
{
private:
    //global id of the node
    unsigned int _id;
    //coordinates of the node
    std::array<double,DIM> _coord;
    //boundary flag
    unsigned short _boundary;
public:
    Node(){};
    //constructor given the coordinates and id (2D case)
    Node(const unsigned int &id, const double &x, const double &y = 0, const unsigned short &boundary = 0 )
        :_id(id),
        _coord(classInit(x,y)),
        _boundary(boundary)
        {};
    //standard id getter
    unsigned int getId() const;
    //strandard coordinates getter
    std::array<double,DIM> getCoord() const;
    double getX() const;
    double getY() const;
    //standard getter for the boundary flag
    unsigned short getBound() const;
    //standar setters
    void setX(const double &r);
    void setY(const double &s);
    void setBound(const unsigned short &b);
    //standard printer
    void printNode(std::ostream &out = std::cout) const;
    //a method to return a vector containing the coordinates of all the internal nodes 
    //along a certain dimension of the element based off the degree of the polinomial used for such dimension
    std::vector<std::array<double,DIM>> nodes(const unsigned int &nx, const Node &v);

    //a static functio to initialize the class memebers
    static std::array<double,DIM> classInit(const double &x, const double &y = 0);

    //standard destructor
    ~Node() = default;
};


//______________________________________________________________


//Mesh element, 1D case

class Segment
{
private:
    //id of the element
    unsigned int _element_id;
    // global indexes of the nodes that compose the element
    std::array<unsigned int,2> _vert;
    // vector of Node elements corresponding to vert indexes
    std::array<Node,2> _nodes;
public:
    //default constructor
    Segment() = default;
    //constructor given the index of vertexes
    Segment(const unsigned int &id, const unsigned int &vert1, const unsigned int &vert2)
        :_element_id(id),
        _vert({vert1, vert2})
        {}
    //constructor given a vector of vertexes
    Segment(const unsigned int &id, const std::array<unsigned int,2> vertexes)
        :_element_id(id),
        _vert(vertexes)
        {}
    //constructor given a vector of Node elements and index element
    Segment(const unsigned int &id, std::array<Node,2> nodes)
        :_element_id(id),
        _vert({nodes[0].getId(), nodes[1].getId()}),
        _nodes(nodes)
        {}
    //getter for the element id
    unsigned int getId() const;
    //getter for the ids of the vertexes
    std::array<unsigned int,2> getVert() const;
    //getter:gets the nodes of the element
    std::array<Node,2> getNodes() const;
    //setter for the nodes array (to get the coordinates)
    void setNodes(const std::array<Node,2> &nodes);
    // printer: prints the id and vertexes of the element
    void printElem(std::ostream &out = std::cout) const;
    // a method to compute the number of internal nodes for each element, given the degree of the polunamilas used
    // in each dimension for the spectral elements
    static unsigned int internalN(const unsigned int &nx);
    // a method to compute the jacobian (double in 1D)
    double jacobian() const;
    // a member to compute the inverse transformation onto the reference element 
    double inverse_map(const double &x) const;
    
    //default destructor
    ~Segment() = default;

};


//______________________________________________________________

// Base class, 2D element

class Element_2D
{
protected:
    //id of the element
    unsigned int _element_id;
    // global indexes of the nodes that compose the element
    std::vector<unsigned int> _vert;
    // vector of Node elements corresponding to vert indexes
    std::vector<Node> _nodes;

public:
    //default constructor
    Element_2D(){};
    //constructor given the global indexes of the vertexes
    Element_2D(const unsigned int &id, const unsigned int &vert1, const unsigned int &vert2, const unsigned int &vert3)
        :_element_id(id),
        _vert({vert1, vert2, vert3})
        {}
    // //constructor given the global indexes of the vertexes
    Element_2D(const unsigned int &id, const unsigned int &vert1, const unsigned int &vert2, const unsigned int &vert3, const unsigned int &vert4)
        :_element_id(id),
        _vert({vert1, vert2, vert3, vert4})
        {}
    //constructor by index vector-passing
    Element_2D(const unsigned int &id, const std::vector<unsigned int> vertexes)
        :_element_id(id),
        _vert(vertexes)
        {}
    //constructor  by Node vector-passing
    Element_2D(const unsigned int &id,const std::vector<Node> nodes)
        :_element_id(id),
        _nodes(nodes)
        {for(auto i : nodes)
            _vert.emplace_back(i.getId());
        }
    //getter:gets the id of the element
    unsigned int getId() const;
    //getter: gets the ids of the vertexes
    std::vector<unsigned int> getVert() const;
    //getter:gets the nodes of the element
    std::vector<Node> getNodes() const;
    //setter for the nodes array (to get the coordinates)
    void setNodes(const std::vector<Node> &nodes);
    //printer: prints the ids and coordinates of the vertexes
    void printElem(std::ostream &out = std::cout) const;
    // a method to compute the number of internal nodes for each element, given the degree of the polunamilas used
    // in each dimension for the spectral elements
    static unsigned int internalN(const unsigned int &nx);
    // a method to compute the jacobian of a specific element, describing the affine tranformation, using eigen classes
    Matrix2d jacobian() const;
    // a method to compute the inverse transformation onto the refernce element [-1,1]x[-1,1]
    std::tuple<double,double> inverse_map(const double &x, const double &y) const;
    
    //virtual default destructor
    virtual ~Element_2D(){};
};





class Triangle :public Element_2D
{
public:
    Triangle()
        :Element_2D()
        {}
    //constructor given the global indexes of the vertexes
    Triangle(const unsigned int &elem, const unsigned int &vert1, const unsigned int &vert2, const unsigned int &vert3)
        :Element_2D(elem, vert1, vert2, vert3)
        {}
    
    Triangle(const unsigned int &elem, std::vector<unsigned int> &vertexes)
        :Element_2D(elem, vertexes)
        {}

    Triangle(const unsigned int &elem, std::vector<Node> &nodes)
        :Element_2D(elem, nodes)
        {}
    // virtual destructor overridden
    ~Triangle() = default;
};





class Rectangle : public Element_2D
{

// the following is a viual representation of the ordering of vertices adopted for the construction
// of the rectanglular element to be used for the SEM mesh
//
//                 side 3
//            V3 __________ V4
//              |          |
//              |          |
//      side 4  |          |  side 2
//              |          |
//              |__________|
//            V1            V2
//                 side 1
public:
    Rectangle()
        :Element_2D()
        {}
    //constructor given the global indexes of the vertexes
    Rectangle(const unsigned int &elem,const unsigned int &vert1, const unsigned int &vert2, const unsigned int &vert3, const unsigned int &vert4)
        :Element_2D(elem, vert1, vert2, vert3, vert4 )
        {}
    
    Rectangle(const unsigned int &elem, std::vector<unsigned int> &vertexes)
        :Element_2D(elem, vertexes)
        {}
    
    Rectangle(const unsigned int &elem, std::vector<Node> &nodes)
        :Element_2D(elem,  nodes)
        {}
    //default destructor overridden
    ~Rectangle() = default;
};


//______________________________________________________________




//Mesh

class Mesh
{
private:
    //number of total nodes of the mesh
    unsigned int _nNodes;
    //number of total elements of the mesh
    unsigned int _nElems;
    //flag for the type of elements used
    bool _triatrue;
    //vectors containing the nodes and elements of the mesh
    std::vector<Node> _nodes;
    std::vector<Segment> _segments;
    std::vector<Triangle> _triangles;
    std::vector<Rectangle> _rectangles;
public:
    //standard constructor
    Mesh(const bool &typ = true) 
        :_triatrue(typ)
        {}
    //constructor given the number of vertexes and Elements (passed as arguments in main)
    Mesh(const unsigned int &V, const unsigned int &E, const bool &typ = true )
        :_nNodes(V),
        _nElems(E),
        _triatrue(typ)
        {}
    //getter for the number of Nodes
    unsigned int get_nNodes() const;
    //getter for the number of Elements
    unsigned int get_nElems() const;
    // getter for the flag
    bool get_triatrue() const;
    //standard printer
    void printMesh(std::ostream& out = std::cout) const;
    //setter for the number of elements and nodes for a txt file
    void setNum();
    //setter for the nodes and elements of the mesh
    //this method reads the number of elements and nodes from a csv file,
    //the nodes coordinates composing the mesh form a csv and finally the elements (to determine how to nodes are connected)
    //for another csv file    
    void setMesh_csv();
    //member function to map the global indexes of the nodes in the messh to the local indexes of the nodes inside in the elements
    //  THE GLOBAL INDEX OF THEE i LOCAL NODE OF A SPECIFIC ie ELEMENT CAN BE COMPUTED BY ACCESSING TO THE (i-1,ie-1) ELEMENT OF THE 
    //  VECTOR (MULTIDIMENSIONAL) RETUNED BY THE FUNCTION


    //             __________________________
    //            |      |      |     |     |
    //            |  3   |  6   |  9  | 12  |      Omega and spectral elements 
    //            |      |      |     |     |      ordering
    //            __________________________
    //            |      |      |     |     |
    //            |  2   |  5   |  8  | 11  |    
    //            |      |      |     |     |
    //            __________________________
    //            |      |      |     |     |
    //            |  1   |  4   |  7  | 10  |
    //            |      |      |     |     |
    //            __________________________



    static std::vector<std::vector<unsigned int>>indexMapping(const unsigned int &degreeX/*degree of the polynomial along x axis*/, 
                                                const unsigned int &nElementsX/*number of elements along x axis*/,
                                                const unsigned int &degreeY = 0/*degree of the polynomial along y axis*/,
                                                const unsigned int &nElementsY = 0/*number of elements along y axis*/);
    ~Mesh() = default;
};



#endif