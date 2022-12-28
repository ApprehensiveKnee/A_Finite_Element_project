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

// Consider 2D problems
inline constexpr unsigned DIM=1;


class Node
{
private:
    //global id of the node
    unsigned int _id;
    //coordinates of the node
    std::vector<double> _coord;
public:
    Node(){};
    //constructor given the coordinates and id (1D case)
    Node(const unsigned int &id, const double &x)
        :_id(id),
        _coord({x})
        {};
    //constructor given the coordinates and id (2D case)
    Node(const unsigned int &id, const double &x, const double &y)
        :_id(id),
        _coord({x,y})
        {};
    //standard id getter
    unsigned int getId() const;
    //strandard coordinates getter
    std::vector<double> getCoord() const;
    //standar setters
    void setX(const double &r);
    void setY(const double &s);
    //standard printer
    void printNode(std::ostream &out = std::cout) const;
    //some methods to return a vector containing the coordinates of all the internal nodes 
    //along a certain dimension of the element based off the degree of the polinomial used for such dimension

    //the method for Segment(1D)
    // in this case, nx represents the degree of the polynomial used while v is the other
    //extremis of the segment for which we have to determine the coordinates of the internal nodes
    std::vector<std::vector<double>> nodes1D(const unsigned int &nx, const Node &v);
    //the methods for Element_2D(2D)
    //in this case the method can be used either to compute the coordinates of the internal nodes over x
    //or y directions: n is the degree of the polynomial chosen over a certain dimension, while v represent
    // the other extremis of the side for which we want to determine the internal nodes
    std::vector<std::vector<double>> nodes2D(const unsigned int &n/*either nx or ny*/, const Node &v/*either v1 or v3-v4 (triangles or rectangles)*/);


    //standard destructor
    ~Node() = default;
};


//Mesh element, 1D case

class Segment
{
private:
    //id of the element
    unsigned int _element_id;
    // global indexes of the nodes that compose the element
    std::vector<unsigned int> _vert;
public:
    //default constructor
    Segment() = default;
    //constructor given the index of vertexes
    Segment(const unsigned int &id, const unsigned int &vert1, const unsigned int &vert2)
        :_element_id(id),
        _vert({vert1, vert2})
        {}
    //constructor givenn a vector of vertexes
    Segment(const unsigned int &id, const std::vector<unsigned int> vertexes)
        :_element_id(id),
        _vert(vertexes)
        {}
    //getter for the element id
    unsigned int getId() const;
    //getter for the ids of the vertexes
    std::vector<unsigned int> getVert() const;
    // printer: prints the id and vertexes of the element
    void printElem(std::ostream &out = std::cout) const;
    // a method to compute the number of internal nodes for each element, given the degree of the polunamilas used
    // in each dimension for the spectral elements
    static unsigned int internalN(const unsigned int &nx);
    
    //default destructor
    ~Segment() = default;

};

// Base class, 2D element

class Element_2D
{
protected:
    //id of the element
    unsigned int _element_id;
    // global indexes of the nodes that compose the element
    std::vector<unsigned int> _vert;

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
    //constructor trough a vector
    Element_2D(const unsigned int &id, const std::vector<unsigned int> vertexes)
        :_element_id(id),
        _vert(vertexes)
        {}
    //getter:gets the id of the element
    unsigned int getId() const;
    //getter: gets the ids of the vertexes
    std::vector<unsigned int> getVert() const;
    //printer: prints the ids and coordinates of the vertexes
    void printElem(std::ostream &out = std::cout) const;
    // a method to compute the number of internal nodes for each element, given the degree of the polunamilas used
    // in each dimension for the spectral elements
    static unsigned int internalN(const unsigned int &nx);
    
    
    //virtual default destructor
    virtual ~Element_2D(){};
};



class Triangle :public Element_2D
{
public:
    //constructor given the global indexes of the vertexes
    Triangle(const unsigned int &elem, const unsigned int &vert1, const unsigned int &vert2, const unsigned int &vert3)
        :Element_2D(elem, vert1, vert2, vert3)
        {}
    
    Triangle(const unsigned int &elem, std::vector<unsigned int> &vertexes)
        :Element_2D(elem, vertexes)
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
//            V4 __________ V3
//              |          |
//              |          |
//      side 4  | Omega    |  side 2
//              |          |
//              |__________|
//            V1            V2
//                 side 1
public:
    //constructor given the global indexes of the vertexes
    Rectangle(const unsigned int &elem,const unsigned int &vert1, const unsigned int &vert2, const unsigned int &vert3, const unsigned int &vert4)
        :Element_2D(elem, vert1, vert2, vert3, vert4 )
        {}
    
    Rectangle(const unsigned int &elem, std::vector<unsigned int> &vertexes)
        :Element_2D(elem, vertexes)
        {}
    //default destructor overridden
    ~Rectangle() = default;
};



// Mesh 1D case


class Mesh_1D
{
private:
    //number of total nodes of the mesh
    unsigned int _nNodes;
    //number of total elements of the mesh
    unsigned int _nElems;
    //vector containing the nodes and elements of the mesh
    std::vector<Node> _nodes;
    std::vector<Segment> _segments;
public:
    //standard constructor
    Mesh_1D(){};
    //constructor given the number of vertexes and elements
    Mesh_1D(const unsigned int &V, const unsigned int &E)
        :_nNodes(V),
        _nElems(E)
        {}
    //getter for the number of Nodes
    unsigned int get_nNodes() const;
    //getter fotr the number of Elements
    unsigned int get_nElems() const;
    //standard printer
    void printMesh(std::ostream& out = std::cout) const;
    //setter for the number of elements and nodes for a txt file
    void setNum();
    //setter for the nodes and elements of the mesh
    //this method reads the number of elements and nodes from a txt file,
    //the nodes coordinates composing the mesh form a csv and finally the elements (to determine how to nodes are connected)
    //for another csv file    
    void setMesh();
    //member function to map the global indexes of the nodes in the messh to the local indexes of the nodes inside in the elements
    //  THE GLOBAL INDEX OF THEE i LOCAL NODE OF A SPECIFIC ie ELEMENT CAN BE COMPUTED BY ACCESSING TO THE (i-1,ie-1) ELEMENT OF THE 
    //  VECTOR (MULTIDIMENSIONAL) RETUNED BY THE FUNCTION
    static std::vector<std::vector<unsigned int>> indexMapping(const unsigned int &degree/*degree of the polynomial*/, const unsigned int &nElements/*number of elements*/);
    //default destructor
    ~Mesh_1D() = default;



};

//Mesh 2D case

class Mesh_2D
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
    std::vector<Triangle> _triangles;
    std::vector<Rectangle> _rectangles;
public:
    //standard constructor
    Mesh_2D(const bool &typ = true):_triatrue(typ){};
    //constructor given the number of vertexes and Elements (passed as arguments in main)
    Mesh_2D(const unsigned int &V, const unsigned int &E, const bool &typ = true )
        :_nNodes(V),
        _nElems(E),
        _triatrue(typ)
        {};
    //getter for the number of Nodes
    unsigned int get_nNodes() const;
    //getter for the number of Elements
    unsigned int get_nElems() const;
    //standard printer
    void printMesh(std::ostream& out = std::cout) const;
    //setter for the number of elements and nodes for a txt file
    void setNum();
    //setter for the nodes and elements of the mesh
    //this method reads the number of elements and nodes from a txt file,
    //the nodes coordinates composing the mesh form a csv and finally the elements (to determine how to nodes are connected)
    //for another csv file    
    void setMesh();
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
                                                const unsigned int &degreeY/*degree of the polynomial along y axis*/,
                                                const unsigned int &nElementsY/*number of elements along y axis*/);
    ~Mesh_2D() = default;
};



#endif