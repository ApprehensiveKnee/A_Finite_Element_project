#ifndef ELEM
#define ELEM
 
#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <iterator>
#include <fstream>
#include <sstream>

// Consider 2D problems
inline constexpr unsigned DIM=2;


class Node
{
private:
    //global id of the node
    unsigned int _id;
    //coordinates of the node
    std::array<double, DIM > _coord;
public:
    Node(){};
    //constructor given the coordinates, dimensions and id
    Node(const unsigned int &id, const double &x, const double &y)
        :_id(id),
        _coord ({x,y})
        {};
    //standard id getter
    unsigned int getId() const;
    //strandard coordinates getter
    std::array<double, DIM > getCoord() const;
    //standar setters
    void setX(const double &r);
    void setY(const double &s);
    //standard printer
    void printNode(std::ostream &out = std::cout) const;
    //standard destructor
    ~Node() = default;
};



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
public:
    //constructor given the global indexes of the vertexes
    Rectangle(const unsigned int &elem,const unsigned int &vert1, const unsigned int &vert2, const unsigned int &vert3, const unsigned int &vert4)
        :Element_2D(elem, vert1, vert2, vert3, vert4 )
        {}
    
    Rectangle(const unsigned int &elem, std::vector<unsigned int> &vertexes)
        :Element_2D(elem, vertexes)
        {}
    
    ~Rectangle() = default;
};





//----------------------------

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
    Mesh_2D(bool typ = true):_triatrue(typ){};
    //constructor give the number of verteces and Elements (passed in code)
    Mesh_2D(const unsigned int &V, const unsigned int &E, bool typ = true )
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
};



#endif