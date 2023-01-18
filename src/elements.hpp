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

// Consider xD problems
inline constexpr unsigned DIM=2;
inline constexpr double tol= 1.e-9;
inline constexpr double inf = -5.e2;
using namespace Eigen;

//______________________________________________________________


//a class to store the coordinates of the relevant points (quadrature coordinates)
class Point
{
protected:
    std::array<double, DIM> _coord;
public:
    Point(){};
    Point(const double &x, const double &y = 0):
        _coord(classInit(x,y))
        {};
    //standard getters
    const std::array<double,DIM>& getCoord() const;
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
    Node(const unsigned int &id, const double &x, const double &y = 0)
        :Point(x,y),
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
    std::vector<std::array<double,DIM>> nodes(const unsigned int &nx, const Node &v) const;

    //standard destructor
    ~Node() = default;
};


class MeshElement
{
protected:
    //id of the element
    unsigned int _element_id;
public:
    MeshElement(){};
    MeshElement(const unsigned int &id)
        :_element_id(id)
        {}
    //getter:gets the id of the element
    const unsigned int& getId() const;
    // printer: prints the id and vertexes of the element
    //virtual void printElem(std::ostream &out = std::cout) const;
    // a method to compute the number of internal nodes for each element, given the degree of the polunamilas used
    // in each dimension for the spectral elements
    static unsigned int internalN(const unsigned int &nx);
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
    Element_1D(const unsigned int &id, const Node &node1, const Node &node2)
        :MeshElement(id),
        _nodes({node1,node2})
        {}
    //constructor given a vector of Node elements and index element
    Element_1D(const unsigned int &id, const std::array<Node,2> & nodes)
        :MeshElement(id),
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
    //default constructor
    Element_2D()
        :MeshElement()
        {}
    //constructor taking the nodes of the element
    Element_2D(const unsigned int &id, const Node &node1, const Node &node2, const Node &node3, const Node &node4)
        :MeshElement(id),
        _nodes({node1, node2, node3, node4})
        {}
    //constructor by Node vector-passing
    Element_2D(const unsigned int &id,const std::vector<Node>& nodes)
        :MeshElement(id),
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




//Mesh
template<class ElementType>
class Mesh
{
private:
    //number of total nodes of the mesh
    unsigned int _nNodes;
    //number x,y elements of the mesh
    std::array<unsigned int, DIM> _nElems;
    //vectors containing the nodes and elements of the mesh
    std::vector<Node> _nodes;
    std::vector<ElementType> _elems;

public:
    //standard constructor
    Mesh():
        _nodes(0),
        _elems(0)
        {}
    //getter for the number of Nodes
    const unsigned int& get_nNodes() const{return _nNodes;};
    //getter for number of Elements
    unsigned int get_nElems() const{unsigned int sum = _nElems[0]; for(unsigned int i= 1; i < DIM; i++)sum*=_nElems[i] ; return sum;};
    const unsigned int& get_nEx() const{return _nElems[0];};
    const unsigned int& get_nEy() const{if constexpr (DIM == 2)return _nElems[1];else return 0;};
    // getter for _nodes vector
    const std::vector<Node>& getNodes() const
    {
        return _nodes;
    };
    // getter for _elems vector
    const std::vector<ElementType>& getElems() const
    {
        return _elems;
    };
    //standard printer
    void printMesh(std::ostream& out = std::cout) const
    {
        out <<"_____________________________________"<< std::endl;
        out << "The number of nodes of the mesh is: " << this->get_nNodes() << std::endl; 
        out << "The number of elements of the mesh is: " << this->get_nElems() << std::endl; 
        out << "The coordinates of the nodes are:"<<std::endl;

        out << "Node: Coordinates, Boundary"<< std::endl;
        for(auto i : _nodes)
        {
            i.printNode(out);
        }
        out << "The elements are:"<<std::endl;
        
        for(auto i : _elems)
        {
            i.printElem(out);
        }

    };
    //setter for the number of elements and nodes for a txt file
    void setNum()
    {
        //open a text file and read the parameters form it
        // std::ifstream global_numbers("../mesh/numbers.txt");
        // global_numbers >> _nNodes >> _nElems;
        // global_numbers.close();

        std::string linestr,word;
        std::fstream global_numbers ("../mesh/matlab_mesh_generator/2D/numbers.csv", std::ios::in);
        if(global_numbers.is_open())
        {
            std::getline(global_numbers, linestr); //header - throw away
            std::getline(global_numbers, linestr);
            std::stringstream str(linestr);
            //get the number of nodes
            std::getline(str, word,',');
            std::istringstream(word)>>_nNodes;
            //get the number of elems along the x axis
            std::getline(str,word,',');
            std::istringstream(word)>>_nElems[0];
            if constexpr(DIM == 2)
            {
                //get the number of elems along the y axis
                std::getline(str,word,',');
                std::istringstream(word)>>_nElems[1];
            }
            


        }else
        {
            std::cout << "Could not open the numbers.csv" << std::endl;
            _nNodes = 0;
            std::iota(_nElems.begin(), _nElems.end(), 0);
        }

            

        _nodes.reserve(this->get_nNodes());
        _elems.reserve(this->get_nElems());
        return;

        
        
    };
    //setter for the nodes and elements of the mesh
    //this method reads the number of elements and nodes from a csv file,
    //the nodes coordinates composing the mesh form a csv and finally the elements (to determine how to nodes are connected)
    //for another csv file    
    void setMesh_csv()
    {

        if constexpr (DIM == 1)
        {
            //generate and fill the mesh

            //first set the number of nodes and elements by calling the setNum() method
            setNum();
            
            //define some temp variables
            std::string linestr, word;
            unsigned int holder_id;
            double holderX;

            
            //read nodes coordinates from file
            std::fstream nodes ("../mesh/matlab_mesh_generator/1D/nodes_coordinates.csv", std::ios::in);
            if(nodes.is_open())
            {
                std::getline(nodes, linestr); // header row - throw away
                while(std::getline(nodes, linestr))
                {
                    //move the input date form one line to a stream
                    std::stringstream str(linestr);
                    //get the index of the nodes
                    std::getline(str, word,',');
                    std::istringstream(word)>>holder_id;
                    //get the first coordinate
                    std::getline(str,word,',');
                    std::istringstream(word)>>holderX;
                    //create a temp Node and store it in the vector
                    Node tempnode(holder_id, holderX);
                    _nodes.emplace_back(tempnode);
                }
                //set Dirichelet bc at the extremis
                nodes.close();
            }
            else
            std::cout<<"Could not open the nodes file (nodes_coordinates)\n";


            unsigned int holder1, holder2;
            //read the elements from file
            std::fstream elems ("../mesh/matlab_mesh_generator/1D/elements_vertexes.csv", std::ios::in);
            if(elems.is_open())
            {
                std::getline(elems, linestr); // header row - throw away
                while(std::getline(elems, linestr))
                {
                    //move the input date form one line to a stream
                    std::stringstream str(linestr);
                    //get the index of the nodes
                    std::getline(str, word,',');
                    std::istringstream(word)>>holder_id;
                    //get the first vertex
                    std::getline(str,word,',');
                    std::istringstream(word)>>holder1;
                    //get the second vertex
                    std::getline(str,word,',');
                    std::istringstream(word)>>holder2;
                    
                    //create the element to be added to the list of elements using the paramethers parsed
                    // in the csv file
                    Element_1D temp(holder_id, {_nodes[holder1-1], _nodes[holder2-1]});
                    
                    //finally, we insert the new element in the element vector
                    _elems.emplace_back(temp);
                    
                }
                elems.close();
                
            }
            else
            std::cout<<"Could not open the elements file (element_vertices)\n";

        }
        else if constexpr (DIM == 2)
        {
            //generate and fill the mesh

            //first set the number of nodes and elements by calling the setNum() method
            setNum();
            
            //define some temp variables
            std::string linestr, word;
            unsigned int holder_id;
            double holderX,holderY;

            
            //read nodes coordinates from file
            std::fstream nodes ("../mesh/matlab_mesh_generator/2D/nodes_coordinates.csv", std::ios::in);
            if(nodes.is_open())
            {
                std::getline(nodes, linestr); // header row - throw away
                while(std::getline(nodes, linestr))
                {
                    //move the input date form one line to a stream
                    std::stringstream str(linestr);
                    //get the index of the nodes
                    std::getline(str, word,',');
                    std::istringstream(word)>>holder_id;
                    //get the first coordinate
                    std::getline(str,word,',');
                    std::istringstream(word)>>holderX;
                    //get the second coordinate
                    std::getline(str,word,',');
                    std::istringstream(word)>>holderY;
                    //create a temp Node and store it in the vector
                    Node tempnode(holder_id, holderX, holderY);
                    _nodes.emplace_back(tempnode);
                }
                nodes.close();
            }
            else
            std::cout<<"Could not open the nodes file (nodes_coordinates)\n";


            unsigned int holder1, holder2, holder3, holder4;
            //read the elements from file
            std::fstream elems ("../mesh/matlab_mesh_generator/2D/elements_vertexes.csv", std::ios::in);
            if(elems.is_open())
            {
                std::getline(elems, linestr); // header row - throw away
                while(std::getline(elems, linestr))
                {
                    //move the input date form one line to a stream
                    std::stringstream str(linestr);
                    //get the index of the nodes
                    std::getline(str, word,',');
                    std::istringstream(word)>>holder_id;
                    //get the first vertex
                    std::getline(str,word,',');
                    std::istringstream(word)>>holder1;
                    //get the second vertex
                    std::getline(str,word,',');
                    std::istringstream(word)>>holder2;
                    //get the third vertx
                    std::getline(str,word,',');
                    std::istringstream(word)>>holder3;
                    //get the eventual forth vertex
                    std::getline(str,word,',');
                    std::istringstream(word) >> holder4;
                    
                    //create the element to be added to the list of elements using the paramethers parsed
                    //in the csv file
                    Element_2D temp(holder_id,{_nodes[holder1-1], _nodes[holder2-1],_nodes[holder3-1],_nodes[holder4-1]});
                    //finally we insert the new element in the element vectors
                    _elems.emplace_back(temp);

                    
                }
                elems.close();
                
            }
            else
            std::cout<<"Could not open the elements file (element_vertices)\n";
        }
        else
        {
            std::cout << "DIM is out of range" << std::endl;
            return;
        }
        

        
    };
    
    //member function to generate a basic omogeneous mesh over the domain [a,b]. We must pass some paramethers to the function:
    // ->limits of the domain
    // -> number of segment elements
    void genMesh(const double &xa , const double &xb , const unsigned int &nex )
    {
        _nodes.clear();
        _elems.clear();
        _nNodes = nex +1;
        _nElems[0] = nex;
        if constexpr(DIM == 2)
        {
            _nElems[1] = 0;
        }

        {
            for(unsigned int i = 0; i<nex +1; ++i)
            {
                double x = xa + i*(xb-xa)/nex;
                Node current_node(i+1,x, 0);
                _nodes.emplace_back(current_node);        
            }

            std::vector<std::vector<unsigned int>> temp = this->indexMapping(1,nex);
            for(unsigned int i = 0; i< temp[0].size(); ++i)
            {
                //set the vertexes indexes and coordinates
                Element_1D current_elem(i+1,{_nodes[temp[0][i]-1], _nodes[temp[1][i]-1]});
                _elems.emplace_back(current_elem);
                
            }
            
        }

    };
    //member function to generate a basic omogeneous mesh over the domain [a,b]x[c,d]. We must pass some paramethers to the function
    // to set the mesh generator
    // -> limits of the domain
    // -> number of rectangular elements along the two dimensions
    void genMesh(const double &xa , const double &xb , const double &ya , const double &yb , const unsigned int &nex, const unsigned int &ney)
    {

        _nodes.clear();
        _elems.clear();
            
        _nNodes = (nex + 1)*(ney + 1);
        _nElems[0] = nex;
        if constexpr (DIM == 2)
            _nElems[1] = ney;

        {

            //first generate all the nodes coordinates and later associate index
            //then scan the vector created and attribute a global index to a certain pair of coordinates
            //so as to create a vector and enlarge the vector of nodes and 2D elements
            double x = xa;
            double y = ya;
            std::vector<std::pair<double,double>> temp(0);
            unsigned int global_index(1);
            unsigned int index(0);
            for(unsigned int row = 0; row < ney+1; row ++)
            {
                x = xa;
                y = ya + row*(yb-ya)/ney;
                for(unsigned int col = 0; col < nex +1; col ++)
                {
                    x =xa + col*(xb-xa)/nex;
                    temp.emplace_back(x,y);
                    
                    if(col < 2)
                    {
                        global_index = 1+ col + 2*row;
                    }
                    else
                    {
                        global_index = 1+ col*(ney+1) + row;
                    }

                    Node current_node(global_index, temp[index].first, temp[index].second);
                    _nodes.emplace_back(current_node);
                    
                    

                    index ++;
                    
                }
            }
            std::sort(_nodes.begin(),_nodes.end(),[](const Node &a, const Node &b){return a.getId() < b.getId();});

        }
        

        //now create the elements using indexMapping
        std::vector<std::vector<unsigned int>> temp = this->indexMapping(1,nex,1,ney);
        for(unsigned int i = 0; i< temp[0].size(); ++i)
        {
            //set the vertexes indexes and coordinates
            Element_2D current_elem(i+1,{_nodes[temp[0][i]-1], _nodes[temp[1][i]-1],_nodes[temp[2][i]-1],_nodes[temp[3][i]-1]});
            _elems.emplace_back(current_elem);
            
        }

        
    };
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
                                                const unsigned int &nElementsY = 0/*number of elements along y axis*/)
    {


        if constexpr (DIM == 1)
        {


            unsigned int k (0);
            //define the mapping vector as multidimensional vector having number of rows equal to the number of points per element,
            //and each colum represents one element
            std::vector<std::vector<unsigned int>> nov(degreeX+1);
            for(unsigned int j = 0; j <nElementsX; ++j)
            {
                for(unsigned int i = 1; i <= degreeX + 1; ++i){
                    nov[i-1].emplace_back(k+i);
                }
                k += degreeX;
            }

            return nov;

        }
        else if constexpr (DIM == 2)
        {

            unsigned int npdx = degreeX +1;
            unsigned int npdy = degreeY +1;
            unsigned int ldnov=npdx*npdy;
            unsigned ne=nElementsX*nElementsY;
            std::vector<std::vector<unsigned int>> nov(ldnov, std::vector<unsigned int>(ne, 0));

            //element 1

            for(unsigned int i = 1; i <= ldnov; ++i)
            {
                nov[i-1][0]= i;
            }

            //elements first column

            unsigned int k = ldnov - npdx;
            for(unsigned int ie = 2; ie <= nElementsY; ie++){
                for(unsigned int i = 1; i <= ldnov; ++i)
                {
                    nov[i-1][ie-1] = i+k;
                }
            }
            unsigned int kmax=k+npdx;

            // other columns

            unsigned int nxm1=npdx-1;
            for (unsigned int iex = 2; iex <= nElementsX; ++iex){

                //other rows, bottom elements
                unsigned int ie = (iex - 1)*nElementsY + 1;
                for(unsigned int j = 1; j <= npdy; ++j ){
                    k = (j-1)*npdx;
                    nov[k][ie-1] = nov[j*npdx-1][ie-nElementsY-1];
                    for(unsigned int i = 1; i <= nxm1; ++i)
                    {
                        nov[k+i][ie-1]= kmax+i;
                    }
                    kmax=kmax+nxm1;
                }

                //other elements
                for(unsigned int iey = 2; iey <= nElementsY; ++iey){
                    ie=(iex-1)*nElementsY+iey;

                    //first row
                    for(unsigned int i = 1; i <= npdx; ++i)
                    {
                        nov[i-1][ie-1]=nov[ldnov-npdx+i-1][ie-2]; 
                    }
                    for(unsigned int j = 2; j <= npdy; ++j){
                        k=(j-1)*npdx;
                        nov[k][ie-1]=nov[j*npdx-1][ie-nElementsY-1];
                        for(unsigned int i = 1; i <= nxm1; ++i)
                        {
                            nov[k+i][ie-1]=kmax+i; 
                        }
                        kmax=kmax+nxm1;
                    }
                    
                }

            }

            return nov;

        }
        else
        {
            std::vector<std::vector<unsigned int>> nov;
            return nov; 
        }
        

    };

    // a method to determine wheter an element is on the border
    bool onBorder(ElementType elem)
    {
        if constexpr (DIM == 1)
        {
            if(elem.getId() == 1 || elem.getId() == this->get_nElems()) return true;
            else return false;
        }
        else if constexpr (DIM == 2)
        {
            if(elem.getId() < this->get_nEy()+1  || // elem is on the the left side
               elem.getId() % this->get_nEy() == 0  || //elem is on the top side
               elem.getId() % this-> get_nEy() == 1 || //elem is on the bottom side
               elem.getId() > this->get_nElems() - this->get_nEy() // elem is on the right side
            )
                return true;
            else return false;
        }
    }

    // a method to determine wether a point is on the border
    unsigned short onBorder(Point node)
    {
        if constexpr (DIM == 1)
        {
            if(node.getX() == _nodes[0].getX()|| node.getX() == _nodes[this->get_nNodes()-1].getX())return 1;
            else return 0;
        }
        if constexpr (DIM == 2)
        {
            if(node.getY() == _nodes[0].getY()) return 1;// on bottom side
            else if (node.getY() == _nodes[this->get_nNodes()-1].getY()) return 3;// on top side
            else if (node.getX() == _nodes[0].getX()) return 4; // on left side
            else if (node.getX() == _nodes[this->get_nNodes()-1].getX()) return 2; // on right side
            else return 0;
        }
    }

    ~Mesh() = default;
};



#endif