//===========================================================
// HEADER FILE FOR THE CLASS MESH AND DOF HANDLER DEFINITIONS
//===========================================================
#ifndef MESH
#define MESH


#include "elements.hpp"


//In this case we use forward declaration to get rid of cyclic dependences between classes

//forward declaration for class mesh
// Mesh class definition
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
    const unsigned int& get_nNodes() const;
    //getter for number of Elements
    unsigned int get_nElems() const;
    const unsigned int& get_nEx() const;
    const unsigned int& get_nEy() const;
    // getter for _nodes vector
    const std::vector<Node>& getNodes() const;
    // getter for _elems vector
    const std::vector<ElementType>& getElems() const;
    //standard printer
    void printMesh(std::ostream& = std::cout ) const;
    //setter for the number of elements and nodes for a txt file
    void setNum();
    //setter for the nodes and elements of the mesh
    //this method reads the number of elements and nodes from a csv file,
    //the nodes coordinates composing the mesh form a csv and finally the elements (to determine how to nodes are connected)
    //for another csv file    
    void setMesh_csv();
    
    //member function to generate a basic omogeneous mesh over the domain [a,b]. We must pass some paramethers to the function:
    // ->limits of the domain
    // -> number of segment elements
    void genMesh(const double & , const double & , const unsigned int & );
    //member function to generate a basic omogeneous mesh over the domain [a,b]x[c,d]. We must pass some paramethers to the function
    // to set the mesh generator
    // -> limits of the domain
    // -> number of rectangular elements along the two dimensions
    void genMesh(const double & , const double & , const double & , const double & , const unsigned int &, const unsigned int &);
    
    ~Mesh() = default;
};



//declaration of class DofHandler
class DoFHandler
{
private:
    // a member to store the degree of FE space used along the directions
    std::array<unsigned int, DIM> _r;
    // a vector to store the nodes of the global mesh;
    std::vector<Point> _gnodes;
    // a vector to store the ID array: local_index --> global_index, computed using the indexMapping method of the mesh element
    std::vector<std::vector<unsigned int>> _map;
public:
    //standard constructor
    DoFHandler(const unsigned int &degX = r, const unsigned int &degY = r )
        :_r(classInit(degX, degY)),
        _gnodes(0)
        {};
    // a member to get the degree of the FE space
    const std::array<unsigned int, DIM> getDeg() const;
    // a member to get the vector of points
    const std::vector<Point>& getPoints() const;
    // a member to get the ID array
    const std::vector<std::vector<unsigned int>>& getMap()const ;
    // a member function to generate the ID array
    //->to map the global indexes of the nodes in the messh to the local indexes of the nodes inside in the elements
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
    //

    static const std::vector<std::vector<unsigned int>>indexMapping(const unsigned int &degreeX /*degree of the polynomial along x axis*/, 
                                                const unsigned int &nElementsX/*number of elements along x axis*/,
                                                const unsigned int &degreeY = 0/*degree of the polynomial along y axis*/,
                                                const unsigned int &nElementsY = 0/*number of elements along y axis*/);


    // a method to compute the number of degree of freedom per element
    unsigned int dof_per_cell() const;
        

    //a member to generate the coordiates of all the points relevant to compute the solution.
    //this computed information is stored inside a vector of points
    template <class ElementType>
    const std::vector<Point>& genPoints(const Mesh<ElementType>& mesh)
    {
        
        // Reset the global nodes of the mesh...
        _gnodes.clear();

        // ...and initialise the map
        _map =DoFHandler::indexMapping(this->getDeg()[0], mesh.get_nEx(),this->getDeg()[DIM-1], mesh.get_nEy());

        if constexpr(DIM ==1)
        {
            
            const unsigned int degX = this->getDeg()[0];

            //Generate the coordinates of internal points over the reference element...
            Node first(1,-1.);
            Node second(2,1.);
            const std::vector<std::array<double,DIM>> coordinates(first.nodes(degX, second));
            //Fill the first degX + 1 elements of the vector of points
            for(unsigned int i = 0; i <= degX; ++i)
            {
                // to do so we apply the direct mapping operator of the first element on the coordinates of reference points
                _gnodes.emplace_back(Point(mesh.getElems()[0].direct_map(coordinates[i][0]),0, i == 0 ? 1:0));
            }

            //Then loop over all the other elements, filling starting from the second coordinate
            //(The first coordinate of the internal point of an element is given as the last coordinate for the previous element)
            unsigned int index = degX;
            for(auto elem = mesh.getElems().begin()+1; elem < mesh.getElems().end(); ++elem)
            {
                for(unsigned int i = 1; i <=degX; ++i)
                {   
                    _gnodes.emplace_back(Point(elem->direct_map(coordinates[i][0]), 0,0));
                }
            }
            //Finally update the boundary condition for the last element
            _gnodes[_gnodes.size()-1].setBound(2);

            return _gnodes;


        }
        else if constexpr(DIM == 2)
        {

            auto [degX,degY] = this->getDeg();

            // Get the number of elements using the corresponing methods of the mesh class
            unsigned int nEx = mesh.get_nEx();
            unsigned int nEy = mesh.get_nEy();

            // Generate the ID array
            std::vector<std::vector<unsigned int>> mapping(this->indexMapping(degX, nEx, degY, nEy));

            // Reserve the space for the elements of the vector:
            _gnodes.resize(((nEx+1)+(degX-1)*(nEx))*((nEy+1)+(degY-1)*(nEy)),Point());

            //Generate the vector of coordinates of the nodes on the reference element
            double h = 2./degY;
            std::vector<std::array<double,DIM>> coordinates (0);
            for (unsigned int i = 0; i < degY; ++i)
            {
                //for each new value of i we 'move up' over the rows of the internal nodes
                Node v1(1, -1. , -1. + h*i );
                Node v2(2, 1. , -1. + h*i );
                std::vector<std::array<double,DIM>> coordinatesX = v1.nodes(degX, v2); // generate the coordinates of the inernal nodes on a certain value of y
                // and fill the vector
                for(unsigned int j = 0; j <= degX; ++j)
                {
                    coordinates.emplace_back(coordinatesX[j]);
                }
            }

            //and the last row of internal points is dealt with outside of the loop cycle to avoid round off error
            {
                Node v1(1, -1. , 1. );
                Node v2(2, 1. , 1. );
                std::vector<std::array<double,DIM>> coordinatesX = v1.nodes(degX, v2);
                for(unsigned int j = 0; j <= degX; ++j)
                {
                    coordinates.emplace_back(coordinatesX[j]);
                }

            }

            
            
            // In the 2D the setting of the boundary condition is not as simple as the 1D case: in other to compute those,
            // we extract the information on the verteces coordinates and eventually impose the boundary flag if the coordinates of the nodes computed
            // are closer to at least on of those coordinates by less than a se tolerance
            Node v1 =mesh.getElems()[0].getNodes()[0]; // v1
            Node v2 =mesh.getElems()[nEy*(nEx-1)].getNodes()[1];// v2
            Node v3 =mesh.getElems()[nEy*nEx-1].getNodes()[3]; //v3
            //Node v4 =getElems()[nEy-1].getNodes[2]; //v4
            
            // Loop over all the elements of the mesh...
            for(ElementType elem: mesh.getElems())
            {
                //... and compute the coordinates of the internal nodes of the elements using mapping operator for each point in the coordiantes vector

                for(unsigned int i = 0; i<coordinates.size(); ++i)
                {
                    unsigned int index = mapping[i][elem.getId()-1]; //get the global index of the point 
                    // then update the global mesh after computing the coordinates of the point using the direct map operator of the vector
                    double x,y;
                    std::tie(x,y)=elem.direct_map(coordinates[i][0], coordinates[i][1]);
                    // now assign the boundary flag on the basis of the comparison with the coordinates of vertexes:
                    // this is done only if the point considered is not a vertex of the single element considered, for which the boundary flag
                    // are known already
                    unsigned short b(0);
                    if(i == 0/*v1_local*/ || i == degX /*v2_local*/|| i == (degX+1)*(degY) /*v3_local*/|| i == (degX+1)*(degY+1)-1 /*v4_local*/)
                    {
                        unsigned int local_vert = (i == 0)*0 + (i == degX)*1 + (i == (degX+1)*(degY))*2 + (i == (degX+1)*(degY+1)-1)*3;
                        b = elem.getNodes()[local_vert].getBound();
                    }
                    else
                    {
                        if(std::abs(x - v1.getX()) < tol) // the point is on side 4
                            b = 4;
                        else if(std::abs(x - v2.getX()) < tol) // the point is on side 2
                            b = 2;
                        else if(std::abs(y - v1.getY()) < tol) // the point is on side 1
                            b = 1;
                        else if(std::abs(y - v3.getY()) < tol) // the point is on side 3
                            b = 3;
                        else
                            b = 0;
                    }

                    //Finally update the vector representing the list of nodes in the global mesh

                    _gnodes.at((index-1)) = Point(x,y,b);
                    
                }

            }

            return _gnodes;

        }

        
        
    }

    void printPoints(std::ostream& = std::cout) const;

    std::array<unsigned int, DIM> classInit(const unsigned int& rx, const unsigned int& ry);
};



// ============== MESH CLASS METHODS DEFINITIONS =============  
template<class ElementType>
const unsigned int& Mesh<ElementType>::get_nNodes() const{return _nNodes;};

template<class ElementType>
unsigned int Mesh<ElementType>::get_nElems() const{unsigned int sum = _nElems[0]; for(unsigned int i= 1; i < DIM; i++)sum*=_nElems[i] ; return sum;};

template<class ElementType>
const unsigned int& Mesh<ElementType>::get_nEx() const{return _nElems[0];};

template<class ElementType>
const unsigned int& Mesh<ElementType>::get_nEy() const{if constexpr (DIM == 2)return _nElems[1];else return 0;};


template<class ElementType>
const std::vector<Node>& Mesh<ElementType>::getNodes() const
{
    return _nodes;
};

template<class ElementType>
const std::vector<ElementType>& Mesh<ElementType>::getElems() const
{
    return _elems;
};


template<class ElementType>
void Mesh<ElementType>::printMesh(std::ostream& out) const
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


template<class ElementType>
void Mesh<ElementType>::setNum()
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

template<class ElementType>
void Mesh<ElementType>::setMesh_csv()
{
    // Reset the members of the mesh
    _nNodes= 0;
    _nElems.fill(0);
    _nodes.clear();
    _elems.clear();


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
            //set boundary flag at the extremis of the considered interval
            _nodes[0].setBound(1);
            _nodes[this->get_nNodes()-1].setBound(2);
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
            // Set the element's boundary flags at the extremis
            _elems[0].setBound(1);
            _elems[this->get_nElems()-1].setBound(2);
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
        unsigned short boundary;

        
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
                // get the boundary flag
                std::getline(str,word,',');
                std::istringstream(word)>>boundary;
                //create a temp Node and store it in the vector
                Node tempnode(holder_id, holderX, holderY, boundary);
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
                //get the forth vertex
                std::getline(str,word,',');
                std::istringstream(word) >> holder4;
                //get the boundary flag
                std::getline(str,word,',');
                std::istringstream(word) >> boundary;

                
                //create the element to be added to the list of elements using the paramethers parsed
                //in the csv file
                Element_2D temp(holder_id,{_nodes[holder1-1], _nodes[holder2-1],_nodes[holder3-1],_nodes[holder4-1]}, boundary);
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

template<class ElementType>
void Mesh<ElementType>::genMesh(const double &xa , const double &xb , const unsigned int &nex )
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

        //once the nodes have been created, impose the boundary flag onto the nodes at the extremis
        _nodes[0].setBound(1);
        _nodes[this->get_nNodes() -1].setBound(2);

        std::vector<std::vector<unsigned int>> temp = DoFHandler::indexMapping(1,nex);
        for(unsigned int i = 0; i< temp[0].size(); ++i)
        {
            //set the vertexes indexes and coordinates
            Element_1D current_elem(i+1,{_nodes[temp[0][i]-1], _nodes[temp[1][i]-1]});
            _elems.emplace_back(current_elem);
            
        }

        //onche the elements have been created, impose the boundary flag onto the elements at the extremis
        _elems[0].setBound(1);
        _elems[this->get_nElems() -1].setBound(2);
    }

};


template<class ElementType>
void Mesh<ElementType>::genMesh(const double &xa , const double &xb , const double &ya , const double &yb , const unsigned int &nex, const unsigned int &ney)
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
                
                Node current_node(global_index, temp[index].first, temp[index].second );

                if(temp[index].second == ya)// on bottom side
                {
                    // check if it is vertex V1
                    if(temp[index].first == xa)
                    {
                        current_node.setBound(31-1);
                    }
                    // check if it is vertex V2
                    else if(temp[index].first == xb)
                    {
                        current_node.setBound(31-2);
                    }
                    // otherwise
                    else
                    {
                        current_node.setBound(1);
                    }

                }
                else if (temp[index].second == yb)// on top side
                {
                    // check if it is vertex V3
                    if(temp[index].first == xb)
                    {
                        current_node.setBound(31-3);
                    }
                    // check if it is vertex V4
                    else if(temp[index].first == xa)
                    {
                        current_node.setBound(31-4);
                    }
                    // otherwise
                    else
                    {
                        current_node.setBound(3);
                    }

                }
                else if (temp[index].first == xa) // on left side
                {
                    current_node.setBound(4);
                }
                else if (temp[index].first == xb) // on right side
                {
                    current_node.setBound(2);
                }
                _nodes.emplace_back(current_node);
                

                index ++;
                
            }
        }
        std::sort(_nodes.begin(),_nodes.end(),[](const Node &a, const Node &b){return a.getId() < b.getId();});

    }
    

    //now create the elements using indexMapping
    std::vector<std::vector<unsigned int>> temp =DoFHandler::indexMapping(1,nex,1,ney);
    for(unsigned int i = 0; i< temp[0].size(); ++i)
    {
        unsigned int id = i+1;
        Element_2D current_elem(id,{_nodes[temp[0][i]-1], _nodes[temp[1][i]-1],_nodes[temp[2][i]-1],_nodes[temp[3][i]-1]});
        if(id < this->get_nEy()+1) // elem is on the the left side
        {
            // check for the vetrex V1
            if(id % this-> get_nEy() == 1)
                current_elem.setBound(31-1);
            // check for the vetrex V4
            else if(id % this->get_nEy() == 0)
                current_elem.setBound(31-4);
            // otherwise
            else
                current_elem.setBound(4);
        }
        else if(id > this->get_nElems() - this->get_nEy()) // elem is on the the right side
        {
            // check for the vetrex V2
            if(id % this-> get_nEy() == 1)
                current_elem.setBound(31-2);
            // check for the vetrex V3
            else if(id % this->get_nEy() == 0) 
                current_elem.setBound(31-3);
            // otherwise
            else
                current_elem.setBound(2);
        }
        else if(id % this-> get_nEy() == 1) // elem is on the bottom side
        {
            current_elem.setBound(1);
        }
        else if(id % this->get_nEy() == 0) //elem is on the top side
        {
            current_elem.setBound(3);
        }
    
        _elems.emplace_back(current_elem);
        
    }

    
};
    


#endif
