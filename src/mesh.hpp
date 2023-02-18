//===========================================================
// HEADER FILE FOR THE CLASS MESH AND DOF HANDLER DEFINITIONS
//===========================================================

#ifndef MESH
#define MESH


#include "elements.hpp"
#include "quadrature.hpp"


//In this case we use forward declaration of the mesh class to get rid of cyclic dependences between
//the dof handler and the mesh class.


//========================================================
//          FORWARD DECLARATION OF CLASS MESH
//========================================================
// Mesh class definition
template<unsigned int DIM>
class Mesh
{
private:
    //number of total nodes of the mesh
    unsigned int _nNodes;
    //number x,y elements of the mesh
    std::array<unsigned int, DIM> _nElems;
    //vectors containing the nodes and elements of the mesh
    std::vector<Node<DIM>> _nodes;
    std::vector<Element<DIM>> _elems;

public:
    //standard constructor
    Mesh():
        _nodes(0),
        _elems(0)
        {}
    // Getter for the number of Nodes
    const unsigned int& get_nNodes() const;
    // Getter for number of Elements
    unsigned int get_nElems() const;
    const unsigned int& get_nEx() const;
    const unsigned int& get_nEy() const;
    // Getter for _nodes vector
    const std::vector<Node<DIM>>& getNodes() const;
    // Getter for _elems vector
    const std::vector<Element<DIM>>& getElems() const;
    // Default printer
    void printMesh(std::ostream& = std::cout ) const;
    // Setter for the number of elements and nodes for a txt file
    void setNum(const std::string&);
    // Setter for the nodes and elements of the mesh
    //this method reads the number of elements and nodes from a csv file,
    //the nodes coordinates composing the mesh form a csv and finally the elements (to determine how to nodes are connected)
    //for another csv file 
    // The first two parameters are used to perform a first resriction on the generality of the Elment class: they set the number of quadrature points
    // along the x and y directions to be the same for all the elements of the mesh (in principly you could still choose different number of quadrature points for the two dimensions) 
    void setMesh_csv(const std::string&,const unsigned int& = r+1, const unsigned int& = r+1);
    // Member function to generate a basic omogeneous mesh over the domain [a,b] (in 1D or [a,b]x[c,d] in 2D). We must pass some parameters to the function:
    // ->limits of the domain
    // -> number of segment elements
    // -> number of quadrature nodes per each element
    void genMesh(const std::array<double,DIM> & , const std::array<double,DIM> & ,const std::array<unsigned int,DIM>&, const std::array<unsigned int,DIM>& = r+1 );
    
    ~Mesh() = default;
};



//declaration of class DofHandler
template<unsigned int DIM>
class DoFHandler
{
private:
    // a member to store the degree of the FE space used along the directions:
    // for full generality, the degree of the space along the two directions could be different,
    // yet we will assume the deg of the spaces is the same for every direction
    std::array<unsigned int, DIM> _r;
    // a vector to store the nodes of the global mesh;
    std::vector<Point<DIM>> _gnodes;
    // a vector to store the ID array: local_index --> global_index, computed using the indexMapping method of the mesh element
    std::vector<std::vector<unsigned int>> _map;
public:
    // Constructor
    DoFHandler(const unsigned int &degX, const unsigned int &degY )
        :_r(DIM==1?std::array<unsigned int,DIM>{degX}:std::array<unsigned int,DIM>{degX,degY}),
        _gnodes(0)
        {};
    // A member to get the degree of the FE space
    const std::array<unsigned int, DIM> getDeg() const
    {
        return _r;
    };
    // A member to get the vector of points
    const std::vector<Point<DIM>>& getPoints() const
    {
        return _gnodes;
    };
    // A member to get the ID array
    const std::vector<std::vector<unsigned int>>& getMap()const
    {
        return _map;
    };
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
                k +=ldnov - npdx;
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


    // A method to compute the number of degree of freedom per element
    unsigned int dof_per_cell() const
    {   
        if constexpr (DIM == 2)
            return (this->getDeg()[0] + 1)*(this->getDeg()[DIM-1] + 1);
        else
            return (this->getDeg()[0] + 1);
    };

    // A member to generate the coordiates of all the points relevant to compute the solution.
    //this computed information is stored inside a vector of points
    const std::vector<Point<DIM>>& genPoints(const Mesh<DIM>& mesh)
    {
        
        // Reset the global nodes of the mesh...
        _gnodes.clear();

        // ...and initialise the map
        _map =DoFHandler::indexMapping(this->getDeg()[0], mesh.get_nEx(),this->getDeg()[DIM-1], mesh.get_nEy());

        if constexpr(DIM ==1)
        {
            
            const unsigned int degX = this->getDeg()[0];

            //Generate the coordinates of the internal points over the reference element...
            // To do so, define a quadrature object...
            FETools::Quadrature qrx;
            // And compute the LGL nodes on the reference interval(just once  --> see the comments on the 2D code)
            qrx.LGL_quadratures(degX);
            // For the first element only: map the coordinates onto the specific element and
            // emplace the points on the gnodes vector
            for(unsigned int q = 0; q < qrx.getN().size(); ++q)
            {
                double x = mesh.getElems()[0].template directMap(qrx.getN()[q]);
                _gnodes.emplace_back(x,0,0);
            }

            //Then loop over all the other elements, filling starting from the second coordinate
            //(The first coordinate of the internal point of an element is given as the last coordinate for the previous element)
            for(auto elem = mesh.getElems().begin()+1; elem < mesh.getElems().end(); ++elem)
            {
                for(unsigned int q = 0; q < qrx.getN().size(); ++q)
                {   
                    double x = elem.template directMap(qrx.getN()[q]);
                    _gnodes.emplace_back(Point<DIM>(x,0,0));
                }
            }
            //Finally update the boundary condition for the first and last element
            _gnodes[0].setBound(1);
            _gnodes[_gnodes.size()-1].setBound(2);

            return _gnodes;


        }
        else if constexpr(DIM == 2)
        {

            auto [degX,degY] = this->getDeg();

            // Get the number of elements using the corresponing methods of the mesh class
            unsigned int nEx = mesh.get_nEx();
            unsigned int nEy = mesh.get_nEy();

            // Reserve the space for the elements of the vector:
            _gnodes.resize(((nEx+1)+(degX-1)*(nEx))*((nEy+1)+(degY-1)*(nEy)),Point<DIM>());

            // Generate the vector of coordinates to store the internal nodes on the reference element
            std::vector<std::array<double,DIM>> coordinates (0);
            // Declare a quadrature objects to compute the LGL nodes along the x and y direction
            FETools::Quadrature qrx;
            FETools::Quadrature qry;
            // Compute the quadrature points (the internal nodes are computed based off
            // the degree of the FE space used: since the value of r is the same for all the 
            // elements, we compute the LGL nodes just once)
            qrx.LGL_quadratures(degX+1);
            qry.LGL_quadratures(degY+1);
            
            
            // In the 2D the setting of the boundary condition is not as simple as the 1D case: in other to compute those,
            // we extract the information on the verteces coordinates and eventually impose the boundary flag if the coordinates of the nodes computed
            // are closer to at least on of those coordinates by less than a se tolerance
            Node<DIM> v1 =mesh.getElems()[0].getNodes()[0]; // v1
            Node<DIM> v2 =mesh.getElems()[nEy*(nEx-1)].getNodes()[1];// v2
            Node<DIM> v3 =mesh.getElems()[nEy*nEx-1].getNodes()[3]; //v3
            //Node v4 =getElems()[nEy-1].getNodes[2]; //v4
            
            // Loop over all the elements of the mesh...
            for(const Element<DIM>& elem: mesh.getElems())
            {
                //... and compute the coordinates of the internal nodes of the elements using mapping operator for each point in the coordiantes vector

                for(unsigned int i = 0; i<qrx.getN().size(); ++i)
                {
                    for(unsigned int j = 0; j <qry.getN().size(); j++)
                    {

                        unsigned int q = (j)*(qrx.getN().size()) + i;
                        unsigned int index = this->getMap()[q][elem.getId()-1]; //get the global index of the point 
                        // then update the global mesh after computing the coordinates of the point using the direct map operator of the vector
                        double x,y;
                        std::tie(x,y)=elem.template directMap<std::tuple<double,double>>(qrx.getN()[i], qry.getN()[j]);
                        // now assign the boundary flag on the basis of the comparison with the coordinates of vertexes:
                        // this is done only if the point considered is not a vertex of the single element considered, for which the boundary flag
                        // are known already
                        unsigned short b(0);
                        if(q == 0/*v1_local*/ || q == degX /*v2_local*/|| q == (degX+1)*(degY) /*v3_local*/|| q == (degX+1)*(degY+1)-1 /*v4_local*/)
                        {
                            unsigned int local_vert = (q == 0)*0 + (q == degX)*1 + (q == (degX+1)*(degY))*2 + (q == (degX+1)*(degY+1)-1)*3;
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

                        _gnodes.at((index-1)) = Point<DIM>(x,y,b);

                    }
                    

                    
                    
                }

            }

            return _gnodes;

        }

        
        
    }

    void printPoints(std::ostream& out = std::cout) const
    {
        out << " The points of the global mesh are:\n" << std::endl;
        
        for(const Point<DIM>& p : this->getPoints())
        {
            out << "======================" << std::endl;
            p.printPoint(out);
            out << "======================" << std::endl;
        }
    };

};



// ============== MESH CLASS METHODS DEFINITIONS ============= 
template<unsigned int DIM>  
const unsigned int& Mesh<DIM>::get_nNodes() const{return _nNodes;};

template<unsigned int DIM> 
unsigned int Mesh<DIM>::get_nElems() const{unsigned int sum = _nElems[0]; for(unsigned int i= 1; i < DIM; i++)sum*=_nElems[i] ; return sum;};

template<unsigned int DIM> 
const unsigned int& Mesh<DIM>::get_nEx() const{return _nElems[0];};

template<unsigned int DIM>  
const unsigned int& Mesh<DIM>::get_nEy() const{if constexpr (DIM == 2)return _nElems[1];else return 0;};

template<unsigned int DIM> 
const std::vector<Node<DIM>>& Mesh<DIM>::getNodes() const
{
    return _nodes;
};

template<unsigned int DIM> 
const std::vector<Element<DIM>>& Mesh<DIM>::getElems() const
{
    return _elems;
};

template<unsigned int DIM> 
void Mesh<DIM>::printMesh(std::ostream& out) const
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

template<unsigned int DIM> 
void Mesh<DIM>::setNum(const std::string& file_name)
{
    //open a text file and read the parameters form it
    // std::ifstream global_numbers("../mesh/numbers.txt");
    // global_numbers >> _nNodes >> _nElems;
    // global_numbers.close();

    std::string linestr,word;
    std::fstream global_numbers ("../mesh/matlab_mesh_generator/meshes/numbers_"+file_name+".csv", std::ios::in);
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
        std::cout << "Could not open the numbers_"+file_name+".csv" << std::endl;
        _nNodes = 0;
        std::iota(_nElems.begin(), _nElems.end(), 0);
    }

        

    _nodes.reserve(this->get_nNodes());
    _elems.reserve(this->get_nElems());
    return;

    
    
};

template<unsigned int DIM> 
void Mesh<DIM>::setMesh_csv(const std::string& file_name,const unsigned int& nqx , const unsigned int& nqy)
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
        setNum(file_name);
        
        //define some temp variables
        std::string linestr, word;
        unsigned int holder_id;
        double holderX;

        
        //read nodes coordinates from file
        std::fstream nodes ("../mesh/matlab_mesh_generator/1D/nodes_coordinates_"+file_name+".csv", std::ios::in);
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
                Node<DIM> tempnode(holder_id, holderX, 0);
                _nodes.emplace_back(tempnode);
            }
            //set boundary flag at the extremis of the considered interval
            _nodes[0].setBound(1);
            _nodes[this->get_nNodes()-1].setBound(2);
            nodes.close();
        }
        else
        std::cout<<"Could not open the nodes file (nodes_coordinates_ "+ file_name +")\n";


        unsigned int holder1, holder2;
        //read the elements from file
        std::fstream elems ("../mesh/matlab_mesh_generator/1D/elements_vertexes_"+file_name+".csv", std::ios::in);
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
                Element<DIM> temp(holder_id, {_nodes[holder1-1], _nodes[holder2-1]}, 0, nqx, 0);
                
                //finally, we insert the new element in the element vector
                _elems.emplace_back(temp);
                
            }
            // Set the element's boundary flags at the extremis
            _elems[0].setBound(1);
            _elems[this->get_nElems()-1].setBound(2);
            elems.close();
            
        }
        else
        std::cout<<"Could not open the elements file (element_vertices_"+file_name+".csv)\n";

    }
    else if constexpr (DIM == 2)
    {
        //generate and fill the mesh

        //first set the number of nodes and elements by calling the setNum() method
        setNum(file_name);
        
        //define some temp variables
        std::string linestr, word;
        unsigned int holder_id;
        double holderX,holderY;
        unsigned short boundary;

        
        //read nodes coordinates from file
        std::fstream nodes ("../mesh/matlab_mesh_generator/meshes/nodes_coordinates_"+file_name+".csv", std::ios::in);
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
                Node<DIM> tempnode(holder_id, holderX, holderY, boundary);
                _nodes.emplace_back(tempnode);
            }
            nodes.close();
        }
        else
        std::cout<<"Could not open the nodes file (nodes_coordinates_"+file_name+".csv)\n";


        unsigned int holder1, holder2, holder3, holder4;
        //read the elements from file
        std::fstream elems ("../mesh/matlab_mesh_generator/meshes/elements_vertexes_"+file_name+".csv", std::ios::in);
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
                Element<DIM> temp(holder_id,{_nodes[holder1-1], _nodes[holder2-1],_nodes[holder3-1],_nodes[holder4-1]}, boundary, nqx, nqy);
                //finally we insert the new element in the element vectors
                _elems.emplace_back(temp);

                
            }
            elems.close();
            
        }
        else
        std::cout<<"Could not open the elements file (element_vertice_"+file_name+".csvs)\n";
    }
    else
    {
        std::cout << "DIM is out of range" << std::endl;
        return;
    }
    

    
};

template<unsigned int DIM> 
void Mesh<DIM>::genMesh(const std::array<double,DIM>& x,
                        const std::array<double,DIM>& y,
                        const std::array<unsigned int,DIM>& ne,
                        const std::array<unsigned int,DIM>& nq)
{
    _nodes.clear();
    _elems.clear();
    
    if constexpr(DIM ==1)
    {
        _nNodes = ne[0] +1;
        _nElems[0] = ne[0];
        for(unsigned int i = 0; i<ne[0] +1; ++i)
        {
            double p = x[0] + i*(x[1]-x[0])/ne[0];
            Node<DIM> current_node(i+1,p, 0);
            _nodes.emplace_back(current_node);        
        }

        //once the nodes have been created, impose the boundary flag onto the nodes at the extremis
        _nodes[0].setBound(1);
        _nodes[this->get_nNodes() -1].setBound(2);

        std::vector<std::vector<unsigned int>> temp = DoFHandler<DIM>::indexMapping(1,ne[0]);
        for(unsigned int i = 0; i< temp[0].size(); ++i)
        {
            //set the vertexes indexes and coordinates
            Element<DIM> current_elem(i+1,{_nodes[temp[0][i]-1], _nodes[temp[1][i]-1]}, 0, nq[0]);
            _elems.emplace_back(current_elem);
            
        }

        //onche the elements have been created, impose the boundary flag onto the elements at the extremis
        _elems[0].setBound(1);
        _elems[this->get_nElems() -1].setBound(2);
    }
    else
    {
        _nodes.clear();
        _elems.clear();
        
        _nNodes = (ne[0] + 1)*(ne[1] + 1);
        _nElems[0] = ne[0];
        if constexpr (DIM == 2)
            _nElems[1] = ne[1];

        {

            //first generate all the nodes coordinates and later associate index
            //then scan the vector created and attribute a global index to a certain pair of coordinates
            //so as to create a vector and enlarge the vector of nodes and 2D elements
            double p = x[0];
            double q = y[0];
            std::vector<std::pair<double,double>> temp(0);
            unsigned int global_index(1);
            unsigned int index(0);
            for(unsigned int row = 0; row < ne[1]+1; row ++)
            {
                p = x[0];
                q = y[0] + row*(y[1]-y[0])/ne[1];
                for(unsigned int col = 0; col < ne[0] +1; col ++)
                {
                    p =x[0] + col*(x[1]-x[0])/ne[0];
                    temp.emplace_back(p,q);
                    
                    if(col < 2)
                    {
                        global_index = 1+ col + 2*row;
                    }
                    else
                    {
                        global_index = 1+ col*(ne[1]+1) + row;
                    }
                    
                    Node<DIM> current_node(global_index, temp[index].first, temp[index].second, 0);

                    if(temp[index].second == y[0])// on bottom side
                    {
                        // check if it is vertex V1
                        if(temp[index].first == x[0])
                        {
                            current_node.setBound(31-1);
                        }
                        // check if it is vertex V2
                        else if(temp[index].first == x[0])
                        {
                            current_node.setBound(31-2);
                        }
                        // otherwise
                        else
                        {
                            current_node.setBound(1);
                        }

                    }
                    else if (temp[index].second == y[1])// on top side
                    {
                        // check if it is vertex V3
                        if(temp[index].first == x[0])
                        {
                            current_node.setBound(31-3);
                        }
                        // check if it is vertex V4
                        else if(temp[index].first == x[1])
                        {
                            current_node.setBound(31-4);
                        }
                        // otherwise
                        else
                        {
                            current_node.setBound(3);
                        }

                    }
                    else if (temp[index].first == x[0]) // on left side
                    {
                        current_node.setBound(4);
                    }
                    else if (temp[index].first == x[1]) // on right side
                    {
                        current_node.setBound(2);
                    }
                    _nodes.emplace_back(current_node);
                    

                    index ++;
                    
                }
            }
            std::sort(_nodes.begin(),_nodes.end(),[](const Node<DIM> &a, const Node<DIM> &b){return a.getId() < b.getId();});

        }
        

        //now create the elements using indexMapping of the DoF Handler class
        std::vector<std::vector<unsigned int>> temp =DoFHandler<DIM>::indexMapping(1,ne[0],1,ne[1]);
        for(unsigned int i = 0; i< temp[0].size(); ++i)
        {
            unsigned int id = i+1;
            Element<DIM> current_elem(id,{_nodes[temp[0][i]-1], _nodes[temp[1][i]-1],_nodes[temp[2][i]-1],_nodes[temp[3][i]-1]}, 0, nq[0], nq[1]);
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

    }

};


    


#endif
