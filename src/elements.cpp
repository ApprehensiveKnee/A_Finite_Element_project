#include "elements.hpp"

// ______________ DEFINITIONS FOR NODE MEMBER FUNC. ______________


unsigned int Node::getId() const{return _id;};

std::vector<double> Node::getCoord() const{return _coord;};

void Node::setX(const double &r){_coord[0]=r; return;};

void Node::setY(const double &s){_coord[1]=s; return;};

void Node::printNode(std::ostream &out) const
{
    auto temp = this->getCoord();
    out << _id << ": ";
    std::copy( &temp[0], &temp[temp.size()], std::ostream_iterator<double>(out, ", "));
    out<<std::endl;
    return;
};


std::vector<std::vector<double>> Node::nodes1D(const unsigned int &nx, const Node &v)
{
    std::vector<std::vector<double>> temp(Segment::internalN(nx)); // vector containing the coordinates of the nodes along a certain dimension of the element
    std::vector<double> v1= this->getCoord()/*get the coordinates of vert1*/;
    std::vector<double> v2= v.getCoord()/*get the coordinates of vert2*/;
    double hx = std::abs(v1[0] - v2[0])/nx;// step over x
    for(unsigned int i = 0; i< nx; i ++) /*in this case, the last element v2x is manually inserted as to avoid round-off error discrepancies*/
    {
       temp[i].emplace_back( v1[0] + hx*i);
    }
    temp[nx].emplace_back(v2[0]);
    return temp;
};

std::vector<std::vector<double>> Node::nodes2D(const unsigned int &n, const Node &v)
{
    std::vector<std::vector<double>> temp(Element_2D::internalN(n)); // vector containing the coordinates of the nodes along a certain dimension of the element
    std::vector<double> v1= this->getCoord()/*get the cordinates of vert1*/;
    std::vector<double> v2= v.getCoord()/*get the coordinates of vert2*/;
    double hx = std::abs(v1[0] - v2[0])/n;// step over x
    double hy = std::abs(v1[1] - v2[1])/n;// step over y

    for(unsigned int i = 0; i<n ; ++i) /*in this case, the last element v2x is manually inserted as to avoid round-off error discrepancies*/
    {
       temp[i].emplace_back(v1[0] + hx*i);
       temp[i].emplace_back(v1[1] + hy*i);
    }

    temp[n].emplace_back(v2[0]);
    temp[n].emplace_back(v2[1]);
    return temp;
};



//______________________________________________________________



// ______________ DEFINITIONS FOR 1D ELEMENT MEMBER FUNC. ______________


unsigned int Segment::getId() const{return _element_id;};

std::vector<unsigned int> Segment::getVert() const{return _vert;};

void Segment::printElem(std::ostream &out) const
{
    //print the global index of the element
    out<<"Vertices of Element #"<<_element_id<<" are:" << std::endl;
    out << "----------"<<std::endl;
    for (auto i : _vert)
    {
        //print the ids of the nodes of the element
        out <<" "<< i <<" ";
        
    }
    
    out <<std::endl << "----------"<<std::endl;
};


unsigned int Segment::internalN(const unsigned int &nx){return nx + 1;};



//______________________________________________________________


// ______________ DEFINITIONS FOR 2D ELEMENT MEMBER FUNC. ______________


unsigned int Element_2D::getId() const{return _element_id;};


std::vector<unsigned int> Element_2D::getVert() const{return _vert;};


void Element_2D::printElem(std::ostream &out) const
{
    //print the global index of the element
    out<<"Vertices of Element #"<<_element_id<<" are:" << std::endl;
    out << "----------"<<std::endl;
    for (auto i : _vert)
    {
        //print the ids of the nodes of the element
        out <<" "<< i <<" ";
        
    }
    
    out <<std::endl << "----------"<<std::endl;
};

unsigned int Element_2D::internalN(const unsigned int &nx){return (nx +1) ;};

//______________________________________________________________




// ______________ DEFINITIONS FOR 1D MESH MEMBER FUNC. ______________


unsigned int Mesh_1D::get_nNodes() const{return _nNodes;};

unsigned int Mesh_1D::get_nElems() const{return _nElems;};

void Mesh_1D::printMesh(std::ostream& out) const
{
    out <<"_____________________________________"<< std::endl;
    out << "The number of nodes of the mesh is: " << _nNodes << std::endl; 
    out << "The number of elements of the mesh is: " << _nElems << std::endl; 
    out << "The coordinates of the nodes are:"<<std::endl;
    out << "Node: CoorX"<< std::endl;
    for(auto i : _nodes){
        i.printNode(out);
    }
    out << "The elements are:"<<std::endl;
    for(auto i : _segments){
        i.printElem(out);
    }
    

};

void Mesh_1D::setNum()
{
    //open a text file and read the parameters form it
    std::ifstream global_numbers("../mesh/numbers.txt");
    global_numbers >> _nNodes >> _nElems;
    global_numbers.close();

    _nodes.reserve(this->get_nNodes());
    _segments.reserve(this->get_nElems());
        
    
};
  
void Mesh_1D::setMesh()
{
    //generate and fill the mesh

    //first set the number of nodes and elements by calling the setNum() method
    setNum();
    
    //define some temp variables
    std::string linestr, word;
    unsigned int holder_id;
    double holderX;
    bool boundary;

    
    //read nodes coordinates from file
    std::fstream nodes ("../mesh/nodes_coordinates.csv", std::ios::in);
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
            //get the boundary value: determine whether the node is on the boundary
            std::getline(str,word,',');
            std::istringstream(word)>>boundary;
            //create a temp Node and store it in the vector
            Node tempnode(holder_id, holderX);
            _nodes.emplace_back(tempnode);
        }
        nodes.close();
    }
    else
    std::cout<<"Could not open the nodes file (nodes_coordinates)\n";


    unsigned int holder1, holder2;
    //read the elements from file
    std::fstream elems ("../mesh/elements_vertexes.csv", std::ios::in);
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
            
            Segment temp(holder_id, holder1, holder2);
            _segments.emplace_back(temp);
            
        }
        elems.close();
        
    }
    else
    std::cout<<"Could not open the elements file (element_vertices)\n";

    
};

std::vector<std::vector<unsigned int>> Mesh_1D::indexMapping(const unsigned int &degree, const unsigned int &nElements)
{
    unsigned int k (0);
    //define the mapping vector as multidimensional vector having number of rows equal to the number of elements
    //and for each element we construct the vector mapping the global index
    std::vector<std::vector<unsigned int>> temp (k);
    for(unsigned int j = 0; j <nElements; ++j)
    {
        for(unsigned int i = 1; i <= degree + 1; ++i){
            temp[i-1].emplace_back(k+i);
        }
        k += degree;
    }

    return temp;

};

//______________________________________________________________



// ______________ DEFINITIONS FOR 2D MESH MEMBER FUNC. ______________

unsigned int Mesh_2D::get_nNodes() const {return _nNodes;};

unsigned int Mesh_2D::get_nElems() const {return _nElems;};

void Mesh_2D::printMesh(std::ostream& out) const
{
    out <<"_____________________________________"<< std::endl;
    out << "The number of nodes of the mesh is: " << _nNodes << std::endl; 
    out << "The number of elements of the mesh is: " << _nElems << std::endl; 
    out << "The coordinates of the nodes are:"<<std::endl;
    out << "Node: CoorX, CoorY"<< std::endl;
    for(auto i : _nodes){
        i.printNode(out);
    }
    out << "The elements are:"<<std::endl;
    if(_triatrue)
        for(auto i : _triangles){
            i.printElem(out);
        }
    else
        for(auto i : _rectangles){
            i.printElem(out);
        }

};

void Mesh_2D::setNum()
{
    //open a text file and read the parameters form it
    std::ifstream global_numbers("../mesh/numbers.txt");
    global_numbers >> _nNodes >> _nElems;
    global_numbers.close();


    //set the right elements to use for the mesh
    if(_nNodes/_nElems == 4)_triatrue = false;
    else if(_nNodes/_nElems == 3)_triatrue = true;

    _nodes.reserve(this->get_nNodes());
    if(_triatrue)
    {
        _triangles.reserve((this->get_nElems()));
        return;
    }
    else
    {
        _rectangles.reserve(this->get_nElems());
        return;
    }
        
    
}

void Mesh_2D::setMesh()
{
    //generate and fill the mesh

    //first set the number of nodes and elements by calling the setNum() method
    setNum();
    
    //define some temp variables
    std::string linestr, word;
    unsigned int holder_id;
    double holderX,holderY;

    
    //read nodes coordinates from file
    std::fstream nodes ("../mesh/nodes_coordinates.csv", std::ios::in);
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
            std::getline(str,word);
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
    std::fstream elems ("../mesh/elements_vertexes.csv", std::ios::in);
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
            if(!_triatrue){
                //get the eventual forth vertex
                std::getline(str,word,',');
                std::istringstream(word) >> holder4;

                Rectangle temp(holder_id, holder1, holder2, holder3, holder4);
                _rectangles.emplace_back(temp);
            }
            else{
                Triangle temp(holder_id, holder1, holder2, holder3);
                _triangles.emplace_back(temp);

            }
            
        }
        elems.close();
        
    }
    else
    std::cout<<"Could not open the elements file (element_vertices)\n";

    
}

//member function to map the global indexes of the nodes in the messh to the local indexes of the nodes inside in the elements
//  THE GLOBAL INDEX OF THEE i LOCAL NODE OF A SPECIFIC ie ELEMENT CAN BE COMPUTED BY ACCESSING TO THE (i-1,ie-1) ELEMENT OF THE 
//  VECTOR (MULTIDIMENSIONAL) RETUNED BY THE FUNCTION
std::vector<std::vector<unsigned int>> Mesh_2D::indexMapping(  const unsigned int &degreeX/*degree of the polynomial along x axis*/, 
                                                const unsigned int &nElementsX/*number of elements along x axis*/,
                                                const unsigned int &degreeY/*degree of the polynomial along y axis*/,
                                                const unsigned int &nElementsY/*number of elements along y axis*/)
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

};


//______________________________________________________________