#include "elements.hpp"

// ______________ DEFINITIONS FOR NODE MEMBER FUNC. ______________


unsigned int Node::getId() const{return _id;};

std::array<double,DIM> Node::getCoord() const{return _coord;};

double Node::getX() const{return _coord[0];};

double Node::getY() const{return _coord[1];};

unsigned short Node::getBound() const{return _boundary;};

void Node::setX(const double &r){_coord[0]=r; return;};

void Node::setY(const double &s){_coord[1]=s; return;};

void Node::setBound(const unsigned short &b){_boundary=b; return;};

void Node::printNode(std::ostream &out) const
{
    auto temp = this->getCoord();
    out << this->getId() << ": ";
    std::copy( &temp[0], &temp[temp.size()], std::ostream_iterator<double>(out, ", "));
    out << "boundary: " << this->getBound();
    out<<std::endl;
    return;
};


std::vector<std::array<double,DIM>> Node::nodes(const unsigned int &n, const Node &v)
{
    std::vector<std::array<double,DIM>> temp(Segment::internalN(n)); // vector containing the coordinates of the nodes along a certain dimension of the element
    std::array<double,DIM> v1= this->getCoord()/*get the coordinates of vert1*/;
    std::array<double,DIM> v2= v.getCoord()/*get the coordinates of vert2*/;

    if(DIM == 1)
    {

        double hx = std::abs(v1[0] - v2[0])/n;// step over x
        for(unsigned int i = 0; i< n; i ++) /*in this case, the last element v2 is manually inserted as to avoid round-off error discrepancies*/
        {
        temp[i][0] = v1[0] + hx*i;
        }
        temp[n][0]=v2[0];
        return temp;

    }
    else if(DIM == 2)
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

std::array<double,DIM> Node::classInit(const double &x, const double &y)
{
    std::array<double,DIM> temp;
    if(DIM == 1)
    {
        temp[0]= x;
        return temp;
    }
    else if(DIM == 2)
    {
        temp[0]= x;
        temp[1]= y;
        return temp;
    }

};


//______________________________________________________________



// ______________ DEFINITIONS FOR 1D ELEMENT MEMBER FUNC. ______________


unsigned int Segment::getId() const{return _element_id;};

std::array<unsigned int,2> Segment::getVert() const{return _vert;};

std::array<Node,2> Segment::getNodes() const{return _nodes;};

void Segment::setNodes(const std::array<Node,2> &nodes){_nodes = nodes; return;};

void Segment::printElem(std::ostream &out) const
{
    //print the global index of the element
    out<<"Vertices of Element #"<<this->getId()<<" are:" << std::endl;
    out << "----------"<<std::endl;
    for (auto i : _vert)
    {
        //print the ids of the nodes of the element
        out <<" "<< i <<" ";
        
    }
    
    out <<std::endl << "----------"<<std::endl;
};


unsigned int Segment::internalN(const unsigned int &nx){return nx + 1;};

double Segment::jacobian() const
{
    double J(_nodes[1].getX() - _nodes[0].getX());
    return J;
};


double Segment::inverse_map(const double &x) const
{
    double J(this->jacobian());
    double x_ref = x/J - _nodes[0].getX();
    return x_ref;
};



//______________________________________________________________


// ______________ DEFINITIONS FOR 2D ELEMENT MEMBER FUNC. ______________


unsigned int Element_2D::getId() const{return _element_id;};


std::vector<unsigned int> Element_2D::getVert() const{return _vert;};

std::vector<Node> Element_2D::getNodes() const{return _nodes;};


void Element_2D::setNodes(const std::vector<Node> &nodes)
{
    _nodes = nodes; 
    return;
};


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

Matrix2d Element_2D::jacobian() const
{
    //x_k = B*x_ref + c where c = v1_k and B is the jacobian J = [v2_k - v1_x , v3_k - v1_k], x = [ x1_k , x2_k]^T, x_ref = [x1_ref, x2_ref]^T

    Matrix2d J(2,2);
    J << _nodes[1].getX() - _nodes[0].getX(),
        _nodes[1].getY() - _nodes[0].getY(),
        _nodes[2].getX() - _nodes[0].getX(),
        _nodes[2].getY() - _nodes[0].getY();
    return J;
};


std::tuple<double,double> Element_2D::inverse_map(const double &x, const double &y) const
{
    Matrix2d J(this->jacobian());
    Vector2d c(2), x_ref(2), x_k(2);
    c << _nodes[0].getX(), _nodes[0].getY();
    x_k << x,y;
    x_ref = J.inverse()*(x_k-c);
    return { x_ref(0)*2 - 1.,x_ref(1)*2 - 1.}; 
    //added -1  and multiplied by 2 ,otherwise the inverse transformation would bring us to
    // the refernce element [0,1]x[0,1]

};

//______________________________________________________________





// ______________ DEFINITIONS FOR MESH MEMBER FUNC. ______________

unsigned int Mesh::get_nNodes() const {return _nNodes;};

unsigned int Mesh::get_nElems() const {return _nElems;};

bool Mesh::get_triatrue() const{return _triatrue;};

void Mesh::printMesh(std::ostream& out) const
{
    out <<"_____________________________________"<< std::endl;
    out << "The number of nodes of the mesh is: " << _nNodes << std::endl; 
    out << "The number of elements of the mesh is: " << _nElems << std::endl; 
    out << "The coordinates of the nodes are:"<<std::endl;
    if(DIM == 1){
            out << "Node: CoorX"<< std::endl;
        for(auto i : _nodes){
            i.printNode(out);
        }
        out << "The elements are:"<<std::endl;
        for(auto i : _segments){
            i.printElem(out);
        }
        return;
    }
    else if(DIM == 2)
    {
        out << "Node: CoorX, CoorY"<< std::endl;
        for(auto i : _nodes){
            i.printNode(out);
        }
        out << "The elements are:"<<std::endl;
        if(this->get_triatrue())
            for(auto i : _triangles){
                i.printElem(out);
            }
        else
            for(auto i : _rectangles){
                i.printElem(out);
            }
        return;
    }
    else{
        out << "DIM is out of range"  << std::endl;
        return;
    }

};


void Mesh::setNum()
{
    //open a text file and read the parameters form it
    // std::ifstream global_numbers("../mesh/numbers.txt");
    // global_numbers >> _nNodes >> _nElems;
    // global_numbers.close();

    

    if(DIM == 1)
    {

        std::string linestr,word;
        std::fstream global_numbers ("../mesh/matlab_mesh_generator/1D/numbers.csv", std::ios::in);
        if(global_numbers.is_open())
        {
            std::getline(global_numbers, linestr); //header - throw away
            std::getline(global_numbers, linestr);
            std::stringstream str(linestr);
            //get the index of the nodes
            std::getline(str, word,',');
            std::istringstream(word)>>_nNodes;
            //get the first coordinate
            std::getline(str,word,',');
            std::istringstream(word)>>_nElems;


        }
        _nodes.reserve(this->get_nNodes());
        _segments.reserve(this->get_nElems());
        return;
    }
    else if(DIM == 2)
    {

        std::string linestr,word;
        std::fstream global_numbers ("../mesh/matlab_mesh_generator/2D/numbers.csv", std::ios::in);
        if(global_numbers.is_open())
        {
            std::getline(global_numbers, linestr); //header - throw away
            std::getline(global_numbers, linestr);
            std::stringstream str(linestr);
            //get the index of the nodes
            std::getline(str, word,',');
            std::istringstream(word)>>_nNodes;
            //get the first coordinate
            std::getline(str,word,',');
            std::istringstream(word)>>_nElems;


        }

        //set the right elements to use for the mesh
        if(_nNodes/_nElems == 4)_triatrue = false;
        else if(_nNodes/_nElems == 3)_triatrue = true;

        _nodes.reserve(this->get_nNodes());
        if(this->get_triatrue())
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
    else
    {
        return;
    }

    
    
}
  

void Mesh::setMesh_csv()
{

    if(DIM == 1)
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
            _nodes[0].setBound(1);
            _nodes[(this->get_nNodes())-1].setBound(1);
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
                
                Segment temp(holder_id, holder1, holder2);
                //then we move the data from the list of nodes (coordinates) into the 
                //elements themselves, so that the can be able to compute the jacobian of the element
                temp.setNodes({_nodes[holder1-1], _nodes[holder2-1]});
                //finally, we insert the new element in the element vector
                _segments.emplace_back(temp);
                
            }
            elems.close();
            
        }
        else
        std::cout<<"Could not open the elements file (element_vertices)\n";

    }
    else if(DIM == 2)
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
                //get the boundary value: determine whether the node is on the boundary
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
                if(this->get_triatrue()){
                    
                    Triangle temp(holder_id, holder1, holder2, holder3);
                    //then we move the data from the list of nodes (coordinates) into the 
                    //elements themselves, so that the can be able to compute the jacobian of the element
                    temp.setNodes({_nodes[holder1-1], _nodes[holder2-1],_nodes[holder3-1]});
                    //finally we insert the new element in the element vectors
                    _triangles.emplace_back(temp);

                }
                else{
                    
                    //get the eventual forth vertex
                    std::getline(str,word,',');
                    std::istringstream(word) >> holder4;

                    Rectangle temp(holder_id, holder1, holder2, holder3, holder4);
                    //then we move the data from the list of nodes (coordinates) into the 
                    //elements themselves, so that the can be able to compute the jacobian of the element
                    temp.setNodes({_nodes[holder1-1], _nodes[holder2-1],_nodes[holder3-1],_nodes[holder4-1]});
                    //finally we insert the new element in the element vectors
                    _rectangles.emplace_back(temp);

                }
                
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
    

    
}

//member function to map the global indexes of the nodes in the messh to the local indexes of the nodes inside in the elements
//  THE GLOBAL INDEX OF THE i LOCAL NODE OF A SPECIFIC ie ELEMENT CAN BE COMPUTED BY ACCESSING TO THE (i-1,ie-1) ELEMENT OF THE 
//  VECTOR (MULTIDIMENSIONAL) RETUNED BY THE FUNCTION
std::vector<std::vector<unsigned int>> Mesh::indexMapping(  const unsigned int &degreeX/*degree of the polynomial along x axis*/, 
                                                const unsigned int &nElementsX/*number of elements along x axis*/,
                                                const unsigned int &degreeY/*degree of the polynomial along y axis*/,
                                                const unsigned int &nElementsY/*number of elements along y axis*/)
{


    if(DIM == 1)
    {

        unsigned int k (0);
        //define the mapping vector as multidimensional vector having number of rows equal to the number of elements
        //and for each element we construct the vector mapping the global index
        std::vector<std::vector<unsigned int>> nov (k);
        for(unsigned int j = 0; j <nElementsX; ++j)
        {
            for(unsigned int i = 1; i <= degreeX + 1; ++i){
                nov[i-1].emplace_back(k+i);
            }
            k += degreeX;
        }

        return nov;

    }
    else if(DIM == 2)
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


//______________________________________________________________