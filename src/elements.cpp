#include "elements.hpp"

// ______________ DEFINITIONS FOR NODE MEMBER FUNC. ______________


unsigned int Node::getId() const{return _id;};

std::array<double, DIM> Node::getCoord() const{return _coord;};

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


//______________________________________________________________



// ______________ DEFINITIONS FOR ELEMENT MEMBER FUNC. ______________


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


//______________________________________________________________



// ______________ DEFINITIONS FOR MESH MEMBER FUNC. ______________

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
            //get the last vertx
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


//______________________________________________________________