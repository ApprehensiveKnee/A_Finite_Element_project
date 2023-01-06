#include <iostream>
#include "elements.hpp"
#include "mat_utilities.hpp"


int main(int /*argc*/, char * /*argv*/[]){

    Node a(1, 3., 3.);
    Node b(2, 4., 3.);
    Node c(3, 3., 4.);
    Node d(4, 4., 4.);
    std::vector<Node> nodes ({a,b,c,d});
    Rectangle r(1, nodes);
    std::vector<std::array<double,DIM>> temp;
    temp = a.nodes(2, b);
    for(auto i : temp){
        std::cout << "--------------" << std::endl;
        for(auto j : i){
            std::cout << " " << j << " ";
        }   
        std::cout << std::endl;
    }
    std::vector<std::vector<unsigned int>> temp1;
    // temp1 = Mesh::indexMapping(1,2,1,1);
    // for(auto i : temp1){
    //     for(auto j : i){
    //         std::cout << " " << j << " ";
    //     }
    //     std::cout << std::endl;
    // }
    Mesh mesh(false);
    std::cout << mesh.get_triatrue() << std::endl;
    mesh.setMesh_csv();
    mesh.printMesh();
    std::cout << std::endl;
    Quadrature qr;
    qr.LGL_quadratures(5);
    std::cout << std::endl;
    auto n = qr.getN();
    SpectralFE fe(4);
    fe._update_current(r);
    
    std::cout << "Dofs per cell:"<<fe._dof_per_cell() << std::endl;

    std::cout << "DEG:" <<fe.getDeg() << std::endl;

    for(unsigned int i = 0; i<fe.getQPoints().size(); ++i)
    {
        std::cout << "Quadrature node #" << i << std::endl;
        fe.getQPoints()[i].printPoint();
        std::cout << "Evaluation of basis function over the node" << std::endl;
        std::cout << fe.shape_value(i,fe.getQPoints()[i])<< std::endl;
        
        std::cout << "---------------------"<< std::endl;
    }

    for(unsigned int i = 0 ; i < fe.getINodes().size(); ++i)
    {
        std::cout << " Internal node #" << i << std::endl;
        fe.getINodes()[i].printPoint();
        std::cout << "Evaluation of basis function over the node" << std::endl;
        std::cout << fe.shape_value(i,fe.getINodes()[i]) << std::endl;
        std::cout << "---------------------"<< std::endl;
    }


    //evaluation of basis function over the the internal nodes,

    

    std::cout << std::endl;


    
    return 0;
}