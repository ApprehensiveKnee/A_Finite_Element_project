#include <iostream>
#include "elements.hpp"


int main(int /*argc*/, char * /*argv*/[]){

    Node a(1, 0., 0.);
    Node b(2, 1., 0.);
    Node c(3, 1., 1.);
    Node d(4, 0., 1.);
    Rectangle r(1, 1, 2, 3, 4);
    std::vector<std::vector<double>> temp;
    temp = a.nodes2D(2, b);
    for(auto i : temp){
        std::cout << "--------------" << std::endl;
        for(auto j : i){
            std::cout << " " << j << " ";
        }   
        std::cout << std::endl;
    }
    std::vector<std::vector<unsigned int>> temp1;
    temp1 = Mesh_2D::indexMapping(1,2,1,1);
    for(auto i : temp1){
        for(auto j : i){
            std::cout << " " << j << " ";
        }
        std::cout << std::endl;
    }
    Mesh_2D mesh;
    mesh.setMesh();
    mesh.printMesh();
    return 0;
}