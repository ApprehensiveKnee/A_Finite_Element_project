#include <iostream>
#include "elements.hpp"


int main(int /*argc*/, char * /*argv*/[]){
    
    Node a(1,0.,0.);
    Node b(2,1.,0.);
    Node c(3,0.,1.);
    std::vector<unsigned int> vertices = {1,2,3};
    Triangle A_triangle(1,vertices);
    Mesh_2D mesh;
    mesh.setMesh();
    mesh.printMesh();
    return 0;
}