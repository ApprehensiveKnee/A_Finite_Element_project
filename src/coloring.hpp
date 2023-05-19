#ifndef COLORING.HPP

#define COLORING.HPP

#include "mesh.hpp"



// A function to implement the COLORING ALGORITHM:
// The algortihm is taken after the graph theory coloring algorithm: the main idea 
// is to color the nodes of a graph in such a way that no two adjacent nodes have the same color.
// Similarly, in our case, the algorithm is used to group the element of the mesh, stored in the elems vector of the Mesh object,
// in such a way that no element of the same group shares the same node.
// By doing so, elements of the same group can be processed by different threads concurrently without the need of a mutex
// as no race conflicts will arise.

// After sorting the elements in the different groups, in a way similar to what we have alredy done before in the
// solver class, the algorithm will loop over the elements of the different groups and will compute the local matrices contributes to the 
// global matrix.

// The objective is to keep the number of groups minimal.

// Psudo-code:
// (1.) Loop over the elements e = 1, ne
//      (a) For the element e, find the colour C assigned to neighbours of element e
//          - loop over the nodes of element e, i = 1, n 
//          - C = C + colorsOf(i)
//      (b) Find the unused available colour not in C and assign it to element e
//          - loop over the available colours, c = 1, m 
//          - if c is not in C then EC(e) = c
//
// where EC is the array of colours and colorsOf(i) is the function returning a 
// set of colours assigned to elements sharing the node i


// To implement the funtion, first we have to defined the colorsOf function.

// And before that a function that returns the indexes of the elements sharing the same node.
// For this last one, we try two approaches:  (1) using a map and (2) using a vector to store the indexes of the elements sharing the node.

// (1) Using a map:

// std::unordered_map<unsigned int, std::vector<unsigned int>> getElemIdxSharingNode(const Mesh<DIM> &mesh)
// {
//     std::unordered_map<unsigned int, std::vector<unsigned int>> elemIdxSharingNode;
//     for(const Element<DIM>& elem : mesh.getElems())
//     {
//         for(const unsigned int& point_index : elem.getPoints())
//         {
//             elemIdxSharingNode[point_index].push_back(elem.getId());
//         }
//     }
//     return elemIdxSharingNode;
// }

// (2) Using a vector:
template <unsigned short DIM>
std::vector<std::vector<unsigned int>> getElemIdxSharingNode(const Mesh<DIM> & mesh)
{
    std::vector<std::vector<unsigned int>> elemIdxSharingNode(mesh.nodes.size());
    for (const Element<DIM>& elem : mesh.getElems())
    {
        for (const unsigned int& point_index : elem.getPoints())
        {
            elemIdxSharingNode[point_index].push_back(elem.getId());
        }
    }
    return elemIdxSharingNode;
}


#endif
