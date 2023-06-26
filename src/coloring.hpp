//===========================================================
//     HEADER FILE FOR THE COLORING SUPPORTING FUNCTIONS
//===========================================================

#ifndef COL

#define COL

#include "mesh.hpp"
#include <set>
#include <algorithm>

// A function to implement the COLORING ALGORITHM:
// The algortihm is taken after the graph theory coloring algorithm: the main idea
// is to color the nodes of a graph in such a way that no two adjacent nodes have the same color.
// Similarly, in our case, the algorithm is used to group the element of the mesh, stored in the elems vector of the Mesh object,
// in such a way that no element of the same group shares the same node.
// By doing so, elements of the same group can be processed by different threads concurrently without the need of a mutex
// as no race conflicts arise.

// After sorting the elements in the different groups, in a way similar to what we have alredy done before in the
// solver class, the algorithm will loop over the colors serially and compute the local matrices contributes to  the global matrix
// for the elements belonging to the same group

// The objective is to keep the number of groups minimal.

// Pseudo-code:
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

// To implement the coloring algorithm, first we have to define the colorsOf function.

// And before that a function that returns the indexes of the elements sharing the same node.
// For this last one, we try two approaches:  (1) using a map and (2) using a vector to store the indexes of the elements sharing the node.

// (1) Using a map:

// template <unsigned short DIM>
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
std::vector<std::vector<unsigned int>> getElemIdxSharingPoint(const Mesh<DIM> &mesh, const DoFHandler<DIM> &dof)
{
    std::vector<std::vector<unsigned int>> elemIdxSharingPoint(dof.getPoints().size());
    for (const Element<DIM> &elem : mesh.getElems())
    {
        for (const unsigned int &point_index : elem.getPoints())
        {
            elemIdxSharingPoint[point_index - 1].push_back(elem.getId());
        }
    }
    return elemIdxSharingPoint;
}

// Now we can define the colorsOf function:
template <unsigned short DIM>
std::set<unsigned short> colorsOf(const unsigned int &sharedPoint, const Mesh<DIM> &mesh, const DoFHandler<DIM> &dof)
{
    // Compute the vector to store the indexes of the elements sharing the same node
    std::vector<std::vector<unsigned int>> elemIdxSharingPoint = getElemIdxSharingPoint(mesh, dof);
    // Create the set to store the colors of the elements sharing the same node
    std::set<unsigned short> colors;
    for (const unsigned int &elem_index : elemIdxSharingPoint[sharedPoint - 1])
    {
        colors.insert(mesh.getElems()[elem_index - 1].getColor());
    }
    return colors;
}

// Finally, we can define the greedy coloring function

template <unsigned short DIM>
std::vector<std::vector<unsigned int>> greedyColoring(Mesh<DIM> &mesh, const DoFHandler<DIM> &dof, const unsigned short &maxColors)
{
    // First we have to get the indexes of the elements sharing the same node
    std::vector<std::vector<unsigned int>> elemIdxSharingPoint = getElemIdxSharingPoint(mesh, dof);
    // Then we have to define the vector of colors
    std::vector<std::vector<unsigned int>> colors(maxColors);
    // (1.) Loop over the elements e = 1, ne
    for (Element<DIM> &elem : mesh.modifyElems())
    {
        std::set<unsigned short> colorsOfElem;
        // (a)  For the element e, find the colour C assigned to neighbours of element e
        //- Loop over the nodes of element e, i = 1, n
        for (const unsigned int &point_index : elem.getPoints())
        {
            //- C = C + colorsOf(i)
            std::set<unsigned short> temp1 = colorsOfElem;
            std::set<unsigned short> temp2 = colorsOf(point_index, mesh, dof);
            std::set_union(temp1.begin(), temp1.end(), temp2.begin(), temp2.end(),
                           std::inserter(colorsOfElem, colorsOfElem.begin()));
        }

        // (b)  Find the unused available colour not in C and assign it to element e
        //- loop over the available colours, c = 1, m

        for (unsigned short color = 0; color < maxColors; ++color)
        {
            if (colorsOfElem.find(color) == colorsOfElem.end())
            {
                //- if c is not in C then EC(e) = c
                elem.setColor(color);
                colors[color].push_back(elem.getId());
                break;
            }
        }

        // Final addition to the color vector for color-0-elements
        if (elem.getColor() == 0)
        {
            colors[0].push_back(elem.getId());
        }
    }
    return colors;
}

#endif
