//===========================================================
//      HEADER FILE FOR NODE AND ELEMENT CLASSES
//===========================================================

#ifndef ELEM
#define ELEM

#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <numeric>
#include <iterator>
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>
#include <type_traits>
#include <concepts>
#include <filesystem>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <unsupported/Eigen/SparseExtra>
#include "variables.hpp"

#define COLORING

using namespace Eigen;

//========================================================================================================

// a class to store the coordinates of the relevant points
template <unsigned short DIM>
class Point
{
protected:
    // a member array to store the coordinates
    std::array<double, DIM> _coord;
    // a member to store boundary information:
    // if considering the 2D mesh read form csv or the one generated inside the code, the convetion used is the following:
    //
    //                 side 3
    //            V4 __________ V3
    //              |          |
    //              |          |
    //      side 4  |          |  side 2    // Consider the following mesh
    //              |          |
    //              |__________|
    //            V1            V2
    //                 side 1
    //
    // 0 -> internal node
    // 1 -> node on the bottom border
    // 2 -> node on the right border
    // 3 -> node on the top border
    // 4 -> node on the left border
    // The following numbering is assumed to possibly better handle the conditions on the vertexes
    // 31 -1 -> for V1
    // 31 -2 -> for V2
    // 31 -3 -> for V3
    // 31- 4 -> for V4
    // In 1D the first vertex is simply represented with 1, the last with 2
    unsigned short _boundary;

public:
    Point(){};
    // Constructor
    Point(const double &x, const double &y, const unsigned short &boundary) : _coord(initCoord(x, y)),
                                                                              _boundary(boundary){};
    // Array constructor
    Point(const std::array<double, DIM> coord, const unsigned short &boundary) : _coord(coord),
                                                                                 _boundary(boundary){};

    // Standard getters

    const std::array<double, DIM> &getCoord() const
    {
        return _coord;
    };

    const unsigned short &getBound() const
    {
        return _boundary;
    };

    void setBound(const unsigned short &boundary_flag)
    {
        _boundary = boundary_flag;
        return;
    };

    const double &getX() const
    {
        return _coord[0];
    };

    const double &getY() const
    {
        return _coord[DIM - 1];
    };

    // Default setters
    void setX(const double &x)
    {
        _coord[0] = x;
        return;
    };

    void setY(const double &y)
    {
        if constexpr (DIM == 1)
        {
            std::cout << "The method 'setY' can only be called for 2 dimensional problems" << std::endl;
            return;
        }
        else
        {
            _coord[DIM - 1] = y;
            return;
        }
    };

    // A function to print the info of the point
    void printPoint(std::ostream &out = std::cout) const
    {
        auto temp = this->getCoord();
        if constexpr (DIM == 2)
            out << temp[0] << ", " << temp[1] << std::endl;
        else
            out << temp[0] << std::endl;
        out << "boundary: " << this->getBound() << std::endl;

        return;
    };

    // A static method to correclty initilize the coordinates array
    static std::array<double, DIM> initCoord(const double &x, const double &y)
    {
        if constexpr (DIM == 1)
        {
            return std::array<double, DIM>{x};
        }
        else
        {
            return std::array<double, DIM>{x, y};
        }
    };

    ~Point() = default;
};

//========================================================================================================
template <unsigned short DIM>
class Node : public Point<DIM>
{
private:
    // global id of the node
    unsigned int _id;

public:
    Node(){};
    // Constructor
    Node(const unsigned int &id, const double &x, const double &y, const unsigned short &boundary)
        : Point<DIM>(x, y, boundary),
          _id(id){};
    Node(const unsigned int &id, const std::array<double, DIM> &coord, const unsigned short &boundary)
        : Point<DIM>(coord, boundary),
          _id(id){};
    // Default id getter
    const unsigned int &getId() const
    {
        return _id;
    };

    // Default printer
    void printNode(std::ostream &out = std::cout) const
    {
        out << "==============================" << std::endl;
        out << this->getId() << ": ";
        this->printPoint(out);
        out << "==============================" << std::endl;
        out << std::endl;
        return;
    };

    // A method to return a vector containing the coordinates of all the internal nodes
    // along a certain dimension of the element based off the degree of the polinomial used for such dimension
    std::vector<std::array<double, DIM>> nodes(const unsigned int &n, const Node &v) const
    {
        std::vector<std::array<double, DIM>> temp(n + 1); // vector containing the coordinates of the nodes along a certain dimension of the element
        std::array<double, DIM> v1 = this->getCoord() /*get the coordinates of vert1*/;
        std::array<double, DIM> v2 = v.getCoord() /*get the coordinates of vert2*/;

        if constexpr (DIM == 1)
        {

            double hx = std::abs(v1[0] - v2[0]) / n; // step over x
            for (unsigned int i = 0; i < n; i++)     /*in this case, the last element v2 is manually inserted as to avoid round-off error discrepancies*/
            {
                temp[i][0] = v1[0] + hx * i;
            }
            temp[n][0] = v2[0];
            return temp;
        }
        else
        {
            double hx = std::abs(v1[0] - v2[0]) / n; // step over x
            double hy = std::abs(v1[1] - v2[1]) / n; // step over y

            for (unsigned int i = 0; i < n; ++i) /*in this case, the last element v2 is manually inserted as to avoid round-off error discrepancies*/
            {
                temp[i][0] = v1[0] + hx * i;
                temp[i][1] = v1[1] + hy * i;
            }

            temp[n][0] = v2[0];
            temp[n][1] = v2[1];
            return temp;
        }
    };

    // A corresponding method returning a vector of points:
    // This one also includes additional information on the boundary position on the mesh
    std::vector<Point<DIM>> points(const unsigned int &n, const Node &v) const
    {
        std::vector<Point<DIM>> temp(n + 1); // vector containing the coordinates of the nodes along a certain dimension of the element
        std::array<double, DIM> v1 = this->getCoord() /*get the coordinates of vert1*/;
        unsigned short b1 = this->getBound() /*get the boundary of vert1*/;
        std::array<double, DIM> v2 = v.getCoord() /*get the coordinates of vert2*/;
        unsigned short b2 = v.getBound() /*get the boundary of vert2*/;

        if constexpr (DIM == 1)
        {

            double hx = std::abs(v1[0] - v2[0]) / n; // step over x

            for (unsigned int i = 0; i < n; i++) // in this case, the last element v2 is manually inserted as to avoid round-off error discrepancies
            {
                temp[i] = Point<DIM>(v1[0] + hx * i, 0, (i == 0) ? b1 : 0);
            }
            temp[n] = Point<DIM>(v2, b2);
            return temp;
        }
        else
        {
            double hx = std::abs(v1[0] - v2[0]) / n; // step over x
            double hy = std::abs(v1[1] - v2[1]) / n; // step over y

            for (unsigned int i = 0; i < n; ++i) // in this case, the last element v2 is manually inserted as to avoid round-off error discrepancies
            {
                unsigned short b(0);
                if (b1 == b2)
                {
                    b = b1;
                }
                // ________ SPECIAL CASES FOR FIRST VERTEX ________

                else if (b1 == 31 - 1) // start point == vertex 1
                {
                    // either used to compute points on side 1 or 4
                    if (b2 == 31 - 2)      // last == vertex 2
                        b = 1;             // side 1;
                    else if (b2 == 31 - 4) // last == vertex 4
                        b = 4;             // side 4;
                    else
                        b = b2;
                }
                else if (b1 == 31 - 2) // start point == vertex 2
                {
                    // used to compute points on side 2
                    b = 2;
                }
                else if (b1 == 31 - 4) // start point == vertex 4
                {
                    // used to compute points on side 3
                    b = 3;
                }
                // ________ SPECIAL CASES FOR LAST VERTEX ________

                else if (b2 == 31 - 3) // last point == vertex 3
                {
                    // used to compute either points on side 2 or 3
                    if (b1 == 31 - 4)      // first == vertex 4
                        b = 3;             // side 3
                    else if (b1 == 31 - 3) // first == vertex 2
                        b = 2;             // side 2
                    else
                        b = b1;
                }
                else if (b2 == 31 - 2) // last point == vertex 2
                {
                    // used to compute points on side 1
                    b = 1;
                }
                else if (b2 == 31 - 4) // last point == vertex 4
                {
                    // used to compute points on side 4
                    b = 4;
                }
                else
                {
                    // otherwise b is an internal point
                    b = 0;
                }

                temp[i] = Point<DIM>(v1[0] + hx * i, v1[1] + hy * i, (i == 0) ? b1 : b);
            }

            temp[n] = Point<DIM>(v2, b2);
            return temp;
        }
    };

    // Default destructor
    ~Node() = default;
};

//============================= SOME SIMPLE CONCEPTS ======================================
template <typename T>
concept IsMatrix2 = std::is_same_v<T, Eigen::Matrix2d>;
template <typename T>
concept IsDouble = std::is_same_v<T, double>;
template <typename T>
concept IsTuple = std::is_same_v<T, std::tuple<double, double>>;

//========================================================================================================

template <unsigned short DIM>
class Element
{
private:
    // id of the element
    unsigned int _element_id;
    // vector of Node elements corresponding to vert indexes
    std::vector<Node<DIM>> _nodes;

#ifdef COLORING
    // vector for the internal points of the element, used for the coloring algorithm and computed in the genPoints() method
    // inside the DoF Handler class
    std::vector<unsigned int> _points;

    // a variable to store the color that will be assigned to the element
    unsigned short _color;
#endif

    // the following is a viual representation of the ordering of vertices adopted for the construction
    // of the rectanglular element (2D case) to be used for the SEM mesh
    //
    //
    //            V3 __________ V4
    //              |          |
    //              |          |
    //              |          |
    //              |          |
    //              |__________|
    //            V1            V2
    //

    // a member to store the boundary flag
    // the convetion adopted is the same as the Point boundary flag
    unsigned short _boundary;
    // member to store the number of quadrature points to compute the integrals on the element,both for the x and y dimensions
    // PLEASE NOTE: for the upcoming assumptions taken along the project, the element class would not stictly require to
    // have a member storing the number of quadrature nodes, yet for better clarity and generality, it is
    // assumed that different elements within the same mesh could have different number of quadrature points.
    // This would allow for some elements to be better "represented" than others whilst also making the porcess for
    // further developments on the generality of the code easier.
    std::array<unsigned int, DIM> _nq;

public:
    Element(){};
    // Constructor
    Element(const unsigned int &id,
            const std::vector<Node<DIM>> &nodes,
            const unsigned short &boundary = 0,
            const unsigned int &nx = r + 1,
            const unsigned int &ny = r + 1)
        : _element_id(id),
          _nodes(nodes),
          _boundary(boundary),
          _nq(initNQ(nx, ny))
    {
#ifdef COLORING
        _color = 0;
#endif
    }

    // Getter:gets the id of the element
    const unsigned int &getId() const
    {
        return _element_id;
    };
    // Getter:gets the nodes of the element
    const std::vector<Node<DIM>> &getNodes() const
    {
        return _nodes;
    };

#ifdef COLORING
    // Getter:gets the global index of the internal points of the elemet
    const std::vector<unsigned int> &getPoints() const
    {
        return _points;
    };

    // Setter: allows to modify the points of the element
    std::vector<unsigned int> &modifyPoints()
    {
        return _points;
    };

    // Getter: gets the color assigned to the element
    const unsigned short &getColor() const
    {
        return _color;
    };

    // Setter: sets the color assigned to the element
    void setColor(const unsigned short &color)
    {
        _color = color;
        return;
    };

#endif

    // Getter: gets the boundary flag, to implement a tree search of the boundary nodes
    const unsigned short &getBound() const
    {
        return _boundary;
    };
    // Getter: gets the degree of exactness for the quadrature formula to compute integrals on the element
    const std::array<unsigned int, DIM> &getNQ() const
    {
        return _nq;
    };
    // Setter: sets the boundary flag after the initialization of the element
    void setBound(const unsigned short &boundary_flag)
    {
        _boundary = boundary_flag;
        return;
    };
    // Setter for the nodes array (to set the coordinates)
    void setNodes(const std::vector<Node<DIM>> &nodes)
    {
        _nodes = nodes;
        return;
    };

    // Printer: prints the ids and coordinates of the vertexes
    void printElem(std::ostream &out = std::cout) const
    {
        // print the global index of the element
        out << "Vertices of Element #" << this->getId() << " are:" << std::endl;
        out << "----------" << std::endl;
        for (auto i : _nodes)
        {
            // print the ids of the nodes of the element
            out << " " << i.getId() << " ";
        }

        out << std::endl
            << "----------" << std::endl;
    };

    // a method to compute the jacobian of a specific element, describing the affine tranformation

    template <typename T>
        requires(DIM == 2 && IsMatrix2<T>) || (DIM == 1 && IsDouble<T>)
    T jacobian() const
    {
        if constexpr (DIM == 1)
        {
            double J = (this->getNodes()[1].getX() - this->getNodes()[0].getX()) / 2.;
            return J;
        }
        else
        {
            // x_k = B*x_ref + c where c = v1_k and B is the jacobian J = [v2_k - v1_x , v3_k - v1_k], x = [ x1_k , x2_k]^T, x_ref = [x1_ref, x2_ref]^T

            Matrix2d J(2, 2);
            J << (_nodes[1].getX() - _nodes[0].getX()) / 2.,
                (_nodes[1].getY() - _nodes[0].getY()) / 2.,
                (_nodes[2].getX() - _nodes[0].getX()) / 2.,
                (_nodes[2].getY() - _nodes[0].getY()) / 2.;
            return J;
        }
    };
    // a method to compute the inverse transformation onto the refernce element [-1,1] in 1D or [-1,1]x[-1,1] in 2D
    template <typename T>
        requires(DIM == 2 && IsTuple<T>) || (DIM == 1 && IsDouble<T>)
    T inverseMap(const double &x, const double &y = 0) const
    {
        if constexpr (DIM == 1)
        {
            double J = (this->jacobian<double>());
            double x_ref = x / J - _nodes[0].getX();
            return x_ref - 1;
        }
        else
        {
            Matrix2d J(this->jacobian<Eigen::Matrix2d>());
            Vector2d c(2), x_ref(2), x_k(2);
            c << _nodes[0].getX(), _nodes[0].getY();
            x_k << x, y;
            x_ref = J.inverse() * (x_k - c);
            return {x_ref(0) - 1., x_ref(1) - 1.};
        }
    };
    // a method to compute the direct transformation from the reference element
    template <typename T>
        requires(DIM == 2 && IsTuple<T>) || (DIM == 1 && IsDouble<T>)
    T directMap(const double &x_r, const double &y_r = 0) const
    {
        if constexpr (DIM == 1)
        {
            double J = (this->jacobian<double>());
            double x = J * (x_r + 1.) + this->getNodes()[0].getX();
            return x;
        }
        else
        {
            Matrix2d J(this->jacobian<Eigen::Matrix2d>());
            Vector2d c(2), x_ref(2), x_k(2);
            c << _nodes[0].getX(), _nodes[0].getY();
            x_ref << x_r + 1., y_r + 1.;
            x_k = J * x_ref + c;
            return {x_k(0), x_k(1)};
        }
    };

    // a static member to correcly initialize the quadrature degs
    static std::array<unsigned int, DIM> initNQ(const unsigned int &nx, const unsigned int &ny)
    {
        if constexpr (DIM == 1)
        {
            return {nx};
        }
        else if constexpr (DIM == 2)
        {
            return {nx, ny};
        }
    };

    // Destructor
    virtual ~Element() = default;
};

//______________________________________________________________

#endif