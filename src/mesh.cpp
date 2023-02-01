
#include "mesh.hpp"

//______________________DOF HANDLER CLASS MEMBERS______________________


const std::array<unsigned int, DIM> DoFHandler::getDeg() const{return _r;};

const std::vector<Point>& DoFHandler::getPoints() const{return _gnodes;};

const std::vector<std::vector<unsigned int>>& DoFHandler::getMap()const{return _map;};

unsigned int DoFHandler::dof_per_cell() const
{   
    if constexpr (DIM == 2)
        return (this->getDeg()[0] + 1)*(this->getDeg()[DIM-1] + 1);
    else
        return (this->getDeg()[0]);
};

const std::vector<std::vector<unsigned int>> DoFHandler::indexMapping(const unsigned int &degreeX, 
                                                                                const unsigned int &nElementsX,
                                                                                const unsigned int &degreeY,
                                                                                const unsigned int &nElementsY)
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


void DoFHandler::printPoints(std::ostream& out) const
{
    out << " The points of the global mesh are:\n" << std::endl;
    
    for( Point p : this->getPoints())
    {
        out << "======================" << std::endl;
        p.printPoint(out);
        out << "======================" << std::endl;
    }
};


std::array<unsigned int, DIM> DoFHandler::classInit(const unsigned int& rx, const unsigned int& ry)
{
    if constexpr(DIM == 2)
    {
        return {rx,ry};
    }
    else
    {
        return {rx};
    }
}

//_______________________________________________________________

// ______________ DEFINITIONS FOR MESH MEMBER FUNC. ______________


// the definitions for the mesh class members are contained in the relative header file 
// the class itself was declared as template


//______________________________________________________________