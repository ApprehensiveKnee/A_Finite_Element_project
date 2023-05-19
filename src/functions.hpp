
//===========================================================
//          HEADER FILE FOR THE FUNCTION CLASSES 
//===========================================================
#ifndef FUN
#define FUN
#include "elements.hpp"

// Here some classes that implement mathematical functions
//Pure class for scalar function taking values in 2D domains
    template<unsigned short DIM>
    class Function
    {
    public:
        //standard contructor
        Function(){};
        //define a virtual method that will be overridden by the derived classes
        virtual double value(const Point<DIM> &) const= 0;
        virtual std::array<double, DIM> grad(const Point<DIM>&) const{return std::array<double, DIM>{};}
        //define virtual contructor
        virtual ~Function(){};

    };

    //then some general functions derived from the function base class
    //we hereby consider just scalar functions
    template<unsigned short DIM>
    class DiffusionCoefficient : public Function<DIM>
    {
    public:
        // Constructor.
        DiffusionCoefficient()
        {}

        // Evaluation.
        virtual double
        value(const Point<DIM> & /*p*/) const override
        {
        return 1.0;
        }

    };

    // Reaction coefficient.
    template<unsigned short DIM>
    class ReactionCoefficient : public Function<DIM>
    {
    public:
        // Constructor.
        ReactionCoefficient()
        {}

        // Evaluation.
        virtual double
        value(const Point<DIM>& /*p*/) const override
        {
            if constexpr(DIM == 1)
                return 0.0;
            else
                return 1.0;
        }
    };

    // Forcing term.
    template<unsigned short DIM>
    class ForcingTerm : public Function<DIM>
    {
    public:
        // Constructor.
        ForcingTerm()
        {}

        // Evaluation.
        virtual double
        value(const Point<DIM> & p) const override
        {
            // 1D case
            if(DIM == 1)
            {
                // if(p.getX() <=0.5)
                //         return 0.;
                //     else
                //         return (- std::sqrt(p.getX()-0.5));
                return 4. * M_PI * M_PI * sin(2 * M_PI * p.getX());
            }

            // 2D case
            if(DIM == 2)
                return (20.*M_PI*M_PI + 1.)*std::sin(2*M_PI*p.getX())*std::sin(4*M_PI*p.getY());
        }
    };

    // The exact solution (must be known)
    template<unsigned short DIM>
    class ExactSolution : public Function<DIM>
    {
    public:
        // Constructor.
        ExactSolution()
        {}

        // Evaluation.
        virtual double
        value(const Point<DIM> & p) const override
        {
            // 1D case
            if constexpr(DIM == 1)
            {
                // if(p.getX() <=0.5)
                //     return (-(4./15) * std::pow(0.5,2.5)) * p.getX();
                // else
                //     return (-(4./15) * std::pow(0.5,2.5)) * p.getX() + ((4./15) * std::pow(p.getX() - 0.5,2.5));
                return sin(2 * M_PI * p.getX());
            }
            // 2D case
            else
                return std::sin(2*M_PI*p.getX())*std::sin(4.*M_PI*p.getY());
        }

        std::array<double, DIM> 
        grad(const Point<DIM> &p) const override
        {

            // 1D case
            if constexpr (DIM == 1)
                {
                    // if(p.getX() <=0.5)
                    //     return {(-(4./15) * std::pow(0.5,2.5))};
                    // else
                    //     return {(-(4./15) * std::pow(0.5,2.5)) + ((4./15) * 2.5 * std::pow(p.getX() - 0.5,1.5))};
                    return { 2 * M_PI * cos(2 * M_PI * p.getX())};
                }
                
            // 2D case
            else 
                return{2. * M_PI *std::cos(2.*M_PI*p.getX())*std::sin(4.*M_PI*p.getY()), 4. * M_PI * std::sin(2.*M_PI*p.getX())*std::cos(4.*M_PI*p.getY())};
            
        }
    };

    // Some constant functions to use, for instance, with Dirichelet b.c

    template<unsigned short DIM>
    class functionZero : public Function<DIM>
    {
    public:
        // Constructor.
        functionZero()
        {}

        // Evaluation.
        virtual double
        value(const Point<DIM> & /*p*/) const override
        {
        return 0.0;
        }
    };

    template<unsigned short DIM>
    class functionOne : public Function<DIM>
    {
    public:
        // Constructor.
        functionOne()
        {}

        // Evaluation.
        virtual double
        value(const Point<DIM> & /*p*/) const override
        {
        return 1.0;
        }
    };


#endif

