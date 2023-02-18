
//===========================================================
//          HEADER FILE FOR THE FUNCTION CLASSES 
//===========================================================
#ifndef FUN
#define FUN
#include "elements.hpp"

// Here some classes that implement mathematical functions
//Pure class for scalar function taking values in 2D domains
    template<unsigned int DIM>
    class Function
    {
    public:
        //standard contructor
        Function(){};
        //define a virtual method that will be overridden by the derived classes
        virtual double value(const Point<DIM> &) const= 0;
        virtual std::array<double, DIM> grad(const Point<DIM>&) const{return {0.,0.};};
        //define virtual contructor
        virtual ~Function(){};

    };

    //then some general functions derived from the function base class
    //we hereby consider just scalar functions
    template<unsigned int DIM>
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
    template<unsigned int DIM>
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
        return 1.0;
        }
    };

    // Forcing term.
    template<unsigned int DIM>
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
            return (20.*M_PI*M_PI + 1.)*std::sin(2*M_PI*p.getX())*std::sin(4*M_PI*p.getY());
        }
    };

    // The exact solution (must be known)
    template<unsigned int DIM>
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
            return std::sin(2*M_PI*p.getX())*std::sin(4.*M_PI*p.getY());
        }

        std::array<double, DIM> 
        grad(const Point<DIM> &p) const override
        {
            return{2. * M_PI *std::cos(2.*M_PI*p.getX())*std::sin(4.*M_PI*p.getY()), 4. * M_PI * std::sin(2.*M_PI*p.getX())*std::cos(4.*M_PI*p.getY())};
            
        }
    };

    // Some constant functions to use, for instance, with Dirichelet b.c

    template<unsigned int DIM>
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

    template<unsigned int DIM>
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

