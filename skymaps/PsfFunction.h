/** @file PsfFunction.h
@brief declear PsfFunction.

$Header$
*/

#ifndef skymaps_PsfFunction_h
#define skymaps_PsfFunction_h

#include "astro/SkyDir.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace skymaps{
/**
    @class PsfFunction
    @brief Define the power-law psf function 

*/

class PsfFunction 
{
    public:
        /// @param gamma the power-law factor
        PsfFunction(double g=2.0)
            : m_gamma(g)
            , m_norm(1.-1/m_gamma)
        {}

        /** @brief the function.
            @param r is the reference direction.
            @param r_prime points to a position to be evaluated relative to r. 
            @param sigma angular resolution scale factor 
        */
        double operator () (const astro::SkyDir & r, 
            const astro::SkyDir & r_prime,
            double sigma) const ;

        ///@return the value as a function of the scaled deviation squared
        double operator()(double u)const;
        double gamma()const{return m_gamma;}

        /// @return integral of the function from 0 to umax
        double integral(double umax)const;

        /// @return integral of the square of the function from 0 to umax
        double integralSquare(double umax)const;

        /// @returns a u value generate from the distribution out to umax
        //  @param rand  uniformly distributed random number [0,1)
        //  @param umax  max u value to generate
        double mc(double rand, double umax);

    private:
        double m_gamma;
        double m_norm;
};

}

#endif

