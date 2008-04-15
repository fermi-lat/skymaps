/** @file Background.cxx
    @brief implement Background

$Header$
*/

#include "skymaps/Background.h"

#include <stdexcept>
#include <cmath> 

namespace{  // anonymous namespace for helper classes

}  // anonymous namespace
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using namespace skymaps;

Background::Background(const skymaps::SkySpectrum& diffuse, 
                       const skymaps::SkySpectrum& exposure, int n)
: m_diffuse(diffuse)
, m_exposure(exposure)
, m_n(n)
{
    if( n<2 || (n&2) !=0 || n>20){
        throw std::invalid_argument("Background--bad Simpsons rule count: n must be even, <20");
    }
}

Background::~Background()
{}

double Background::value(const astro::SkyDir& dir, double e)const
{
    double val( m_diffuse(dir, e) * m_exposure( dir, e ) );
    return val;
}

///@brief integral for the energy limits, in the given direction
double Background::integral(const astro::SkyDir& dir, double a, double b)const
{
    double step( log(b/a)/m_n ), // step in log scale
           ratio( exp(step) ), // ratio of energies
           c(2.),              // initial simpsons
           e(a);              //  inital energy

    double result( a*value(dir, a) ); // value at low end
    for( int i = 1; i< m_n; ++i ){
        e *= ratio; // next energy
        c = 6-c;   // toggle Simpsons coefficient
        result += c * e* value( dir, e); 
    }
    
    result += b*value(dir, b);// value at high end.
    return result*step/3.;

}

std::string Background::name()const
{
    return "Background: "+ m_diffuse.name()+"+" + m_exposure.name();
}
