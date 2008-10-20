/** @file Background.cxx
    @brief implement Background

$Header$
*/

#include "skymaps/Background.h"

#include <stdexcept>
#include <cmath> 
#include <iterator>
#include <sstream>
#include <stdexcept>

namespace{  // anonymous namespace for helper classes

}  // anonymous namespace
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
using namespace skymaps;

int Background::s_n(4);
void Background::set_simpson(int n){
    s_n=n;
    if( s_n<2 || (s_n&2) !=0 || s_n>20){
        throw std::invalid_argument("Background--bad Simpsons rule count: n must be even, <20");
    }
}

Background::Background(const skymaps::SkySpectrum& diffuse, double fixedexposure)
: m_diffuse(diffuse)
, m_event_type(0)
, m_fixedexposure(fixedexposure)
{
}


Background::Background(const skymaps::SkySpectrum& diffuse, 
                       const skymaps::SkySpectrum& exposuremap)
: m_diffuse(diffuse)
, m_event_type(0)
{
    m_exposures.push_back(&exposuremap);
}

Background::Background(const skymaps::SkySpectrum& diffuse, 
                       std::vector<const skymaps::SkySpectrum*> exposure_list)
: m_diffuse(diffuse)
, m_event_type(0)
{
    std::copy(exposure_list.begin(), exposure_list.end(), 
        std::back_insert_iterator<SpectrumVector>(m_exposures) );
}
Background::Background(const skymaps::SkySpectrum& diffuse, 
                       const skymaps::SkySpectrum& front, 
                       const skymaps::SkySpectrum& back)
: m_diffuse(diffuse)
, m_event_type(0)
{
    m_exposures.push_back(&front);
    m_exposures.push_back(&back);
}


Background::~Background()
{}

void Background::set_event_class(int n) const
{
    if( n>= m_exposures.size() ){
#if 0 // need a way to declare that is OK
        throw std::invalid_argument("Background:: attempt to set class beyond those available");
#else
        n=0;
#endif
    }
    m_event_type = n; // ok, since mutable

}


double Background::value(const astro::SkyDir& dir, double e)const
{
    double val( m_diffuse(dir, e) );
    if( m_exposures.size()==0){
        val *= m_fixedexposure;
    }else{
        val *= m_exposures[m_event_type]->value( dir, e);
    }
    return val;
}

///@brief integral for the energy limits, in the given direction
double Background::integral(const astro::SkyDir& dir, double a, double b)const
{
    double step( log(b/a)/s_n ), // step in log scale
           ratio( exp(step) ), // ratio of energies
           c(2.),              // initial simpsons
           e(a);              //  inital energy

    double result( a*value(dir, a) ); // value at low end
    for( int i = 1; i< s_n; ++i ){
        e *= ratio; // next energy
        c = 6-c;   // toggle Simpsons coefficient
        result += c * e* value( dir, e); 
    }
    
    result += b*value(dir, b);// value at high end.
    return result*step/3.;

}

std::string Background::name()const
{
    std::stringstream text;
    text << "Background: "+ m_diffuse.name();
    if( m_exposures.empty() ){
        text <<", Fixed:  "<< m_fixedexposure;
    }else{
        for( int i(0); i<m_exposures.size(); ++i){
            text << ", eventtype=" << i << ": " << m_exposures[i]->name();
        }
    }
    return text.str();
}
