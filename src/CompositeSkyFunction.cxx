/** @file CompositeSkyFunction.cxx
    @brief implement class CompositeSkyFunction

$Header$

*/

#include "skymaps/CompositeSkyFunction.h"
#include <cassert>

using namespace skymaps;

CompositeSkyFunction::CompositeSkyFunction(const astro::SkyFunction* component, double norm)
{
    if( component!=0 ) add(component, norm);
}

void CompositeSkyFunction::add(const astro::SkyFunction* component, double norm)
{
    push_back(std::make_pair(norm, component));
}


double CompositeSkyFunction::operator()(const astro::SkyDir& dir) const
{
    double ret(0);
    std::vector< std::pair<double, const astro::SkyFunction*> >::const_iterator it = begin();

    for( ; it!=end(); ++it){
        const astro::SkyFunction* sfun(it->second);
        if( sfun==0 ) return 0;
        ret+= it->first * (it->second)->operator ()(dir);
    }
    return ret;
}

double CompositeSkyFunction::average(const astro::SkyDir& dir, double angle, double tolerance)const
{
    double ret(0);
    std::vector< std::pair<double, const astro::SkyFunction*> >::const_iterator it = begin();

    for( ; it!=end(); ++it){
        const astro::SkyFunction* sfun(it->second);
        if( sfun==0 ) return 0;
        ret+= it->first * (it->second)->average(dir, angle, tolerance);
    }
    return ret;

}

