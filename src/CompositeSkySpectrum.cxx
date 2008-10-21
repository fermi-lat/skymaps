/** @file CompositeSkySpectrum.cxx
    @brief implement class CompositeSkySpectrum

$Header$

*/

#include "skymaps/CompositeSkySpectrum.h"
#include <cassert>

using namespace skymaps;

CompositeSkySpectrum::CompositeSkySpectrum(const skymaps::SkySpectrum* component, double norm)
{
    if( component!=0 ) add(component, norm);
}

void CompositeSkySpectrum::add(const skymaps::SkySpectrum* component, double norm)
{
    push_back(std::make_pair(norm, component));
}


double CompositeSkySpectrum::value(const astro::SkyDir& dir, double energy) const
{
    double ret(0);
    std::vector< std::pair<double, const skymaps::SkySpectrum*> >::const_iterator it = begin();

    for( ; it!=end(); ++it){
        ret+= it->first * (it->second)->value(dir, energy);
    }
    return ret;
}

double CompositeSkySpectrum::band_value(const astro::SkyDir& dir, const skymaps::Band& band)const
{
    double ret(0);
    std::vector< std::pair<double, const skymaps::SkySpectrum*> >::const_iterator it = begin();

    for( ; it!=end(); ++it){
        ret+= it->first * (it->second)->band_value(dir, band);
    }
    return ret;

}



///@brief integral for the energy limits, in the given direction
double CompositeSkySpectrum::integral(const astro::SkyDir& dir, double emin, double emax)const
{
    double ret(0);
    std::vector< std::pair<double, const skymaps::SkySpectrum*> >::const_iterator it = begin();

    for( ; it!=end(); ++it){
        if( it->second==0) {

            continue;
        }
        ret+= it->first * (it->second)->integral(dir, emin, emax);
    }
    return ret;

}

std::string CompositeSkySpectrum::name()const {
    std::string ret(m_name);
    if( ret.empty()){
        std::vector< std::pair<double, const skymaps::SkySpectrum*> >::const_iterator it = begin();
        for( ; it!=end(); ++it){
            ret+= ", " + (it->second)->name();
        }
    }
    return ret;
}
