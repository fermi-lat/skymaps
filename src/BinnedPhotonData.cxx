/** @file BinnedPhotonData.cxx
@brief implement class BinnedPhotonData 

$Header$
*/

#include "skymaps/BinnedPhotonData.h"
#include "astro/Photon.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include <algorithm>
#include <functional>

#include <cmath>
#include <utility>
#include <stdexcept>
#include <iomanip>

using astro::SkyDir;
using astro::Photon;

using namespace skymaps;

namespace {
    skymaps::PhotonBinner default_binner;
}
BinnedPhotonData::BinnedPhotonData()
: m_binner(default_binner)
{}

BinnedPhotonData::BinnedPhotonData(const skymaps::PhotonBinner& binner)
: m_binner(binner)
{}

        /// Create  object from a saved fits file
BinnedPhotonData::BinnedPhotonData(const std::string & inputFile, const std::string & tablename)
{
// todo: read FITS file
}

        /// add a photon to the map with the given energy and direction
void BinnedPhotonData::addPhoton(const astro::Photon& gamma)
{   
    // create a emmpty band with this photon's properties
    Band newband (m_binner(gamma));
    int key(newband);
    
    // is it already in our list?
    iterator it=std::lower_bound(begin(), end(), key, std::less<int>());

    if( key!=(*it) ){
        // no, create new entry and copy in the Band
        it = insert(it, newband);
    }

    // now add the entry
    (*it).add(gamma.dir());
    ++m_photons;

}



double BinnedPhotonData::density (const astro::SkyDir & sd) const
{
    return 0;//TODO rewrite this
}

double BinnedPhotonData::value(const astro::SkyDir& dir, double e)const
{
    return 0;//TODO rewrite this
}

double BinnedPhotonData::integral(const astro::SkyDir& dir, double a, double b)const
{
    return 0;//TODO rewrite this
}

void BinnedPhotonData::info(std::ostream& out)const
{
    int total_pixels(0), total_photons(0);
    out << " nside type   emin    emax    sigma   pixels   photons\n";

    for( const_iterator it=begin();  it!=end(); ++it)
    {
        const Band& band = *it;
        int pixels(band.size()), photons(band.photons());
        out 
            <<std::setw(6) << band.nside()
            <<std::setw(4) << band.event_class()
            <<std::setw(8) << int(band.emin()+0.5)
            <<std::setw(8) << int(band.emax()+0.5)
            <<std::setw(8) << int(band.sigma()*180/M_PI*3600+0.5)
            <<std::setw(10)<< pixels
            <<std::setw(10)<< photons 
            <<std::endl;
        total_photons += photons; total_pixels+=pixels;
    }
    out << " total"
        <<std::setw(38)<<total_pixels
        <<std::setw(10)<<total_photons << std::endl;
}

void BinnedPhotonData::write(const std::string & outputFile,
            const std::string & tablename,
            bool clobber) const
{
    throw std::runtime_error("BinnedPhotonData::write -- Not implemented"); 
}

void BinnedPhotonData::writegti(const std::string & outputFile) const
{
    m_gti.writeExtension(outputFile);
}
void BinnedPhotonData::addgti(const skymaps::Gti& other)
{
    m_gti |= other;
}

void BinnedPhotonData::operator+=(const skymaps::BinnedPhotonData& other) {
#if 0 //TODO rewrite this
    for(const_iterator it=other.begin();it!=other.end();++it){
        addPixel(it->first,it->second);
    }
#endif
    addgti(other.gti());
}

