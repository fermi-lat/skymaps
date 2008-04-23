/** @file BinnedPhotonData.cxx
@brief implement class BinnedPhotonData 

$Header$
*/

#include "skymaps/BinnedPhotonData.h"
#include "skymaps/BinnedPhoton.h"
#include "astro/Photon.h"
#include "healpix/Healpix.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include <cmath>
#include <utility>
#include <stdexcept>
#include <errno.h>
#include <iomanip>

using healpix::Healpix;
using astro::SkyDir;
using astro::Photon;

using namespace skymaps;


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
    BinnedPhoton b(m_binner(gamma));
//    if( b.invalid() ) return; // check that photon was within bins
    (*this)[b.band()][b.index()]++;  // create and/or increment the bin
    ++m_photons;
}

int BinnedPhotonData::extract(const BinnedPhoton& bin, double radius,
          std::vector<std::pair<unsigned int, unsigned int> > & vec) const
{

    radius *= (M_PI / 180); // convert to radians
    int total(0);
    vec.clear();

    // Get pixels for nside that are within radius
    std::vector<int> v;
    const Band* band (bin.band());
    healpix::Healpix hpx(band->nside(), Healpix::RING, SkyDir::GALACTIC);

    hpx.query_disc(bin.dir(), radius, v);  

    // select band data to search
    std::map<const Band*, std::map<unsigned int, unsigned int> >::const_iterator
        subit = find(band);
    if( subit==end()){
        throw std::invalid_argument("BinnedPhotonData: band not found");
    }
    const std::map<unsigned int, unsigned int>& submap = subit->second;


    // Add select level pixels to return vector
    for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it)
    {

        Healpix::Pixel hp(*it, hpx);
        std::map<unsigned int, unsigned int>::const_iterator it2 = submap.find(hp.index());
        if (it2 != submap.end()) // Not in PhotonMap
        {
            int count = it2->second;
            vec.push_back(std::make_pair(it2->first, count));
            total += count;
        }
    }
    return total;
}

double BinnedPhotonData::density (const astro::SkyDir & sd) const
{
    return 0;
}

double BinnedPhotonData::value(const astro::SkyDir& dir, double e)const
{
    return 0;
}

double BinnedPhotonData::integral(const astro::SkyDir& dir, double a, double b)const
{
    return 0;
}

void BinnedPhotonData::info(std::ostream& out)const
{
    int total_pixels(0), total_photons(0);
    out << "  band    pixels    photons\n";
    for( std::map<const Band*, std::map<unsigned int, unsigned int> >::const_iterator it=begin();
        it!=end(); ++it)
    {
        const std::map<unsigned int, unsigned int>& pixel_data = it->second;
        int pixels(pixel_data.size()), photons(0);
        for( std::map<unsigned int, unsigned int>::const_iterator it2=pixel_data.begin(); 
            it2!=pixel_data.end();++it2){
            photons += it2->second;
        }

            out 
            <<std::setw(6) <<it->first->nside() 
            <<std::setw(10)<<pixels
            <<std::setw(10)<<photons << std::endl;
        total_photons += photons; total_pixels+=pixels;
    }
    out << " total"
        <<std::setw(10)<<total_pixels
        <<std::setw(10)<<total_photons << std::endl;
}

