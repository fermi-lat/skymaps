/** @file BinnedPhotonData.cxx
@brief implement class BinnedPhotonData 

$Header$
*/

#include "skymaps/BinnedPhotonData.h"
#include "astro/Photon.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

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

    m_binner.add(gamma);
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
    m_binner.info(out);
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

