/** @file BinnedPhoton.cxx
@brief implement class BinnedPhoton 

$Header$
*/

#include "skymaps/BinnedPhoton.h"
#include "skymaps/Band.h"

using namespace skymaps;

BinnedPhoton::BinnedPhoton(const skymaps::Band* band, const astro::SkyDir dir)
: m_band(band)
{
    m_index = m_band->index(dir);

}

astro::SkyDir BinnedPhoton::dir()const
{
    return m_band->dir(m_index);
}