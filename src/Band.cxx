/** @file Band.cxx
@brief implement class Band 

$Header$
*/

#include "skymaps/Band.h"

#include "healpix/Healpix.h"
using namespace skymaps;

astro::SkyDir Band::dir(unsigned int index)const
{
    using healpix::Healpix;
    Healpix hpx(m_nside, Healpix::RING,astro::SkyDir::GALACTIC);
    Healpix::Pixel pix(index,hpx);

    return pix();
}

unsigned int Band::index(const astro::SkyDir& dir)const
{
    using healpix::Healpix;
    Healpix hpx(m_nside, Healpix::RING,astro::SkyDir::GALACTIC);
    Healpix::Pixel pix(dir,hpx);

    return pix.index();
}
