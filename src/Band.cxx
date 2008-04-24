/** @file Band.cxx
@brief implement class Band 

$Header$
*/

#include "skymaps/Band.h"

#include "healpix/Healpix.h"
using namespace skymaps;
using healpix::Healpix;
using astro::SkyDir;

astro::SkyDir Band::dir(unsigned int index)const
{
    Healpix hpx(m_nside, Healpix::RING, SkyDir::GALACTIC);
    Healpix::Pixel pix(index,hpx);

    return pix();
}

unsigned int Band::index(const astro::SkyDir& dir)const
{
    Healpix hpx(m_nside, Healpix::RING, SkyDir::GALACTIC);
    Healpix::Pixel pix(dir,hpx);

    return pix.index();
}


void Band::query_disk(const astro::SkyDir&dir, double radius, std::vector<int>& v)const
{
    Healpix hpx(m_nside, Healpix::RING, SkyDir::GALACTIC);
    hpx.query_disc( dir, radius, v);  

}
