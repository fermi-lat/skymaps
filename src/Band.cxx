/** @file Band.cxx
@brief implement class Band 

$Header$
*/

#include "skymaps/Band.h"

#include "healpix/Healpix.h"
using namespace skymaps;
using healpix::Healpix;
using astro::SkyDir;

void Band::add(const astro::SkyDir& dir)
{
    (*this)[index(dir)]++;
}
void Band::add(int i, int count)
{
    (*this)[i]+=count;
}
double Band::operator()(const astro::SkyDir& dir)const
{
    const_iterator it = find(index(dir));
    return it == end() ? 0 : it->second;
}


astro::SkyDir Band::dir( int index)const
{
    Healpix hpx(m_nside, Healpix::RING, SkyDir::GALACTIC);
    Healpix::Pixel pix(index,hpx);

    return pix();
}

int Band::index(const astro::SkyDir& dir)const
{
    Healpix hpx(m_nside, Healpix::RING, SkyDir::GALACTIC);
    Healpix::Pixel pix(dir,hpx);

    return pix.index();
}


int Band::query_disk(const astro::SkyDir&sdir, double radius, 
                      std::vector<std::pair<int,int> > & vec)const
{
    Healpix hpx(m_nside, Healpix::RING, SkyDir::GALACTIC);
    std::vector<int> v;
    hpx.query_disc( sdir, radius, v); 
    int total(0);
    // Add select level pixels to return vector
    for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it) {
        const_iterator it2 = find(*it);
        if( it2 != end() )  {
            int count = it2->second;
            vec.push_back( std::make_pair(it2->first, count));
            total += count;
#if 0 // option to return empty pixels within range
        }else{
            vec.push_back(std::make_pair(*it, 0));
#endif
        }
    }
    return total;
}

double Band::pixelArea()const
{
   Healpix hpx(m_nside, Healpix::RING, SkyDir::GALACTIC);
   return hpx.pixelArea();
}

int Band::photons()const
{
    int count(0);
    for( const_iterator it=begin(); it!=end(); ++it){
        count += it->second;
    }
    return count;
}
