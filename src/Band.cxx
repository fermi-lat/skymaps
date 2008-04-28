/** @file Band.cxx
@brief implement class Band 

$Header$
*/

#include "skymaps/Band.h"

#include "healpix/Healpix.h"
using namespace skymaps;
using healpix::Healpix;
using astro::SkyDir;

Band::Band(int nside, int event_class, double emin,double emax,
            double sigma, double gamma)
            : m_nside(nside)
            , m_event_class(event_class)
            , m_emin(emin)
            , m_emax(emax)
            , m_sigma(sigma)
            , m_gamma(gamma)
            , m_healpix(new healpix::Healpix(m_nside,Healpix::RING, SkyDir::GALACTIC))
        {}

void Band::add(const astro::SkyDir& dir, int count)
{
    (*this)[index(dir)]+=count;
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
    Healpix::Pixel pix(index,*m_healpix);

    return pix();
}

int Band::index(const astro::SkyDir& dir)const
{
    Healpix::Pixel pix(dir,*m_healpix);

    return pix.index();
}


int Band::query_disk(const astro::SkyDir&sdir, double radius, 
                     std::vector<std::pair<astro::SkyDir,int> > & vec)const
{
    std::vector<int> v;
    m_healpix->query_disc( sdir, radius, v); 
    int total(0);
    // Add selectedpixels to return vector
    for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it) {
        const_iterator it2 = find(*it);
        if( it2 != end() )  {
            int count = it2->second;
            vec.push_back( std::make_pair(dir(it2->first), count));
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
   return m_healpix->pixelArea();
}

int Band::photons()const
{
    int count(0);
    for( const_iterator it=begin(); it!=end(); ++it){
        count += it->second;
    }
    return count;
}
