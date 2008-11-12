/** @file Band.cxx
@brief implement class Band 

$Header$
*/

#include "skymaps/Band.h"

#include "healpix/Healpix.h"
#include <stdexcept>

using namespace skymaps;
using healpix::Healpix;
using astro::SkyDir;

Band::Band(int nside)
            : m_nside(nside)
            , m_event_class(0)
            , m_emin(0)
            , m_emax(0)
            , m_sigma(0)
            , m_gamma(0)
            , m_healpix(new healpix::Healpix(m_nside,Healpix::RING, SkyDir::GALACTIC))
        {}


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



#if 0
Band::Band(const Band& other)
: m_nside(1)
{
    add(other);
}
#endif


void Band::add(const astro::SkyDir& dir, int count)
{
    m_pixels[index(dir)]+=count;
}
void Band::add(int i, int count)
{
    m_pixels[i]+=count;
}
double Band::operator()(const astro::SkyDir& dir)const
{
    PixelMap::const_iterator it = m_pixels.find(index(dir));
    return it == m_pixels.end() ? 0 : it->second;
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

        PixelMap::const_iterator it2 = m_pixels.find(*it);
        if( it2 != m_pixels.end() )  {
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


int Band::query_disk(const astro::SkyDir&sdir, double radius, 
                     std::vector<std::pair<int,int> > & vec)const
{
    std::vector<int> v;
    m_healpix->query_disc( sdir, radius, v); 
    int total(0);
    // Add selectedpixels to return vector
    for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it) {

        PixelMap::const_iterator it2 = m_pixels.find(*it);
        if( it2 != m_pixels.end() )  {
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
   return m_healpix->pixelArea();
}

int Band::photons()const
{
    int count(0);
    for( PixelMap::const_iterator it=m_pixels.begin(); it!=m_pixels.end(); ++it){
        count += it->second;
    }
    return count;
}

void Band::findNeighbors(int index, std::vector<int> &neighbors)const
{
    m_healpix->findNeighbors(index, neighbors);
}
 
void Band::add(const Band& other)
{
    if( m_nside==1){
        // default ctor: making a copy
        m_nside = other.nside();
        m_healpix = new healpix::Healpix(m_nside,Healpix::RING, SkyDir::GALACTIC);
        m_sigma = other.sigma();
        m_gamma = other.gamma();
        m_emin = other.emin();
        m_emax = other.emax();
    }

    for( Band::const_iterator it= other.begin(); it!=other.end(); ++it){
        add(it->first, it->second);
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WeightedSkyDirList::WeightedSkyDirList(const Band& band, const astro::SkyDir& sdir, double radius)
: m_band(band)
{
    std::vector<std::pair<astro::SkyDir,int> > vec;
    band.query_disk(sdir, radius, vec);
    for( std::vector<std::pair<astro::SkyDir,int> >::const_iterator it=vec.begin(); it!=vec.end(); ++it){
        push_back(WeightedSkyDir(it->first, it->second));
    }
}

double WeightedSkyDirList::operator()(const astro::SkyDir& sdir)const
{
    int index(m_band.index(sdir)); // look up
    const_iterator it = this->begin();
    for(; it!=this->end(); it++){
        // this is much slower than saving the pixel indices at the start.
        if( m_band.index(*it) == index) return (*it).weight(); 
    }
    return 0;

}

