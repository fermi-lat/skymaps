/** @file Band.cxx
@brief implement class Band 

$Header:
*/

#include "skymaps/WeightedSkyDir.h"

#include <stdexcept>

using namespace skymaps;
using astro::SkyDir;

void BaseWeightedSkyDirList::arclength(const astro::SkyDir& sdir, std::vector<double>& output)const
{
    output.clear();
    const_iterator it = this->begin();
    for (; it!=this->end(); it++) {
        output.push_back( sdir.difference(*it) );
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WeightedSkyDirList::WeightedSkyDirList(const Band& band, const astro::SkyDir& sdir, double radius, bool not_empty)
: m_band(band)
{
    std::vector<std::pair<astro::SkyDir,int> > vec;
    m_counts = band.query_disk(sdir, radius, vec, not_empty);
    m_pix = band.cache_pix(); // save this to avoid a possibly-expensive query_disc call in future
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



