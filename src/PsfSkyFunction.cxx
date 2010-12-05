/** @file PsfSkyFunction.cxx
@brief implement PsfSkyFunction.

$Header$
*/

#include "skymaps/PsfSkyFunction.h"
#include "skymaps/WeightedSkyDir.h"
#include "healpix/Healpix.h"
#include "healpix/HealPixel.h"

#include <map>
#include <vector>

using namespace skymaps;


PsfSkyFunction::PsfSkyFunction(const astro::SkyDir & roi_dir, double gamma, double sigma)
: m_dir(roi_dir)
, m_psf(gamma)
, m_sigma(sigma)
{}

double PsfSkyFunction::operator () (const astro::SkyDir & r)const
{
    return m_psf(r, m_dir, m_sigma);
}

std::vector<double> PsfSkyFunction::wsdl_vector_value(skymaps::WeightedSkyDirList& dirs)const
{
    std::vector<double> rvals;
    for (std::vector<skymaps::WeightedSkyDir>::const_iterator it = dirs.begin(); it != dirs.end(); ++it) {
        rvals.push_back(m_psf(*it,m_dir,m_sigma));
    }
    return rvals;
}

//copied from SkySpectrum.cxx by M. Kerr 14 January 2009
double PsfSkyFunction::average(const astro::SkyDir& dir, double angle, double tolerance)const
{
    using astro::SkyDir;
    using healpix::Healpix;

    static std::map<double, int> width_map;
    static bool map_built(false);

    int level, min_level = 6, max_level = 13;
    double result(0.0);

    // Get value for one point at center
    double previous = (*this) (dir);
    if (tolerance >= 0.5)  // If tolerance is higher than this, just return value at center.
        return previous;

    /* Build map of pixel widths by level, if not done yet.  Store twice the healpixel
       width for easy comparison. */
    if (!map_built)
    {
        width_map.clear();
        for (level = min_level; level <= max_level; ++level)
        {
            int nside(1 << level);
            int npix(12 * nside * nside);
            double width = sqrt(4 * M_PI / npix);  // Width of healpixel in radians
            width_map[2 * width] = level;
        }
        map_built = true;
    }

    // Use map to determine starting pixel level
    std::map<double, int>::iterator it = width_map.lower_bound(angle);
    if (it == width_map.end() || (it->first > angle && it != width_map.begin()))
        --it;
    level = it->second;

    // Get value for starting pixel level
    result = level_ave(dir, angle, level);

    // Iterate until result changes less than tolerance
    for(level += 1 ; fabs(result/previous -1.) > tolerance && level < max_level; ++ level)
    {
        previous = result;
        result = level_ave(dir, angle, level);
    }

    return result;

}

// Calculate average for a given level
double PsfSkyFunction::level_ave(const astro::SkyDir& dir, double angle, int level) const
{   

    int nside(1 << level);
    std::vector<int> v;
    healpix::Healpix hpx(nside, healpix::Healpix::NESTED, astro::SkyDir::GALACTIC);
    hpx.query_disc(dir, angle, v); 
    double av(0);

    for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
        healpix::HealPixel hp(*it, level);
        av += (*this) (hp());
    }

    return av/v.size();
}
