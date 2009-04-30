/** @file SkyIntegator.cxx

$Header$

*/

#include "skymaps/SkyIntegrator.h"
#include "healpix/Healpix.h"
#include "healpix/HealPixel.h"

#include <map>
#include <deque>
#include <iostream>

using namespace skymaps;
using astro::SkyDir;

bool  SkyIntegrator::verbose(false);
int   SkyIntegrator::s_n(20); //make sufficient to integrate over a broad energy band
float SkyIntegrator::tolerance(0.05); //default 5% tolerance when testing for convergence

void SkyIntegrator::set_simpson(int n){s_n=n;}
void SkyIntegrator::set_verbose(bool q){verbose=q;}
void SkyIntegrator::set_tolerance(float tol){tolerance=tol;}

double SkyIntegrator::ap_int(const astro::SkyFunction& func, const astro::SkyDir& center, float radius) {
    return 3.141592654*radius*radius*average(func,center,radius);
}

double SkyIntegrator::average(const astro::SkyFunction& func, const astro::SkyDir& dir, double radius)
{
    using astro::SkyDir;
    using healpix::Healpix;

    static std::map<double, int> width_map;
    static bool map_built(false);

    int level, min_level = 6, max_level = 13;
    double result(0.0);

    // Get value for one point at center
    double previous = func (dir);
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
    // This just finds a starting level that has half the scale of the radius -- would it be better to use an estimated scale?
    // e.g, integrating a PSF over a region of interest large to the PSF scale
    std::map<double, int>::iterator it = width_map.lower_bound(radius);
    if (it == width_map.end() || (it->first > radius && it != width_map.begin()))
        --it;
    level = it->second;

    // Get value for starting pixel level
    result = level_ave(func,dir, radius, level);

    // Iterate until result changes less than tolerance
    for(level += 1 ; fabs(result/previous -1.) > tolerance && level < max_level; ++ level)
    {
        previous = result;
        result = level_ave(func,dir, radius, level);
    }
    return result;
}
// Calculate average for a given level
double SkyIntegrator::level_ave(const astro::SkyFunction& func, const astro::SkyDir& dir, double radius, int level)
{   

    int nside(1 << level);
    std::vector<int> v;
    healpix::Healpix hpx(nside, healpix::Healpix::NESTED, astro::SkyDir::GALACTIC);
    hpx.query_disc(dir, radius, v); 
    double av(0);

    for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
        healpix::HealPixel hp(*it, level);
        av += func(hp());
    }

    return av/v.size();
}
double SkyIntegrator::level_pix_int(const astro::SkyFunction& func, const astro::SkyDir& pix_center, int level) {
    return pix_int(func,pix_center,1<<level);
}

double SkyIntegrator::pix_int(const astro::SkyFunction& func, const astro::SkyDir& pix_center, int nside){
    //basic idea here is to use the nested structure to quickly estimate the integral; go to the next level
    //(nside * 2) and calculate the integral, compare to previous value, return when meet tolerance

    if ( (nside & (nside - 1)) != 0) {
        std::cout << "\nnside must be a power of 2! returning coarse average\n";
        healpix::Healpix hpx(nside, healpix::Healpix::RING, astro::SkyDir::GALACTIC);
	    return func(pix_center)*hpx.pixelArea();
    }    

    healpix::Healpix hpx(nside, healpix::Healpix::NESTED, astro::SkyDir::GALACTIC);
	healpix::Healpix::Pixel pix(pix_center,hpx);
    int index (pix.index()); //Healpix index for integral region
    double result = func(pix_center)*hpx.pixelArea(); //initial estimate for integral

    std::deque<int> indices(1,index);
    //indices.push_back(index);

    float eps = 0.5; //initial value for difference, to compare to tolerance
    int counter(1);
	int max_iter(5);

    while (eps > tolerance){

        int nnside = nside << counter;
        if ( (counter > max_iter) || nnside >= 8192) {break;}

        // replace existing indices with ones corresponding to twice finer gridding
		int current_size (indices.size());
        for (int i = 0; i < current_size; i++){
            int ind = indices.front();
			indices.pop_front();
            for (int j = 0; j < 4; j++) {
                indices.push_back( (ind << 2) + j);
            }
        }
        
        double new_result(0);
        healpix::Healpix new_hpx( nnside, healpix::Healpix::NESTED, astro::SkyDir::GALACTIC);       
        for (std::deque<int>::const_iterator it = indices.begin(); it != indices.end(); ++it){
            healpix::Healpix::Pixel hp(*it,new_hpx);
            new_result += func(hp());
		}
        new_result *= (float(hpx.pixelArea())/indices.size());
        eps = abs( result - new_result )/result;
        result = new_result;
        counter++;        
    }
    return result;
}