/** @file SkyIntegrator.h
    @brief declare SkyIntegrator class

$Header$
*/

#ifndef pointlike_SkyIntegrator_h
#define pointlike_SkyIntegrator_h

#include "astro/SkyDir.h"
#include "astro/SkyFunction.h"

#include "skymaps/SkySpectrum.h"

#include <map>


namespace skymaps {


/** @class SkyIntegator
    @brief integrate SkyFunctions and SkySpectrums over a ROI

*/

class SkyIntegrator {

public:

    SkyIntegrator(){};

    virtual ~SkyIntegrator(){};

    // Integrate functions over a radial aperture
    static double ap_int(const astro::SkyFunction& func, const astro::SkyDir& center, float radius);
    //TODO: SkySpectrum integration

    // Helper functions for radial integration
    static double average(const astro::SkyFunction& func, const astro::SkyDir& dir, double radius);
    static double ss_average(const skymaps::SkySpectrum& func, const astro::SkyDir& dir, double radius);

    // Average a SkyFunction or a SkySpectrum using a specified Healpix level (nside = 2^level)
    static double level_ave(const astro::SkyFunction& func, const astro::SkyDir& dir, double radius, int level);
    static double ss_level_ave(const skymaps::SkySpectrum& func, const astro::SkyDir& dir, double radius, int level);

    // Do the averaging using RING scheme -- purportedly more efficient, also allows better matching of structure and Healpix resolution
    // But, testing shows it's ever so slightly slower?
    static double nside_ave(const astro::SkyFunction& func, const astro::SkyDir& dir, double radius, int nside);
    static double ss_nside_ave(const skymaps::SkySpectrum& func, const astro::SkyDir& dir, double radius, int nside);

    // Integrate a function over a (presumably relatively large) pixel specified by its position and width (nside=power of 2)
    static double pix_int(const astro::SkyFunction& func, const astro::SkyDir& pix_center, int nside);
    static double level_pix_int(const astro::SkyFunction& func, const astro::SkyDir& pix_center, int level);
    
    static void set_simpson(int n){s_n=n;} //for future SkySpectrum integration
    static void set_verbose(bool q=true){verbose=q;}
    static void set_tolerance(float tol){tolerance=tol;}
    static void set_adaptive(bool q=true){adaptive=q;}

private:
    static int s_n;
    static bool verbose;
    static float tolerance;
    static bool adaptive;


};

}

#endif
