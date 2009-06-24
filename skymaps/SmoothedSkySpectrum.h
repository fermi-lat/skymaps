/** @file SmoothedSkySpectrum.h
@brief smooth a SkySpectrum with an arbitrary kernel (also a SkySpectrum)

@author M. Kerr 

$Header$
*/
#ifndef skymaps_SmoothedSkySpectrum_h
#define skymaps_SmoothedSkySpectrum_h

#include "skymaps/SkySpectrum.h"
#include "astro/SkyDir.h"


namespace skymaps {

    class SmoothedSkySpectrum : public skymaps::SkySpectrum {
       
    public:
        /*
        * Smoothing/convolution of a sky map with the PSF (extend to arbitrary in a later version)
        * @param sf function to be smoothed
        */
        SmoothedSkySpectrum(const skymaps::SkySpectrum &sf,double energy=1000);

        double value(const astro::SkyDir &dir, double e) const;
        double integral(const astro::SkyDir &dir, double a, double b) const;
        static void set_pixel_scale(double p);
        static void set_smoothing_radius(double s);
        static void set_event_class(int e);
        static void set_max_resolution(double r);

    private:
        static double pixel_scale;
        static double smoothing_radius;
        static double max_resolution;
        static int    event_class;
        const skymaps::SkySpectrum* m_sf;
    };


}
#endif
