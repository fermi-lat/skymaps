/** @file SmoothedSkySpectrum.cxx
@brief Convolves healpix maps with a sky function

@author M. Kerr

$Header$
*/
#include "skymaps/SmoothedSkySpectrum.h"
#include "skymaps/IParams.h"
#include "skymaps/PsfSkyFunction.h"
#include "healpix/Healpix.h"

#include <vector>

using namespace skymaps;

double SmoothedSkySpectrum::pixel_scale(2*M_PI/3);
void   SmoothedSkySpectrum::set_pixel_scale(double p) {pixel_scale = p;}

double SmoothedSkySpectrum::smoothing_radius(0.95);
void   SmoothedSkySpectrum::set_smoothing_radius(double s){smoothing_radius = s;}

int    SmoothedSkySpectrum::event_class(0);
void   SmoothedSkySpectrum::set_event_class(int e){event_class = e;}

double SmoothedSkySpectrum::max_resolution(0.01); //degrees
void   SmoothedSkySpectrum::set_max_resolution(double r) {max_resolution = r;}

SmoothedSkySpectrum::SmoothedSkySpectrum (const skymaps::SkySpectrum &sf,double energy) 
: SkySpectrum(energy),
  m_sf(&sf)
{}

double SmoothedSkySpectrum::value(const astro::SkyDir &dir, double e) const{

    //PSF params
    double g( IParams::gamma(e,event_class) );
    double s( IParams::sigma(e,event_class) );
    PsfSkyFunction psf(dir,g,s);

    //Convert PSF completion radius to smoothing radius in radians; maximum 99%
    double sr(0.99);
    sr = smoothing_radius > sr? sr : smoothing_radius;
    double radius = s * sqrt( 2*g* ( pow(1 - sr,1./(1-g)) - 1) );

    //Find the appropriate pixel size for the PSF scale
    int nside(pixel_scale / s);
    int max_nside(58.6/max_resolution);
    nside=nside>max_nside?max_nside:nside;
    nside=nside<1?1:nside;

    //Get pixels within smoothing radius
    healpix::Healpix h(nside,healpix::Healpix::RING);
    std::vector<int> v;
    h.query_disc(dir,radius,v);

    //Perform PSF-weighted average over pixels
    double w(0);
    double r(0);
    double fval,wval;

    for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it) {
        healpix::Healpix::Pixel pix(*it,h);
        astro::SkyDir d = pix();
        fval = m_sf->value(d,e);
        wval = psf(d);
        r += fval*wval;
        w += wval;
    }

    return r/w;
}

double SmoothedSkySpectrum::integral(const astro::SkyDir &dir, double a, double b) const{
    return value( dir, sqrt(a*b) );
}


