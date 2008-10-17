/** @file Convolution.h
@brief Convolves healpix maps with a skyfunction

@author M. Roth 

$Header$
*/
#ifndef skymaps_Convolution_h
#define skymaps_Convolution_h

#include "skymaps/SkySpectrum.h"
#include "healpix/Map.h"
#include <map>

namespace skymaps {

    class Convolution : public skymaps::SkySpectrum, public std::map<int,healpix::Map<double> >{

       
    public:
        /*
        * Convolution of a sky map with arbitrary sky function
        * @param sf Map as a sky function
        * @param ker convolution function
        * @param level resolution of healpix map, level = log2(nside)
        */
        Convolution(const skymaps::SkySpectrum &sf, const skymaps::SkySpectrum &ker, double energy, int level);

        /*
        * Convolution of a sky map with the LAT psf a particular energy band (assumes azithmuthal symmetry of psf)
        * @param sf SkyFunction
        * @param energy determines shape of PSF (from IRF)
        * @param level max resolution (healpix bins) 6:fast-coarse -> 11:slow-fine
        */
        Convolution(const skymaps::SkySpectrum& ss, double energy, int level);

        void createConv(const skymaps::SkySpectrum&sf, const skymaps::SkySpectrum& ker, double energy);
        
        void createConv(const skymaps::SkySpectrum&sf, double energy);

        int layer(double e, bool front) const;

        virtual double value(const astro::SkyDir& dir)const;
        ///@brief interpolate table 
        ///@param e energy in MeV
        virtual double value(const astro::SkyDir& dir, double e)const;

        ///@brief integral for the energy limits, in the given direction
        virtual double integral(const astro::SkyDir& dir, double a, double b)const;

        virtual std::string name()const;//{return std::string("Convolution");}

    private:
        int m_level;
        const skymaps::SkySpectrum* m_sf;
        const skymaps::SkySpectrum* m_ker;
    };


}
#endif
