/** @file IsotropicSpectrum.h
    @brief declare class IsotropicSpectrum

$Header$

*/
#ifndef skymaps_IsotropicPowerLaw_h
#define skymaps_IsotropicPowerLaw_h

#include "skymaps/SkySpectrum.h"
#include "skymaps/SkyImage.h"

#include "astro/SkyFunction.h"
#include "astro/SkyDir.h"

#include <map>

namespace skymaps {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class IsotropicSpectrum
    @brief a SkyFunction that adapts a diffuse map. also includes extragal diffuse

*/

class IsotropicSpectrum : public skymaps::SkySpectrum {
public:
 
    /** @brief ctor to just define an isotropic function
        @param filename text file with two columns: enrgy in MeV, flux in cm^-2 s^-1 
        */
    IsotropicSpectrum(const std::string& filename);
    virtual ~IsotropicSpectrum();


    ///@brief interpolate table 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;



private:
    double F(double)const;
    std::map<double, double> m_table;

};
} // namespace skymaps
#endif


