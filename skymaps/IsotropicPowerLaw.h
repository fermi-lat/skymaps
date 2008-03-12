/** @file IsotropicPowerLaw.h
    @brief declare class IsotropicPowerLaw

$Header$

*/
#ifndef skymaps_IsotropicPowerLaw_h
#define skymaps_IsotropicPowerLaw_h

#include "skymaps/SkySpectrum.h"
#include "skymaps/SkyImage.h"

#include "astro/SkyFunction.h"
#include "astro/SkyDir.h"

#include <vector>
#include <cassert>

namespace skymaps {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class IsotropicPowerLaw
    @brief a SkyFunction that adapts a diffuse map. also includes extragal diffuse

*/

class IsotropicPowerLaw : public skymaps::SkySpectrum {
public:
 
    /** @brief ctor to just define an isotropic power law
        @param flux [1.5e-5] flux integral above 100 MeV in cm**-2 s**-1 sr**-1 
        @param index [2.1] spectral index
        */
    IsotropicPowerLaw(double flux=1.5e-5, double index=2.1);
    virtual ~IsotropicPowerLaw();


    ///@brief interpolate table 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;



private:
    double m_flux, m_index;
    double F(double e)const;

};
} // namespace skymaps
#endif

