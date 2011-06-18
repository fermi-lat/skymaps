/** @file IsotropicConstant.h
    @brief declare class IsotropicConstant

$Header$

*/
#ifndef skymaps_IsotropicConstant_h
#define skymaps_IsotropicConstant_h

#include "skymaps/SkySpectrum.h"
#include "skymaps/SkyImage.h"

#include "astro/SkyFunction.h"
#include "astro/SkyDir.h"

#include <vector>
#include <cassert>

namespace skymaps {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class IsotropicConstant
    @brief a SkyFunction that adapts a diffuse map. also includes extragal diffuse

*/

class IsotropicConstant : public skymaps::SkySpectrum {
public:
 
    /** @brief ctor to just define a constant value
        @param value [1] the constant value
        */
    IsotropicConstant(double constant=1);
    virtual ~IsotropicConstant();

    virtual double value(const astro::SkyDir& dir, double e)const;

    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    double constant()const{return m_constant;}


private:
    double m_constant;
};
} // namespace skymaps
#endif

