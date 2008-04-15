/** @file Background.h
    @brief declare class Background

$Header$

*/
#ifndef skymaps_Background_h
#define skymaps_Background_h

#include "skymaps/SkySpectrum.h"


namespace skymaps {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class Background
    @brief a SkySpectrum that represents the background count density for observing
    any point in the sky, basically a product of two skyspectrum objects

    The value function is the product of the diffuse and exposure
  
*/

class Background : public skymaps::SkySpectrum {
public:

    /** @brief ctor
    @param diffuse  a SkySpectrum representing the diffuse background flux
    @param exposure a SkySpectrum corresponding to the exposure
    @param n [4]    number of points to use with extended Simpson's rule integral
    
    */
    Background(const skymaps::SkySpectrum& diffuse, 
        const skymaps::SkySpectrum& exposure, int n=4);
    ~Background();

    ///@brief a single energy 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    ///@param a lower limit
    ///@param b upper limit
    /// The integral is evaluated with extended Simpson's rule, in the log space.
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    std::string name()const;

private:
    const skymaps::SkySpectrum& m_diffuse;
    const skymaps::SkySpectrum& m_exposure;
    int m_n; ///< simpsons rule 
    
};


}
#endif

