/** @file Exposure.h
    @brief declare class Exposure

$Header$

*/
#ifndef skymaps_Exposure_h
#define skymaps_Exposure_h

#include "skymaps/SkySpectrum.h"
#include "skymaps/EffectiveArea.h"
#include "skymaps/LivetimeCube.h"

#include <string>


namespace skymaps {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class Exposure
    @brief a SkySpectrum that represents the exposure over the sky

    It combines a LivetimeCube with an EffectiveArea 

  
*/

class Exposure : public skymaps::SkySpectrum {
public:

    /** @brief ctor
    @param fits_file create from a FITS exposure cube
    @param tablename [Exposure]
    */
    Exposure(const skymaps::LivetimeCube& ltcube, 
        const skymaps::EffectiveArea & aeff);

    virtual ~Exposure();

    ///@brief a single energy 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    ///@param a lower limit
    ///@param b upper limit
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    virtual std::string name()const;

    /// @brief access to gti from the livetime cube
    const skymaps::Gti& gti()const{return m_ltcube.gti();}

    /// @brief set the number of points used for Simpson's rule integration
    static void set_simpson(int n);

private:
    const skymaps::LivetimeCube& m_ltcube;
    const skymaps::EffectiveArea& m_aeff;

    static double s_cutoff;
    static int s_n;
};

class ExposureMap : public astro::SkyFunction {
public:
    ExposureMap(const Exposure& exp, double energy): m_exp(exp), m_energy(energy){}
    double operator()(const astro::SkyDir& sdir)const{
        return m_exp(sdir, m_energy);
    }
private:
    const Exposure& m_exp;
    double m_energy;
};


}
#endif

