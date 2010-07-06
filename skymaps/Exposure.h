/** @file Exposure.h
    @brief declare class Exposure

$Header$

*/
#ifndef skymaps_Exposure_h
#define skymaps_Exposure_h

#include "skymaps/SkySpectrum.h"
#include "skymaps/EffectiveArea.h"
#include "skymaps/LivetimeCube.h"

//#include "healpix/Healpix.h"

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
    @param ltcube create from a livetime cube
    @param aeff an effective area object
    */
    Exposure(const skymaps::LivetimeCube& ltcube, 
        const skymaps::EffectiveArea & aeff);
    
    Exposure(const skymaps::LivetimeCube& ltcube,
             const skymaps::LivetimeCube& weighted_ltcube,
             const skymaps::EffectiveArea & aeff);

    virtual ~Exposure();

    ///@brief a single energy 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    ///@param a lower limit
    ///@param b upper limit
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    ///@brief return differential (in cosine incidence) value for PSF averaging
    double diff_value(const astro::SkyDir& dir, double e, double costh) const;

    virtual std::string name()const;

    /// @brief access to gti from the livetime cube
    const skymaps::Gti& gti()const{return m_ltcube.gti();}

    std::vector<double> vector_value(const astro::SkyDir& dir, std::vector<double>& energies)const;

    /// @brief set the number of points used for Simpson's rule integration
    static void set_simpson(int n);

    /// @brief set the cosine of the cutoff angle (default = 0.4 <=> 66.4 degrees)
    static void set_cutoff(double cutoff);

    const healpix::Healpix healpix() const {return m_ltcube.data().healpix();}

private:
    const skymaps::LivetimeCube& m_ltcube;
    const skymaps::EffectiveArea& m_aeff;
    const skymaps::LivetimeCube* m_weighted_ltcube;

    static double s_cutoff;
    static int s_n;
};

class ExposureMap : public astro::SkyFunction {
public:
    ExposureMap(const Exposure& exp, double energy): m_exp(exp), m_energy(energy){}

    double operator()(const astro::SkyDir& sdir)const{
        return m_exp.value(sdir,m_energy);
    }
private:
    const Exposure& m_exp;
    double m_energy;
};

class CacheExposureMap : public astro::SkyFunction {
public:
    CacheExposureMap(const Exposure& exp, double energy): m_exp(exp), m_energy(energy), m_hp(exp.healpix()) {}
  
    double operator()(const astro::SkyDir& sdir)const{
        healpix::Healpix::Pixel pix = m_hp.pixel(sdir);
        int index(pix.index());
        std::map<int,double>::const_iterator it = m_cache.find(index);
        if( it != m_cache.end() ){
            return (*it).second;
        }
        double val (m_exp.value(sdir, m_energy));
        m_cache[index] = val;
        return val;
    }
private:
    const Exposure& m_exp;
    double m_energy;
    mutable std::map<int,double> m_cache;
    const healpix::Healpix m_hp;
};

}
#endif

