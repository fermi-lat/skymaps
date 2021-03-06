/** @file Background.h
    @brief declare class Background

$Header$

*/
#ifndef skymaps_Background_h
#define skymaps_Background_h

#include "skymaps/SkySpectrum.h"
#include "skymaps/Band.h"
#include "skymaps/WeightedSkyDir.h"
#include <vector>
#include <algorithm> // for pair, make_pair


namespace skymaps {
// forward declarations for implementation
    class DiffuseFunction;
    class EffectiveArea;
    class LivetimeCube;
    class IsotropicSpectrum;
    class CompositeSkySpectrum;
    class Exposure;
    class CacheExposureMap;

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
    @param fixedexposure constant exposure factor to use (3e10 or so for a year)
    */
    Background(const skymaps::SkySpectrum& diffuse, double fixedexposure);

    /** @brief ctor
    @param diffuse  a SkySpectrum representing the diffuse background flux
    @param exposuremap a SkySpectrum corresponding to the exposure
    */
    Background(const skymaps::SkySpectrum& diffuse, 
        const skymaps::SkySpectrum& exposuremap);

    /** @brief ctor
    @param diffuse  a SkySpectrum representing the diffuse background flux
    @param exposure_list a vector of pointers to SkySpectrum objects corresponding to the exposure,
                    indexed accordign to the event type (usually front/back)
    */
    Background(const skymaps::SkySpectrum& diffuse, 
                       std::vector<const skymaps::SkySpectrum*> exposure_list);
    /** @brief ctor
    @param diffuse  a SkySpectrum representing the diffuse background flux
    @param front a SkySpectrum corresponding to the first event class, probably front
    @param back  a SkySpectrum corresponding to the first event class, probably back
    */
    Background(const skymaps::SkySpectrum& diffuse, 
                       const skymaps::SkySpectrum& front, 
                       const skymaps::SkySpectrum& back);

    /** @brief new ctor
        @param irfname  name of an IRF component for the effective area
        @param livetimefile a livetime cube fits file
        @param galactic a fits file with a galactic diffuse model
        @param isotropic [1.5e-5, 2.2]

    */
    Background(const std::string& irfname, 
        const std::string& livetimefile,
        const std::string& galactic, 
        const std::string& isotropic);

    ~Background();

    ///@brief a single energy 
    ///@param dir direction
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief implement integral over an energy band, selecting the event class
    virtual double band_value(const astro::SkyDir& dir, const skymaps::Band& band)const;

    ///@brief integral for the energy limits, in the given direction
    ///@param dir direction
    ///@param a lower limit
    ///@param b upper limit
    /// The integral is evaluated with extended Simpson's rule, in the log space.
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    ///@brief implement a SkyFunction-like interface
    double operator() (const astro::SkyDir& dir)const;

    ///@brief set parameters for SkyFunction-like interface
    void set_skyfun(int conversion_type, double energy)const;

    ///@brief set the event type (index into the exposures arrray) to use
    /// If not valid, silently use fixed or first entry
    /// note that this sets the event type, which is mutable
    void set_event_class(unsigned int n)const;

    ///#brief set the energy for the SkyFunction interface
    // Here for backwards compatibility
    void setEnergy(double e)const;

    ///@brief number of exposure objects
    int exposures()const{ return m_exposures.size();}

    std::string name()const;

    /// @brief set the number of points used for Simpson's rule integration
    static void set_simpson(int n);

    std::vector<double> wsdl_vector_value (skymaps::WeightedSkyDirList & wsdl) const;

    /// @brief get the rotated lons and lats for a list of positions
    static void rot_grid (std::vector<double>& rlons, std::vector<double>& rlats,
        const std::vector<double>& lons, const std::vector<double>&lats, 
        const astro::SkyDir& roi_center);

    /// @brief fill a grid on the rotated grid of lons and lats
    static void val_grid (std::vector<double>& rvals,
        const std::vector<double>& lons, const std::vector<double>&lats, 
        const astro::SkyDir& roi_center, const astro::SkyFunction& myfunc);

    void grid_values(std::vector<double>& rvals, const std::vector<double>& lons, const std::vector<double>&lats, const astro::SkyDir::CoordSystem coordsys) const;

private:
    const skymaps::SkySpectrum* m_diffuse;
    typedef std::vector<const skymaps::SkySpectrum*> SpectrumVector;
    SpectrumVector m_exposures;
    mutable int m_event_type;
    double m_fixedexposure; ///< fixed exposure to use if not a map
    static int s_n; ///< simpsons rule 
    mutable double m_energy;

    // new implementation
    const EffectiveArea* m_aeff_front;
    const EffectiveArea* m_aeff_back;
    const LivetimeCube*  m_ltcube;
    const Exposure* m_front;
    const Exposure* m_back;
    const DiffuseFunction* m_galaxy;
    const IsotropicSpectrum* m_isotropic;
    CompositeSkySpectrum* m_total_diffuse;

    // cache stuff
    mutable CacheExposureMap* m_exposure_map;
    
};


}
#endif

