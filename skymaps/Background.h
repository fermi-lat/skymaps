/** @file Background.h
    @brief declare class Background

$Header$

*/
#ifndef skymaps_Background_h
#define skymaps_Background_h

#include "skymaps/SkySpectrum.h"
#include <vector>


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
    Background::Background(const skymaps::SkySpectrum& diffuse, 
                       std::vector<const skymaps::SkySpectrum*> exposure_list);
    /** @brief ctor
    @param diffuse  a SkySpectrum representing the diffuse background flux
    @param front a SkySpectrum corresponding to the first event class, probably front
    @param back  a SkySpectrum corresponding to the first event class, probably back
    */
    Background(const skymaps::SkySpectrum& diffuse, 
                       const skymaps::SkySpectrum& front, 
                       const skymaps::SkySpectrum& back);
    ~Background();

    ///@brief a single energy 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    ///@param a lower limit
    ///@param b upper limit
    /// The integral is evaluated with extended Simpson's rule, in the log space.
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    ///@brief set the event type (index into the exposures arrray) to use
    /// If not valid, silently use fixed or first entry
    /// note that this sets the event type, which is mutable
    void set_event_class(int n)const;


    std::string name()const;

    /// @brief set the number of points used for Simpson's rule integration
    static void set_simpson(int n);

private:
    const skymaps::SkySpectrum& m_diffuse;
    typedef std::vector<const skymaps::SkySpectrum*> SpectrumVector;
    SpectrumVector m_exposures;
    mutable int m_event_type;
    double m_fixedexposure; ///< fixed exposure to use if not a map
    static int s_n; ///< simpsons rule 
    
};


}
#endif

