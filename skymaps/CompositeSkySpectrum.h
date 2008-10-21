/** @file CompositeSkySpectrum.h
    @brief declare class CompositeSkySpectrum

$Header$

*/
#ifndef skymaps_CompositeSkySpectrum_h
#define skymaps_CompositeSkySpectrum_h

#include "skymaps/SkySpectrum.h"

#include <vector>
#include <algorithm> 
namespace skymaps {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class CompositeSkySpectrum
    @brief a SkySpectrum that combines weighted SkySpectrum objects

*/

class CompositeSkySpectrum: public skymaps::SkySpectrum,
    public  std::vector< std::pair<double, const skymaps::SkySpectrum*> >{
public:
    /// @brief ctor - default, or same as an add
    CompositeSkySpectrum(const skymaps::SkySpectrum* diffuse=0, double norm=1.);

    /// @brief add a new diffuse component
    /// @param component pointer to a SkySpectrum
    /// @param norm[1] the normalization factor
    void add(const skymaps::SkySpectrum* component, double norm=1.);

    ///@brief return differential value 
    ///@param energy energy in MeV
    virtual double value(const astro::SkyDir& dir, double energy)const;

    ///@brief use a band to select interval. A subclass can use it to select subclass
    /// base class implemention integrates over band
    virtual double band_value(const astro::SkyDir& dir, const skymaps::Band& band)const;
   
    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    ///@brief name for identification
    /// default is to make a list of the names of the components
    virtual std::string name()const;

    ///@brief override name
    void setName(const std::string & name){m_name=name;}
private:
    std::string m_name;

};

}
#endif
