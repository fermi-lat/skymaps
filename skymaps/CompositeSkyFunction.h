/** @file CompositeSkyFunction.h
    @brief declare class CompositeSkyFunction

$Header$

*/
#ifndef skymaps_CompositeSkyFunction_h
#define skymaps_CompositeSkyFunction_h

#include "astro/SkyFunction.h"

#include <vector>
#include <algorithm> 
namespace skymaps {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class CompositeSkyFunction
    @brief a SkySpectrum that combines weighted SkySpectrum objects

*/

    class CompositeSkyFunction: public astro::SkyFunction,
        public  std::vector< std::pair<double, const astro::SkyFunction*> >{
public:
    /// @brief ctor - default, or same as an add
    /// @param component [0] pointer to component to add: if zero, the function will evaluate to zero
    CompositeSkyFunction(const astro::SkyFunction* component, double norm=1.);

    /// @brief add a new diffuse component
    /// @param component pointer to a SkySpectrum
    /// @param norm[1] the normalization factor
    void add(const astro::SkyFunction* component, double norm=1.);

    //! @brief  coordinates of a point in the sky
    //! @return value at that point
    virtual double operator()(const astro::SkyDir& bincenter)const;
    
    //! @brief evaluate average over the opening angle, with tolerance
    virtual double average(const astro::SkyDir& dir, double angle, double tolerance)const;

private:

};

}
#endif
