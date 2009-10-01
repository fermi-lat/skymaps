/** @file ExposureWeighter.h
    @brief declare class ExposureWeighter

$Header$

*/
#ifndef skymaps_ExposureWeighter_h
#define skymaps_ExposureWeighter_h

#include "astro/SkyDir.h"
#include <string>


namespace skymaps {

class EffectiveArea;
class LivetimeCube;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class ExposureWeighter
@helper class to provide proper weighting of PSF over effective area and/or livetime

*/

class ExposureWeighter {
public:

    ExposureWeighter(const std::string& faeff_str, const std::string& baeff_str, const std::string& livetimefile);

    ~ExposureWeighter();

    double operator()(double c_lo, double c_hi, double e_lo, double e_hi, int event_class, astro::SkyDir& dir);

private:

    EffectiveArea* m_faeff;
    EffectiveArea* m_baeff;
    LivetimeCube* m_lt;
    bool m_uselt;
};

class TopHat {
public:
    /**
    @param c_lo lower cosine limit
    @param c_hi upper cosine limit
    */
    TopHat( double c_lo, double c_hi ) : m_clo(c_lo) , m_chi(c_hi) {}

    double operator()(double costh) const
    {
        return ( (costh > m_clo) && (costh <= m_chi) ) ? 1 : 0;
    }
    double m_clo, m_chi;
};

}
#endif

