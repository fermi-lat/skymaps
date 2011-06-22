/** @file SkyImageSpectrum.h
    @brief declare class SkyImageSpectrum

$Header$

*/
#ifndef skymaps_SkyImageSpectrum_h
#define skymaps_SkyImageSpectrum_h

#include "skymaps/SkySpectrum.h"
#include "skymaps/SkyImage.h"

#include "astro/SkyFunction.h"
#include "astro/SkyDir.h"

#include <vector>
#include <cassert>

namespace skymaps {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class SkyImageSpectrum
    @brief a SkyFunction that adapts a diffuse map. also includes extragal diffuse

*/

class SkyImageSpectrum : public skymaps::SkySpectrum {
public:
 
    /** @brief ctor to just define a constant value
        @param value [1] the constant value
        */
    SkyImageSpectrum(const std::string& filename, const std::string& extension="", bool interpolate=false);
    virtual ~SkyImageSpectrum();

    virtual double value(const astro::SkyDir& dir, double e)const;

    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    skymaps::SkyImage skyimage()const{return m_skyimage;}
    std::string filename()const{return m_filename;}
    std::string extension()const{return m_extension;}
    bool interpolate()const{return m_skyimage.interpolate();}


private:
    skymaps::SkyImage m_skyimage;
    std::string m_filename;
    std::string m_extension;
};
} // namespace skymaps
#endif

