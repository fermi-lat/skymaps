/** @file SkyImageSpectrum.cxx
    @brief implement SkyImageSpectrum
    $Header: 
*/

#include "skymaps/SkyImageSpectrum.h"
using namespace skymaps;

SkyImageSpectrum::SkyImageSpectrum(const std::string& filename, const std::string& extension, bool interpolate)
:  m_skyimage(skymaps::SkyImage(filename,extension,interpolate))
,  m_filename(filename)
,  m_extension(extension)
{
    setName("constant");
}

SkyImageSpectrum::~SkyImageSpectrum() {}

double SkyImageSpectrum::value(const astro::SkyDir& dir, double /*energy*/)const
{
    return m_skyimage(dir);
}

double SkyImageSpectrum::integral(const astro::SkyDir& dir, double a, double b)const
{
    return m_skyimage(dir)*(b-a);
}
