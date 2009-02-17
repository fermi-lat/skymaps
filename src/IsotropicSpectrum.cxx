/** @file IsotropicSpectrum.cxx
    @brief implement IsotropicSpectrum
$Header$
*/

#include "skymaps/IsotropicSpectrum.h"
using namespace skymaps;




IsotropicSpectrum::IsotropicSpectrum(const std::string& filename)
{
    setName("isotropic");
}

IsotropicSpectrum::~IsotropicSpectrum()
{
}

double IsotropicSpectrum::value(const astro::SkyDir& dir, double energy)const
{
    return 0; //TODO
}

double IsotropicSpectrum::integral(const astro::SkyDir& /*dir*/, double a, double b)const
{
    return F(a) - F(b);
}

double IsotropicSpectrum::F(double e)const
{
    return 0; //TODO

}
