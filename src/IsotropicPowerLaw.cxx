/** @file IsotropicPowerLaw.cxx
    @brief implement IsotropicPowerLaw
$Header$
*/

#include "skymaps/IsotropicPowerLaw.h"
using namespace skymaps;

namespace {
    double E0(100);
};


IsotropicPowerLaw::IsotropicPowerLaw(double flux, double index)
:  m_flux(flux)
, m_index(index)
{
    setName("isotropic");
}

IsotropicPowerLaw::~IsotropicPowerLaw()
{
}

double IsotropicPowerLaw::value(const astro::SkyDir& dir, double energy)const
{
    return m_flux*(m_index-1)*pow(E0/energy, m_index-1)/energy;
}

double IsotropicPowerLaw::integral(const astro::SkyDir& /*dir*/, double a, double b)const
{
    return F(a) - F(b);
}

double IsotropicPowerLaw::F(double e)const
{
    return m_flux*pow(E0/e, m_index-1);

}
