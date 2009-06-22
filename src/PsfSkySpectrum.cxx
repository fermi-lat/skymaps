/** @file PsfSkyFunction.cxx
@brief implement PsfSkyFunction.

$Header$
*/

#include "skymaps/PsfSkySpectrum.h"
#include "skymaps/IParams.h"
#include "skymaps/PsfFunction.h"

using namespace skymaps;

void PsfSkySpectrum::set_event_class(int ec) {m_ec = ec;}
void PsfSkySpectrum::set_skydir(const astro::SkyDir &nd) {m_dir = nd;}
void PsfSkySpectrum::set_jacobean(bool b) {m_jac = b;}

PsfSkySpectrum::PsfSkySpectrum(const astro::SkyDir & roi_dir)
: m_dir(roi_dir),
 m_ec(0),
 m_jac(true)
{}

double PsfSkySpectrum::value(const astro::SkyDir &r, double energy) const{
    double gamma = IParams::gamma(energy,m_ec);
    double sigma = IParams::sigma(energy,m_ec);
    PsfFunction psf(gamma);
    if (m_jac) {
        return psf(r,m_dir,sigma) / (2 * M_PI * sigma * sigma);
    }
    else {
        return psf(r,m_dir,sigma);
    }
}

double PsfSkySpectrum::integral(const astro::SkyDir &dir, double a, double b) const {
    return value(dir,sqrt(a*b));
}