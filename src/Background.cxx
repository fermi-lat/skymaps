/** @file Background.cxx
    @brief implement Background

$Header$
*/

#include "skymaps/Background.h"
#include "skymaps/Band.h"

#include "skymaps/DiffuseFunction.h"
#include "skymaps/CompositeSkySpectrum.h"
#include "skymaps/BinnedPhotonData.h"
#include "skymaps/EffectiveArea.h"
#include "skymaps/IsotropicSpectrum.h"
#include "skymaps/LivetimeCube.h"
#include "skymaps/Exposure.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "healpix/Healpix.h"


#include <stdexcept>
#include <cmath> 
#include <iterator>
#include <sstream>
#include <stdexcept>

namespace{  // anonymous namespace for helper classes

}  // anonymous namespace
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
using namespace skymaps;

int Background::s_n(4);
void Background::set_simpson(int n){
    if( n<2 || (n&1) !=0 || n>20){
        throw std::invalid_argument("Background--bad Simpsons rule count: n must be even, <=20");
    }
    s_n=n;
}

Background::Background(const skymaps::SkySpectrum& diffuse, double fixedexposure)
: m_diffuse(&diffuse)
, m_event_type(0)
, m_fixedexposure(fixedexposure)
, m_energy(1000)
, m_exposure_map(0)
{
    set_skyfun(m_event_type,m_energy);
}


Background::Background(const skymaps::SkySpectrum& diffuse, 
                       const skymaps::SkySpectrum& exposuremap)
: m_diffuse(&diffuse)
, m_event_type(0)
, m_energy(1000)
, m_exposure_map(0)
{
    m_exposures.push_back(&exposuremap);
    set_skyfun(m_event_type,m_energy);
}

Background::Background(const skymaps::SkySpectrum& diffuse, 
                       std::vector<const skymaps::SkySpectrum*> exposure_list)
: m_diffuse(&diffuse)
, m_event_type(0)
, m_energy(1000)
, m_exposure_map(0)
{
    std::copy(exposure_list.begin(), exposure_list.end(), 
        std::back_insert_iterator<SpectrumVector>(m_exposures) );
    set_skyfun(m_event_type,m_energy);
}
Background::Background(const skymaps::SkySpectrum& diffuse, 
                       const skymaps::SkySpectrum& front, 
                       const skymaps::SkySpectrum& back)
: m_diffuse(&diffuse)
, m_event_type(0)
, m_energy(1000)
, m_exposure_map(0)
{
    m_exposures.push_back(&front);
    m_exposures.push_back(&back);
    set_skyfun(m_event_type,m_energy);
}

Background::Background(const std::string& irfname, const std::string& livetimefile,
        const std::string& galactic, 
        const std::string& isotropic)
{
    m_aeff_front= new EffectiveArea(irfname+"_front");
    m_aeff_back = new EffectiveArea(irfname+"_back");
    m_ltcube=     new LivetimeCube(livetimefile);
    m_front =     new Exposure(*m_ltcube, *m_aeff_front);
    m_back =      new Exposure(*m_ltcube, *m_aeff_back);

    m_galaxy    = new DiffuseFunction(galactic);
    m_isotropic=  new IsotropicSpectrum(isotropic);
    m_total_diffuse = new CompositeSkySpectrum(m_galaxy, 1.0);
    m_total_diffuse->add(m_isotropic, 1.0);

    std::cout << "Set up background with components:"
        << "\n\tAeff:      " << irfname << "(normal 1 GeV: front="
                             << (*m_aeff_front)(1000.)
                             << ", back="<<(*m_aeff_back)(1000.)<<")"
        << "\n\tlivetime:  " << livetimefile
        << "\n\tgalactic:  " << galactic
        << "\n\tisotropic: " << isotropic
        << std::endl;

    m_diffuse = m_total_diffuse;
    m_exposures.push_back(m_front);
    m_exposures.push_back(m_back);
    m_energy = 1000;
    m_event_type=0;
    m_exposure_map = 0;
    set_skyfun(m_event_type,m_energy);
}


Background::~Background()
{
    if (m_exposure_map != 0) {
        delete m_exposure_map;
    }
}

void Background::set_event_class(unsigned int n) const
{
    if( n>= m_exposures.size() ){
#if 0 // need a way to declare that is OK
        throw std::invalid_argument("Background:: attempt to set class beyond those available");
#else
        n=0;
#endif
    }
    set_skyfun(n,m_energy); // make new cache

}

void Background::setEnergy(double e)const
{
    set_skyfun(m_event_type,e); // make new cache
}

double Background::value(const astro::SkyDir& dir, double e)const
{
    double val( (*m_diffuse)(dir, e) );
    if( m_exposures.size()==0){
        val *= m_fixedexposure;
    }else{
        val *= m_exposures[m_event_type]->value( dir, e);
    }
    return val;
}

double Background::operator ()(const astro::SkyDir& dir)const
{
    return (*m_diffuse)(dir,m_energy)*(*m_exposure_map)(dir);
}

void Background::set_skyfun(int conversion_type,double energy) const
{
    m_energy = energy;
    m_event_type = conversion_type;
    if (m_exposure_map != 0) {
        delete m_exposure_map;
    }
    const skymaps::Exposure& my_exp = reinterpret_cast<const skymaps::Exposure&> (*(m_exposures[m_event_type]));
    m_exposure_map = new skymaps::CacheExposureMap(my_exp,m_energy);
}

double Background::band_value(const astro::SkyDir& dir, const skymaps::Band& band)const
{
    set_event_class(band.event_class());
    return integral(dir, band.emin(), band.emax());
}




///@brief integral for the energy limits, in the given direction
double Background::integral(const astro::SkyDir& dir, double a, double b)const
{
    double step( log(b/a)/s_n ), // step in log scale
           ratio( exp(step) ), // ratio of energies
           c(2.),              // initial simpsons
           e(a);              //  inital energy

    double result( a*value(dir, a) ); // value at low end
    for( int i = 1; i< s_n; ++i ){
        e *= ratio; // next energy
        c = 6-c;   // toggle Simpsons coefficient
        result += c * e* value( dir, e); 
    }
    
    result += b*value(dir, b);// value at high end.
    return result*step/3.;

}

std::string Background::name()const
{
    std::stringstream text;
    text << "Background: "+ m_diffuse->name();
    if( m_exposures.empty() ){
        text <<", Fixed:  "<< m_fixedexposure;
    }else{
        for(unsigned int i(0); i<m_exposures.size(); ++i){
            text << ", eventtype=" << i << ": " << m_exposures[i]->name();
        }
    }
    return text.str();
}

std::vector<double> Background::wsdl_vector_value(skymaps::WeightedSkyDirList& dirs)const
{
    std::vector<double> rvals;
    for (std::vector<skymaps::WeightedSkyDir>::const_iterator it = dirs.begin(); it != dirs.end(); ++it) {
        rvals.push_back(this->operator()(*it) );
    }
    return rvals;
}

// This will rotate a grid of lons/lats from the ROI frame to the
// flat grid frame on the equator.  (Presumably the lons/lats
// are in the ROI.)  They are in the Galactic coordinate system.
void Background::rot_grid (std::vector<double>& rlons, std::vector<double>& rlats,
                           const std::vector<double>& lons, const std::vector<double>&lats, 
                           const astro::SkyDir& roi_center)
{
    rlons.clear(); rlats.clear();
    rlons.reserve(lons.size()); rlats.reserve(lats.size());
    astro::SkyDir rot_axis = astro::SkyDir(roi_center.l()+90,0);
    Hep3Vector rot_axisv(rot_axis.dir());
    astro::SkyDir sd;
    double rot_extent = roi_center.b()*(M_PI/180);
    for (unsigned int i=0; i < lons.size(); ++i) {
        sd = astro::SkyDir(lons[i],lats[i]);
        sd().rotate(rot_extent,rot_axisv);
        rlons.push_back(sd.ra());
        rlats.push_back(sd.dec());
    }
}

// This will take the flat grid at the equator, rotate it to the
// ROI frame, and evaluate a SkyFunction over it.
void Background::val_grid (std::vector<double>& rvals,
                           const std::vector<double>& lons, const std::vector<double>&lats, 
                           const astro::SkyDir& roi_center, const astro::SkyFunction& myfunc)
{
    rvals.clear();
    rvals.reserve(lons.size()*lats.size());
    astro::SkyDir rot_axis = astro::SkyDir(roi_center.l()-90,0);
    Hep3Vector rot_axisv(rot_axis.dir());
    astro::SkyDir sd;
    double rot_extent = roi_center.b()*(M_PI/180);
    for (std::vector<double>::const_iterator it_lon = lons.begin(); it_lon != lons.end(); ++it_lon){
        for (std::vector<double>::const_iterator it_lat = lats.begin(); it_lat != lats.end(); ++it_lat){
            sd = astro::SkyDir(*it_lon,*it_lat);
            sd().rotate(rot_extent,rot_axisv);
            rvals.push_back(myfunc(astro::SkyDir(sd.ra(),sd.dec(),astro::SkyDir::GALACTIC)));
        }
    }
}
void Background::grid_values(std::vector<double>& rvals, const std::vector<double>& lons, const std::vector<double>&lats, const astro::SkyDir::CoordSystem coordsys) const
{
    rvals.clear();
    rvals.reserve(lons.size()*lats.size());
    for (std::vector<double>::const_iterator it_lon = lons.begin(); it_lon != lons.end(); ++it_lon){
        for (std::vector<double>::const_iterator it_lat = lats.begin(); it_lat != lats.end(); ++it_lat){
            rvals.push_back(this->operator ()(astro::SkyDir((*it_lon),(*it_lat),coordsys)));
        }
    }
}