/** @file SpectralFunction.cxx
@brief implementation of class SpectralFunction

$Header$
*/
#include "skymaps/SpectralFunction.h"
#include "skymaps/Band.h"
#include "skymaps/SkySpectrum.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <sstream>

using namespace skymaps;


int SpectralFunction::s_n(4);
double SpectralFunction::s_e0(1000.);
double SpectralFunction::s_flux_scale(1.);
std::vector<const skymaps::SkySpectrum*> SpectralFunction::s_exposures;

namespace{

    class Default{ 
    public:
        std::string name; int npar; double param[5];
    }defaults[]=
    { 
        {"PowerLaw",      2, {-11, 2.0}},
        {"ExpCutoff",     3, {-11, 2.0, 5e3}},
        {"BrokenPowerLaw",3, {-11, 2.0, 5e3}},
        {""}
    };

} //anom namespace

void SpectralFunction::set_simpson(int n){
    s_n=n;
    if( s_n<2 || (s_n&2) !=0 || s_n>20){
        throw std::invalid_argument("SpectralFunction::set_simpson--bad Simpsons rule count: n must be even, <20");
    }
}


void SpectralFunction::set_exposures(const skymaps::SkySpectrum* front, const skymaps::SkySpectrum* back)
{
    s_exposures.clear();
    s_exposures.push_back(front);
    s_exposures.push_back(back);
}
const skymaps::SkySpectrum* SpectralFunction::exposure(int n){
    if( s_exposures.size() <2 ) {
        throw std::invalid_argument("SpectralFunction: exposures not set");
    }
    return s_exposures.at(n);
}

SpectralFunction::SpectralFunction(const std::string& name, const std::vector<double>& pars)
: m_name(name)
, m_pars(pars)
{
    if( s_exposures.size() <2 ) {
        throw std::invalid_argument("SpectralFunction: exposures not set");
    }
    setup();
}
void SpectralFunction::setup()
{
    int i(0); 
    while(1){
        if(defaults[i].name.empty()){
            std::stringstream buf;
            buf << "Did not find spectral name \"" << m_name<<"\"";;
            std::cerr << buf.str() << std::endl;
            std::cerr << "Names are: ";
            for( int j(0); ! defaults[j].name.empty(); ++j){
                std::cerr << defaults[j].name << " " ;
            }
            std::cerr << std::endl;

            throw std::invalid_argument( buf.str());
        }
        if( defaults[i].name == m_name) break;
        ++i;
    }
           
    int n(defaults[i].npar);
    m_index = i;
    
    if( m_pars.empty() ){
        // load defaults
        m_pars.reserve(n);
        std::copy(defaults[i].param, defaults[i].param+n, m_pars.begin());   
    }
    if( m_pars.size() < n ){
        throw std::invalid_argument("SpectralFunction: two few parameters set");
    }
}

double SpectralFunction::value(double e)const
{
    switch( m_index ){
        case 0: //PowerLaw:
/*          n0         log10 differential flux at e0 MeV
            gamma      (absolute value of) photon index
            return (10**n0/self.flux_scale)*(self.e0/e)**gamma
*/
            {
            double n0(m_pars[0]), gamma(m_pars[1]);
            return std::pow(10.,n0)/s_flux_scale*std::pow(s_e0/e, gamma); 
            }
        case 1: //ExpCutoff:
/*              n0         log10 differential flux at e0 MeV
                gamma      (absolute value of) photon index
                cutoff     e-folding cutoff energy (MeV)
      """
      if cutoff < 0: return 0
      return (10**n0/self.flux_scale)*(self.e0/e)**gamma*N.exp(-e/cutoff)
 */
            {
            double n0(m_pars[0]), gamma(m_pars[1]), cutoff(m_pars[2]);
            if( cutoff<0) return 0;
            return std::pow(10.,n0)/s_flux_scale*std::pow(s_e0/e, gamma)*exp(-e/cutoff);
            }
        case 2: // BrokenPowerLaw
            {
                return 0;
            }
        default:
            return 0; // temporary
    }
}


double SpectralFunction::expected(const astro::SkyDir& dir,const skymaps::Band& band)const
{
    double a(band.emin()), b(band.emax()); // range of integration
    const SkySpectrum& expose  =* exposure(band.event_class()) ; // exposure object to use

    double step( log(b/a)/s_n ), // step in log scale
           ratio( exp(step) ), // ratio of energies
           c(2.),              // initial simpsons
           e(a);              //  inital energy

    double result( a* expose(dir,a)*value(a) ); // value at low end
    for( int i = 1; i< s_n; ++i ){
        e *= ratio; // next energy
        c = 6-c;   // toggle Simpsons coefficient
        result += c * e* expose(dir,e)*value(  e); 
    }
    
    result += b*expose(dir,b)*value(b);// value at high end.
    return result*step/3.;
}

