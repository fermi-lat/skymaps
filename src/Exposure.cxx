/** @file Exposure.cxx
    @brief implement Exposure

$Header$
*/

#include "skymaps/Exposure.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include <stdexcept>

namespace{  // anonymous namespace for helper classes
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class Aeff
@brief function class implements effective area vs costh, as adapter to EffectiveArea

(currently linear only, no area.)
*/
class Aeff {
public:
    /**
    @param energy energy to evaluate it
    @param cutoff limit for cos(theta)
    */
    Aeff( const skymaps::EffectiveArea& aeff, double energy, double cutoff)
        : m_aeff(aeff)
        , m_energy(energy)
        , m_cutoff(cutoff)
    {}

    double operator()(double costh) const
    {
        return costh<m_cutoff? 0 : m_aeff(m_energy, costh);
 //          return costh<m_cutoff? 0 : (costh-m_cutoff)/(1.-m_cutoff);
    }
    const skymaps::EffectiveArea& m_aeff;
    double m_energy;
    double m_cutoff;
};

   void writeEnergies(const std::string & filename,
                      const std::vector<double> & energies) {
      std::string ext("ENERGIES");

      tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
      fileSvc.appendTable(filename, ext);
      tip::Table * table = fileSvc.editTable(filename, ext);

      table->appendField("Energy", "1D");
      table->setNumRecords(energies.size());

      tip::Table::Iterator row = table->begin();
      tip::Table::Record & record = *row;
      
      std::vector<double>::const_iterator energy = energies.begin();
      for ( ; energy != energies.end(); ++energy, ++row) {
         record["Energy"].set(*energy);
      }

      delete table;
   }

}  // anonymous namespace
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using namespace skymaps;

double Exposure::s_cutoff(0.4);
int Exposure::s_n(4);
void Exposure::set_simpson(int n){
    s_n=n;
    if( s_n<2 || (s_n&2) !=0 || s_n>20){
        throw std::invalid_argument("Exposure--bad Simpsons rule count: n must be even, <20");
    }
}

Exposure::Exposure(const LivetimeCube& ltcube, const EffectiveArea & aeff)
: m_ltcube(ltcube)
, m_aeff(aeff)
{}

Exposure::~Exposure()
{}

double Exposure::value(const astro::SkyDir& dir, double energy)const
{
    Aeff fun( Aeff(m_aeff, energy, s_cutoff) ); 
    const healpix::CosineBinner& bins = m_ltcube.bins(dir);
    double val(  bins(fun) );
    return val;
}

/// integral for the energy limits, in the given direction 
double Exposure::integral(const astro::SkyDir& dir, double a, double b)const
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

std::string Exposure::name()const
{
    return "Exposure";
}
