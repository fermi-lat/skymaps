/** @file IsotropicSpectrum.cxx
    @brief implement IsotropicSpectrum
$Header$
*/

#include "skymaps/IsotropicSpectrum.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
using namespace skymaps;




IsotropicSpectrum::IsotropicSpectrum(const std::string& filename)
{
    setName("isotropic:"+filename);
    std::ifstream input_file(filename.c_str());
    if(!input_file.is_open()) {
        std::cerr << "ERROR:  Unable to open:  " << filename << std::endl;
        throw std::invalid_argument(std::string("IsotropicSpectrum: could not open file ")+filename);
    }
    m_emin=0;
    while (!input_file.eof()){
        std::string line; std::getline(input_file, line);
        std::stringstream buf(line); 
        double e, v;
        buf >> e >> v;
        if( m_emin==0)m_emin=e; 
        m_energies.push_back(e);
        m_data.push_back(v);
        m_emax=e;
    }

}

IsotropicSpectrum::~IsotropicSpectrum()
{
}
// code adapted from DiffuseFunction implementation

size_t IsotropicSpectrum::layer(double e)const
{
    if( e< m_emin  ){
        std::stringstream error; 
        error << "Diffuse function: energy out of range: "<< e ;
        throw std::invalid_argument(error.str());
    }
    if( e>m_emax ) {
        return m_energies.size();
    }
    size_t step(0);
    for( std::vector<double>::const_iterator it( m_energies.begin()); it!=m_energies.end(); ++it, ++step){
        if( (*it) > e ) break;
    }
    return step-1;
}


double IsotropicSpectrum::value(const astro::SkyDir& , double e)const
{
    size_t l(layer(e)); 

    if( l==m_energies.size()) return m_data[l-1];

    double e1(m_energies[l]), e2(m_energies[l+1]);
    double f1( m_data[l] ), f2(m_data[l+1] );
    double alpha ( log(f1/f2)/log(e2/e1) );
    return f1*pow( e1/e, alpha);

}

double IsotropicSpectrum::integral(const astro::SkyDir& dir, double a, double b)const
{
    // estimate integral by assuming power-law from start to end
    double fa(value(dir, a)), fb(value(dir,b));
    double q ( 1. - log(fa/fb)/log(b/a) );
    return fa* a * (pow(b/a, q)-1)/q;

}

