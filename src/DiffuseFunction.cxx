/** @file DiffuseFunction.cxx
    @brief implement DiffuseFunction
$Header$
*/

#include "skymaps/DiffuseFunction.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "healpix/Healpix.h"
#include "healpix/HealPixel.h"
#include "tip/Table.h"
#include "tip/IFileSvc.h"
#include <cmath>
#include <map>
#include <stdexcept>

using namespace skymaps;


DiffuseFunction::DiffuseFunction(std::string diffuse_cube_file, double energy, bool interpolate)
: SkySpectrum(energy)
, m_name(diffuse_cube_file)
, m_data(diffuse_cube_file, "", interpolate)
{
    // expect to find a table with the energies to correspond with the layers
    try {
        using namespace tip;
        Table* table = IFileSvc::instance().editTable(diffuse_cube_file, "ENERGIES");

        for (Table::Iterator itor = table->begin(); itor != table->end(); ++itor) {

            double e = (*itor)["energy"].get();
            m_energies.push_back(e);

        }
    }catch(const std::exception& ){
        throw std::invalid_argument("DiffuseFunction: the diffuse cube is not valid");
    }
    if( layers() != m_energies.size() ){
        throw std::invalid_argument("DiffuseFunction: the diffuse cube is not valid");
    }
    m_emin = m_energies[0];
    m_emax = m_energies.back();
    std::cout << "DiffuseFunction: read file "<< diffuse_cube_file <<", with " 
        << layers() << " energies from " << m_emin << " to " << m_emax << std::endl;
    setEnergy(energy);
}

DiffuseFunction::~DiffuseFunction()
{
}

double DiffuseFunction::energy_bin(int k)const
{
    return m_energies[k];
}

double DiffuseFunction::extraGal(double energy)const
{
    static double E0(100), flux(1.5E-5), alpha(1.1);
    return flux*alpha*pow(E0/energy, alpha)/energy;
}

size_t DiffuseFunction::layer(double e)const
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


double DiffuseFunction::value(const astro::SkyDir& dir, double e)const
{
    size_t l(layer(e)); 
    
    if( l==m_energies.size()) return m_data.pixelValue(dir,l-1);

    double e1(m_energies[l]), e2(m_energies[l+1]);
    double f1( m_data.pixelValue(dir,l) ), f2(m_data.pixelValue(dir,l+1) );
    double alpha ( log(f1/f2)/log(e2/e1) );
    return f1*pow( e1/e, alpha);

}

double DiffuseFunction::integral(const astro::SkyDir& dir, double a, double b)const
{
    // estimate integral by assuming power-law from start to end
    double fa(value(dir, a)), fb(value(dir,b));
    double q ( 1. - log(fa/fb)/log(b/a) );
    return fa* a * (pow(b/a, q)-1)/q;
}

std::vector<double> DiffuseFunction::integrals(const astro::SkyDir& dir, const std::vector<double>&energies)const
{
    std::vector<double> result;
    static double infinity (300000);

    for( std::vector<double>::const_iterator it = energies.begin(); it!=energies.end(); ++it){
        std::vector<double>::const_iterator next(it+1);
        double a = *it;
        double b = next!=energies.end()? *next : infinity;
        result.push_back(integral(dir,a,b));
    }
    return result;
}
