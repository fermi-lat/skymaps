/** @file DiffuseFunction.h
    @brief declare class DiffuseFunction

$Header$

*/
#ifndef skymaps_DiffuseFunction_h
#define skymaps_DiffuseFunction_h

#include "skymaps/SkySpectrum.h"
#include "skymaps/SkyImage.h"

#include "astro/SkyFunction.h"
#include "astro/SkyDir.h"

#include <vector>
#include <cassert>

namespace skymaps {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class DiffuseFunction
    @brief a SkyFunction that adapts a diffuse map. also includes extragal diffuse

*/

class DiffuseFunction : public skymaps::SkySpectrum {
public:

    /** @brief ctor that reads a FITS cube represention of the diffuse, with multple layers for
         eneries. 
         @param diffuse_cube_file Name of file. It must have an extension ENERGIES with the corresponding energy valuse
         @param energy[1000] initial energy for the SkyFunction 
         @param interpolate[true] interpolate the input map

    */
    DiffuseFunction(std::string diffuse_cube_file, double energy=1000., bool interpolate=true);
 
    virtual ~DiffuseFunction();

    double extraGal(double energy) const; ///< access to the extra galactic component

    ///@brief interpolate table 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    virtual std::string name()const{return m_name;}


    //-------------------------------------------------------
    /** @brief set vector of values for set of energy bins
        @param dir  direction 
        @param energies vector of the bin edges.
        @return result vector of the values
    */
    std::vector<double> integrals(const astro::SkyDir& dir, 
        const std::vector<double>&energies)const;

    /// @return number of layers
    size_t layers()const { return m_data.layers();}

    /// @brief access to the contained SkyImage
    const skymaps::SkyImage& image()const { return m_data;}

private:
    double level_ave(const astro::SkyDir& dir, double angle, int level) const;


    std::vector<double> m_energies; ///< list of energies
    size_t layer(double e)const;

    double energy_bin(int k) const;
    std::string m_name;
    skymaps::SkyImage m_data;
    double m_emin, m_emax;

};
} // namespace skymaps
#endif

