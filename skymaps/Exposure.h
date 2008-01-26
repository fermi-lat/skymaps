/** @file Exposure.h
    @brief declare class Exposure

$Header$

*/
#ifndef skymaps_Exposure_h
#define skymaps_Exposure_h

#include "skymaps/SkySpectrum.h"
#include <string>
#include "healpix/HealpixArray.h"
#include "healpix/CosineBinner.h"


namespace skymaps {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class Exposure
    @brief a SkySpectrum that represents the exposure over the sky


  
*/

class Exposure : public skymaps::SkySpectrum {
public:

    /** @brief ctor
    @param fits_file create from a FITS exposure cube
    @param tablename [Exposure]
    
    (need to connect to a discription of the effective area function or functions)

    */
    Exposure(const std::string & fits_file, const std::string& tablename="Exposure");
    ~Exposure();

    ///@brief a single energy 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    ///@param a lower limit
    ///@param b upper limit
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    virtual std::string name()const;


private:
    std::string m_filename;
    healpix::HealpixArray<healpix::CosineBinner> m_exposure;
    
};


}
#endif

