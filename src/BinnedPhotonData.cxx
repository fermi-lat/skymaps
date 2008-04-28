/** @file BinnedPhotonData.cxx
@brief implement class BinnedPhotonData 

$Header$
*/

#include "skymaps/BinnedPhotonData.h"

#ifndef OLD
#include "healpix/HealPixel.h"
#endif
#include "astro/Photon.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include <algorithm>
#include <functional>

#include <cmath>
#include <utility>
#include <stdexcept>
#include <iomanip>

using astro::SkyDir;
using astro::Photon;

using namespace skymaps;

namespace {
    skymaps::PhotonBinner default_binner;
}
BinnedPhotonData::BinnedPhotonData()
: m_binner(default_binner)
{}

BinnedPhotonData::BinnedPhotonData(const skymaps::PhotonBinner& binner)
: m_binner(binner)
{}


BinnedPhotonData::BinnedPhotonData(const std::string & inputFile, const std::string & tablename)
: m_binner(default_binner) // should change, or make flexible
, m_photons(0)
{
    if( tablename != "PHOTONMAP"){
        /// @TODO: implement native scheme
        throw std::invalid_argument("BinnedPhotonData: table "+tablename+" is not supported");
    }
    // read in and convert old-style PHOTONMAP
    const tip::Table & table=*tip::IFileSvc::instance().readTable(inputFile, tablename);
    const tip::Header& hdr = table.getHeader();
    double eratio;
    int stored_photons(0), stored_pixels(0);


    double m_emin, m_logeratio;
    int m_levels, m_minlevel, m_pixels;
    using healpix::HealPixel;

    // Guard against headers not being found in fits file.  Set to default on error

    try	{hdr["EMIN"].get(m_emin);} catch (const std::exception& ) {m_emin = 100.;}
    try
    {
        hdr["ERATIO"].get(eratio);
        m_logeratio = log(eratio);
    }
    catch (const std::exception& ) {m_logeratio = log(2.35);}
    try	{hdr["LEVELS"].get(m_levels);} catch (const std::exception& ) {m_levels = 8;}
    try	{hdr["MINLEVEL"].get(m_minlevel);} catch (const std::exception& ) {m_minlevel = 6;}
    try
    {
        hdr["PHOTONS"].get(stored_photons);
        hdr["PIXELS"].get(stored_pixels);
    }
    catch (const std::exception& ) {}
    tip::Table::ConstIterator itor = table.begin();
    std::cout << "Creating BinnedPhotonData from file " << inputFile << ", table " << tablename << std::endl;

    for(tip::Table::ConstIterator itor = table.begin(); itor != table.end(); ++itor)
    {
        long level, index, count;
        (*itor)["LEVEL"].get(level);
        (*itor)["INDEX"].get(index);
        (*itor)["COUNT"].get(count);
        HealPixel p(index, level,2*(level-m_minlevel));
        // set energy for center of this bin
        double energy( m_emin*pow(eratio, level-m_minlevel+0.5)), time(0.);
        // and its direction
        SkyDir sdir(p());
        
        // create its band, using the binner (assuming consistent!)
        addPhoton( astro::Photon(sdir, energy, time, 0), count);
        m_pixels ++;
    }
    delete &table; 
    std::cout << "Photons available: " << stored_photons 
        << "  Pixels available: " << stored_pixels <<std::endl;
    std::cout << "Photons loaded: " << m_photons 
        << "  Pixels created: " << m_pixels <<std::endl;
    setName("BinnedPhotonData from " +inputFile); // default name is the name of the file
    try{
        // now load the GTI info, if there
        gti() = Gti(inputFile);
        std::cout << "  GTI interval: "
            << int(gti().minValue())<<"-"<<int(gti().maxValue())<<std::endl; 
    }catch(const std::exception&){
        std::cerr << "BinnedPhotonData:: warning: no GTI information found" << std::endl;
    }
}

void BinnedPhotonData::addPhoton(const astro::Photon& gamma, int count)
{   
    // create a emmpty band with this photon's properties
    Band newband (m_binner(gamma));
    int key(newband);
    
    // is it already in our list?
    iterator it=std::lower_bound(begin(), end(), key, std::less<int>());

    if( key!=(*it) ){
        // no, create new entry and copy in the Band
        it = insert(it, newband);
    }

    // now add the counts to the band's pixel
    (*it).add(gamma.dir(), count);

    m_photons+= count;
}


double BinnedPhotonData::density (const astro::SkyDir & sd) const
{
    double result(0);
    static double norm((M_PI/180)*(M_PI/180) ); // normalization factor: 1/degree

    for (const_iterator it = begin(); it!=end(); ++it) {
        const Band& band ( *it);
        int count( band(sd) );
        result += count / band.pixelArea();
    }
    return result*norm;
}

double BinnedPhotonData::value(const astro::SkyDir& dir, double e)const
{
    double result(0);

    for( const_iterator it=begin();  it!=end(); ++it)  {
        const Band& band = *it;
        if( e< band.emin() || e >= band.emax() ) continue;
        result += band(dir);
    }
    return result;

}

double BinnedPhotonData::integral(const astro::SkyDir& dir, double a, double b)const
{

    return value(dir, sqrt(a*b));
}

void BinnedPhotonData::info(std::ostream& out)const
{
    int total_pixels(0), total_photons(0);
    out << "index  nside type   emin    emax   sigma    pixels   photons\n";

    int i(0);
    for( const_iterator it=begin();  it!=end(); ++it, ++i)
    {
        const Band& band = *it;
        int pixels(band.size()), photons(band.photons());
        out <<std::setw(4) << i
            <<std::setw(8) << band.nside()
            <<std::setw(5) << band.event_class()
            <<std::setw(7) << int(band.emin()+0.5)
            <<std::setw(8) << int(band.emax()+0.5)
            <<std::setw(8) << int(band.sigma()*180/M_PI*3600+0.5)
            <<std::setw(10)<< pixels
            <<std::setw(10)<< photons 
            <<std::endl;
        total_photons += photons; total_pixels+=pixels;
    }
    out << " total"
        <<std::setw(44)<<total_pixels
        <<std::setw(10)<<total_photons << std::endl;
}

void BinnedPhotonData::write(const std::string & outputFile,
            const std::string & tablename,
            bool clobber) const
{
    throw std::runtime_error("BinnedPhotonData::write -- Not implemented"); 
}

void BinnedPhotonData::writegti(const std::string & outputFile) const
{
    m_gti.writeExtension(outputFile);
}
void BinnedPhotonData::addgti(const skymaps::Gti& other)
{
    m_gti |= other;
}

void BinnedPhotonData::operator+=(const skymaps::BinnedPhotonData& other) {
#if 0 //TODO rewrite this
    for(const_iterator it=other.begin();it!=other.end();++it){
        addPixel(it->first,it->second);
    }
#endif
    addgti(other.gti());
}

