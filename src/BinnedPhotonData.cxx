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
#include <errno.h>

using astro::SkyDir;
using astro::Photon;

using namespace skymaps;

namespace {
    skymaps::PhotonBinner default_binner;

    std::string band_table("BANDS"),
        pixel_table("PIXELS");

}
BinnedPhotonData::BinnedPhotonData(int bins_per_decade)
: m_binner(bins_per_decade)
{}

BinnedPhotonData::BinnedPhotonData(const skymaps::PhotonBinner& binner)
: m_binner(binner)
{}


BinnedPhotonData::BinnedPhotonData(const std::string & inputFile,  const std::string band_table)
: m_binner(default_binner) // should change, or make flexible
, m_photons(0)
{
    const tip::Table * ptable(0);
    try{
        ptable = tip::IFileSvc::instance().readTable(inputFile, band_table);
    }catch(const std::exception&){
        std::cerr << "BinnedPhotonData: table "<< band_table << " not found in file " << inputFile<< std::endl;
        throw;
    }
    const tip::Table& table(*ptable);  // reference for convenience

    if( band_table == "PHOTONMAP")  // Old style
    {
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
        std::cout << "Creating BinnedPhotonData from file " << inputFile << ", table " << band_table << std::endl;

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
        delete ptable; 
        std::cout << "Photons available: " << stored_photons 
            << "  Pixels available: " << stored_pixels <<std::endl;
        std::cout << "Photons loaded: " << m_photons 
            << "  Pixels created: " << m_pixels <<std::endl;
        setName("BinnedPhotonData from " +inputFile); // default name is the name of the file
    }

    else // New style
    {
        // Get band info
        clear();  // clear band list
        const tip::Header& hdr = table.getHeader();
        int version_number(0), stored_bands(0), stored_pixels(0), stored_photons(0), pixels_loaded(0);

        // Guard against headers not being found in fits file.  Set to default on error

        try
        {
            hdr["VERSION"].get(version_number);
        }
        catch (const std::exception& )
        {
            version_number = 1;
        }

        if (version_number < 1 || version_number > 1)
        {
            throw std::runtime_error(std::string("BinnedPhotonData:: Invalid input file version."));
        }


        hdr["NBRBANDS"].get(stored_bands); 
        hdr["PIXELS"].get(stored_pixels);
        hdr["PHOTONS"].get(stored_photons);

        std::cout << "Creating Bands, pixels from file " << inputFile << ", tables " << band_table; 
        std::cout << ", " << pixel_table <<std::endl;
        
        std::vector<int> counts;


        for(tip::Table::ConstIterator itor = table.begin(); itor != table.end(); ++itor)
        {
            int nside, event_class, count;
            double emin, emax, sigma, gamma;
            (*itor)["NSIDE"].get(nside);
            (*itor)["EVENT_CLASS"].get(event_class);
            (*itor)["EMIN"].get(emin);
            (*itor)["EMAX"].get(emax);
            (*itor)["SIGMA"].get(sigma);
            (*itor)["GAMMA"].get(gamma);
            (*itor)["COUNT"].get(count);
            Band b(nside, event_class, emin, emax, sigma, gamma);
            push_back(b);
            counts.push_back(count);

        }
        delete ptable; 

        // Get pixel info from the pixel table
        ptable = tip::IFileSvc::instance().readTable(inputFile, pixel_table);
        const tip::Table& table2(*ptable);
        m_photons = 0;

        BinnedPhotonData::iterator band_iterator = this->begin();
        tip::Table::ConstIterator itor = table2.begin();
        std::vector<int>::const_iterator citor = counts.begin();

        // now just copy
        for(; band_iterator != this->end(); ++band_iterator, ++citor) // for each band
        {
            for (int i = 0; i < *citor; ++i, ++itor) // for number of pixels stored for this band
            {
                int index, count;
                (*itor)["INDEX"].get(index);
                (*itor)["COUNT"].get(count);
                (*band_iterator).add(index,count);
                m_photons += count;
                ++pixels_loaded;
            }
        }
        delete &table2; 

        std::cout << "\tBands expect:   " << stored_bands 
                  << ", found:  " << size() <<std::endl;
        std::cout << "\tPixels expect:  " << stored_pixels 
                  << ", found:  " << pixels_loaded <<std::endl;
        std::cout << "\tPhotons expect: " << stored_photons 
                  << ", found:  " << m_photons <<std::endl;

    }

    try // now load the GTI info, if there
    {
        gti() = Gti(inputFile);
        std::cout << "  GTI interval: "
            << int(gti().minValue())<<"-"<<int(gti().maxValue())<<std::endl; 
    }
    catch(const std::exception&)
    {
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
    int newkey(*it);

    if( key!=(*it) ){
        // no, create new entry and copy in the Band
        it = insert(it, newband);
    }

    // now add the counts to the band's pixel
    (*it).add(gamma.dir(), count);

    m_photons+= count;
}
void BinnedPhotonData::add(const BinnedPhotonData& other)
{
    iterator mit(begin());
    for(const_iterator oit( other.begin()); oit!=other.end(); ++oit,++mit){
        mit->add(*oit);
    }
    addgti(other.gti());
    m_photons += other.m_photons;

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
    out << "index  emin    emax class   sigma   nside    pixels   photons\n";

    int i(0);
    for( const_iterator it=begin();  it!=end(); ++it, ++i)
    {
        const Band& band = *it;
        int pixels(band.size()), photons(band.photons());
        out <<std::setw(4) << i
            <<std::setw(7) << int(band.emin()+0.5)
            <<std::setw(8) << int(band.emax()+0.5)
            <<std::setw(6) << band.event_class()
            <<std::setw(8) << int(band.sigma()*180/M_PI*3600+0.5) // convert to arcsec
            <<std::setw(8) << band.nside()
            <<std::setw(10)<< pixels
            <<std::setw(10)<< photons 
            <<std::endl;
        total_photons += photons; total_pixels+=pixels;
    }
    out << " total"
        <<std::setw(45)<<total_pixels
        <<std::setw(10)<<total_photons << std::endl;
}

void BinnedPhotonData::write(const std::string & outputFile, bool clobber) const
{

    int version_number(1); /* Use this number to indicate when the layout of the fits file changes.  This will allow
                           the read() function to interpret and input all defined output forrmats. */

    if (clobber)
    {
        int rc = std::remove(outputFile.c_str());
        if( rc == -1 && errno == EACCES ) 
            throw std::runtime_error(std::string(" Cannot remove file " + outputFile));
    }

    unsigned int total_pixels(0), total_photons(0);

    {
    // First, add header table to the file
    tip::IFileSvc::instance().appendTable(outputFile, band_table);
    tip::Table & table = *tip::IFileSvc::instance().editTable( outputFile, band_table);

    table.appendField("NSIDE", "1J");
    table.appendField("EVENT_CLASS", "1J");
    table.appendField("EMIN", "1D");
    table.appendField("EMAX", "1D");
    table.appendField("SIGMA", "1D");
    table.appendField("GAMMA", "1D");
    table.appendField("COUNT", "1J"); // Number of pixels in this band
    table.setNumRecords(size());

    // get iterators for the Table and the Band list
    tip::Table::Iterator itor = table.begin();
    BinnedPhotonData::const_iterator band_iterator = this->begin();


    // now just copy
    for( ; band_iterator != this->end(); ++band_iterator, ++itor)
    {
        (*itor)["NSIDE"].set(band_iterator->nside());
        (*itor)["EVENT_CLASS"].set(band_iterator->event_class());
        (*itor)["EMIN"].set(band_iterator->emin());
        (*itor)["EMAX"].set(band_iterator->emax());
        (*itor)["SIGMA"].set(band_iterator->sigma());
        (*itor)["GAMMA"].set(band_iterator->gamma());
        (*itor)["COUNT"].set(band_iterator->size());
        total_pixels += band_iterator->size();
        total_photons += band_iterator->photons();
    }

    // set the headers (TODO: do the comments, too)
    tip::Header& hdr = table.getHeader();
    hdr["NAXIS1"].set(3 * sizeof(long) + 4 * sizeof(double));
    hdr["NBRBANDS"].set(size()); 
    hdr["PIXELS"].set(total_pixels);
    hdr["PHOTONS"].set(total_photons);
    hdr["VERSION"].set(version_number);

    // close it?
    delete &table;
    }
    {
    // Now, add pixels from the pixel table
    tip::IFileSvc::instance().appendTable(outputFile, pixel_table);
    tip::Table & table = *tip::IFileSvc::instance().editTable( outputFile, pixel_table);

    table.appendField("INDEX", "1J"); // Healpix index for pixel
    table.appendField("COUNT", "1J"); // Number of photons in this pixel
    table.setNumRecords(total_pixels);

    // initialize iterator for the Table 
    tip::Table::Iterator itor = table.begin();

    for(const_iterator band_iterator = begin(); band_iterator != end(); ++band_iterator)  // For each band
    {
        // Output pixel info
        for(std::map<int, int>::const_iterator pitor = band_iterator->begin(); pitor != band_iterator->end(); ++pitor, ++itor)
        {
            (*itor)["INDEX"].set(pitor->first);
            (*itor)["COUNT"].set(pitor->second);
        }
    }

    // set the headers (TODO: do the comments, too)
    tip::Header& hdr = table.getHeader();
    hdr["NAXIS1"].set(2 * sizeof(long));
    hdr["NBRBANDS"].set(size()); 

    // close it?
    delete &table;
    }
    // Now set the gti
    m_gti.writeExtension(outputFile);
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

