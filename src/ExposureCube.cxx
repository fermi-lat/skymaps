/** @file ExposureCube.cxx
    @brief Implementation of class ExposureCube

   $Header$
   copied from: /nfs/slac/g/glast/ground/cvs/map_tools/src/Exposure.cxx,v 1.33 2008/01/16 21:33:59 burnett Exp $
*/
#include "skymaps/ExposureCube.h"
#include "healpix/HealpixArrayIO.h"
#include "tip/Table.h"
#include "astro/EarthCoordinate.h"

#include <memory>
#include <algorithm>

using namespace skymaps;
using healpix::HealpixArrayIO;
using healpix::CosineBinner;
using healpix::Healpix;
using astro::SkyDir;
using CLHEP::Hep3Vector;


ExposureCube::ExposureCube(const std::string& inputfile, const std::string& tablename)
: SkyExposureCube(SkyBinner(2))
, m_zenith_frame(false)
{
   setData( HealpixArrayIO::instance().read(inputfile, tablename));
}

/// return the closest power of 2 for the side parameter
/// 41252 square degrees for the sphere
/// nside=64 is about one degee: 49152 pixels
inline int side_from_degrees(double pixelsize){ 
    int n = 1;
    while( 12*n*n < 41252/(pixelsize*pixelsize) && n < 512){
        n *=2;
    }
    return n; 
} 

ExposureCube::ExposureCube(double pixelsize, double cosbinsize, double zcut)
: SkyExposureCube(
    SkyBinner(Healpix(
      side_from_degrees(pixelsize),  // nside
      Healpix::NESTED, 
      astro::SkyDir::EQUATORIAL) )
  )
, m_zcut(zcut), m_lost(0)
{
    unsigned int cosbins = static_cast<unsigned int>(1./cosbinsize);
    if( cosbins != CosineBinner::nbins() ) {
        SkyBinner::iterator is = data().begin();
        for( ; is != data().end(); ++is){ // loop over all pixels
            CosineBinner & pixeldata= *is; // get the contents of this pixel
            pixeldata.resize(cosbins);
            CLHEP::Hep3Vector pixdir = data().dir(is)();
        }
        CosineBinner::setBinning(0, cosbins);
    }
    create_cache();
}

void ExposureCube::create_cache()
{
    size_t datasize(data().size());
    m_dir_cache.reserve(datasize);

    SkyBinner::iterator is = data().begin();
    for( ; is != data().end(); ++is){ // loop over all pixels
        Simple3Vector pixdir(data().dir(is)());
        m_dir_cache.push_back(std::make_pair(&*is, pixdir));
    }
}
/** @class Filler
    @brief private helper class used in for_each to fill a CosineBinner object
*/
class ExposureCube::Filler {
public:
    /** @brief ctor
        @param deltat time to add
        @param dir direction to use to determine angle (presumably the spacecraft z-axis)
        @param zenith optional zenith direction for potential cut
        @param zcut optional cut: if -1, ignore
    */
    Filler( double deltat, const astro::SkyDir& dir, astro::SkyDir zenith=astro::SkyDir(), double zcut=-1)
        : m_dir(dir())
        , m_zenith(zenith())
        , m_deltat(deltat)
        , m_zcut(zcut)
        , m_total(0), m_lost(0)
    {}
    void operator()( const std::pair<CosineBinner*, Simple3Vector> & x)
    {
        // check if we are making a horizon cut:
        bool ok( m_zcut==-1);
        if( ! ok) {
            double z(x.second.dot(m_zenith));
            ok = z > m_zcut;
        }
        if( ok) {
            // if ok, add to the angle histogram
            x.first->fill(x.second.dot(m_dir), m_deltat);
            m_total += m_deltat;
        }else{
            m_lost += m_deltat;
        }
    }
    double total()const{return m_total;}
    double lost()const{return m_lost;}
private:
    Simple3Vector m_dir, m_zenith;
    double m_deltat, m_zcut, m_total, m_lost;
};

void ExposureCube::fill(const astro::SkyDir& dirz, double deltat)
{
    Filler sum = for_each(m_dir_cache.begin(), m_dir_cache.end(), Filler(deltat, dirz));
    addtotal(deltat);
}


void ExposureCube::fill(const astro::SkyDir& dirz, const astro::SkyDir& zenith, double deltat)
{
    Filler sum = for_each(m_dir_cache.begin(), m_dir_cache.end(), Filler(deltat, dirz, zenith, m_zcut));
    double total(sum.total());
    addtotal(total);
    m_lost += sum.lost();
}


void ExposureCube::write(const std::string& outputfile, const std::string& tablename)const
{
    healpix::HealpixArrayIO::instance().write(data(), outputfile, tablename);
}

void ExposureCube::load(const tip::Table * scData, 
                    const GTIvector& gti, 
                    bool verbose) {
   
   tip::Table::ConstIterator it = scData->begin();
   const tip::ConstTableRecord & row = *it;
   long nrows = scData->getNumRecords();

   for (long irow = 0; it != scData->end(); ++it, ++irow) {
      if (verbose && (irow % (nrows/20)) == 0 ) std::cerr << ".";
      if( processEntry( row, gti) )break;
   }
   if (verbose) std::cerr << "!" << std::endl;
}


bool ExposureCube::processEntry(const tip::ConstTableRecord & row, const GTIvector& gti)
{
#if 0 // enable when use SAA?
    double latGeo, lonGeo;
    row["lat_Geo"].get(latGeo);
    row["lon_Geo"].get(lonGeo);
    astro::EarthCoordinate earthCoord(latGeo, lonGeo);
    if( earthCoord.insideSAA() ) return false;
#endif

    double  start, stop, livetime; 
    row["livetime"].get(livetime);
    row["start"].get(start);
    row["stop"].get(stop);
    double deltat = livetime > 0 ? livetime : stop-start;


    double fraction(1); 
    bool  done(false);
    if( !gti.empty() ) {
        fraction = 0;

        GTIvector::const_iterator it  = gti.begin();
        for ( ; it != gti.end(); ++it) {
            double first = it->first,
                second=it->second;

            if( start < first ) {
                if( stop < first) continue; // history interval before gti
                if( stop <= second){
                    fraction = (stop-first)/(stop-start); // overlap start of gti
                    break;
                }
                fraction = (second-first)/(stop-start); // gti subset of history
                break;
            }else {
                if( start > second) continue; // interval after gti 
                if( stop <= second ) {
                    fraction = 1.0; break;  // fully contained
                }
                fraction = (second-start)/(stop-start);  // overlap end of gti
                break;
            }
 
        }
        done = fraction==0 && start > gti.back().second; 
    }
    if( fraction>0. ) {

        double ra, dec, razenith, deczenith;
        row["ra_scz"].get(ra);
        row["dec_scz"].get(dec);
        row["ra_zenith"].get(razenith);
        row["dec_zenith"].get(deczenith);
        SkyDir zenith(razenith, deczenith);
        SkyDir scz(ra,dec);  

        if(m_zenith_frame){
            // special mode for Earth, or zenith frame
            double theta = scz.difference(zenith); 
            if( fabs(theta)<1e-8) theta=0;

            static Hep3Vector north_pole(0,0,1);
            Hep3Vector east_dir( north_pole.cross(zenith()).unit() ); // east is perp to north_pole and zenith
            Hep3Vector north_dir( zenith().cross(east_dir));

            double azimuth=atan2( scz().dot(east_dir), scz().dot(north_dir) ); // z is north, heading.
            if( azimuth <0) azimuth += 2*M_PI; // to 0-360 deg.
            // now convert to degrees
            azimuth *= 180/M_PI;
            theta *= 180/M_PI;

            // and finally to longitude, latitude in zenith frame
            azimuth +=90; if(azimuth>360)azimuth-=360;
            fill(SkyDir(azimuth, theta-90), deltat * fraction);

        }else{
        
            fill(scz, zenith, deltat* fraction);
        }
    }
    return done; 

}

