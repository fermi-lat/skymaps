/** @file LivetimeCube.cxx
    @brief Implementation of class LivetimeCube

   $Header$
   copied from: /nfs/slac/g/glast/ground/cvs/map_tools/src/Exposure.cxx,v 1.33 2008/01/16 21:33:59 burnett Exp $
*/
#include "skymaps/LivetimeCube.h"
#include "skymaps/SkyImage.h"



#include "healpix/HealpixArrayIO.h"
#include "facilities/commonUtilities.h"
#include "facilities/Util.h"

#include "fitsio.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "astro/EarthCoordinate.h"

#include <memory>
#include <algorithm>
#include <stdexcept>
#include <cmath>

using namespace skymaps;
using healpix::HealpixArrayIO;
using healpix::CosineBinner;
using healpix::Healpix;
using astro::SkyDir;
using CLHEP::Hep3Vector;


namespace {


}  // anonymous namespace
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

LivetimeCube::LivetimeCube
         (const std::string& inputfile
         ,const astro::SkyDir& dir
         ,double cone_angle
         ,double zcut
         ,double pixelsize
         ,double cosbinsize
         ,double quiet
         )
: SkyLivetimeCube(
    SkyBinner(Healpix(
      side_from_degrees(pixelsize),  // nside
      Healpix::NESTED, 
      astro::SkyDir::EQUATORIAL) )
  )
, m_zcut(zcut), m_lost(0)
, m_zenith_frame(false)
, m_cone_angle(cone_angle)
, m_dir(dir)
, m_quiet(quiet)
{
    if( !inputfile.empty() ) {
        static std::string tablename("EXPOSURE");
        setData( HealpixArrayIO::instance().read(inputfile, tablename));
        return;
    }
    unsigned int cosbins = static_cast<unsigned int>(1./cosbinsize);
    if( cosbins != CosineBinner::nbins() ) {
        // only if changed from default
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

void LivetimeCube::create_cache()
{
    double cut( cos(m_cone_angle*M_PI/180) );
    if( m_cone_angle==0 ){ cut=0.9995; }
    size_t datasize(data().size());
    m_dir_cache.reserve(datasize);

    Simple3Vector reference(m_dir()); //
    SkyBinner::iterator is = data().begin();
    for( ; is != data().end(); ++is){ // loop over all pixels
        Simple3Vector pixdir(data().dir(is)());
        double ct(reference.dot(pixdir));
        if( ct> cut){
            m_dir_cache.push_back(std::make_pair(&*is, pixdir));
        }
    }
    if(  m_dir_cache.size() < data().size() && !m_quiet){
        std::cout << "LivetimeCube: Filling " << m_dir_cache.size() << "/" <<data().size() <<" pixels" <<std::endl;
    }
}
/** @class Filler
    @brief private helper class used in for_each to fill a CosineBinner object
*/
class LivetimeCube::Filler {
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

void LivetimeCube::fill(const astro::SkyDir& dirz, double deltat)
{
    Filler sum = for_each(m_dir_cache.begin(), m_dir_cache.end(), Filler(deltat, dirz));
    addtotal(deltat);
}


void LivetimeCube::fill_zenith(const astro::SkyDir& dirz, const astro::SkyDir& zenith, double deltat)
{
    Filler sum = for_each(m_dir_cache.begin(), m_dir_cache.end(), Filler(deltat, dirz, zenith, m_zcut));
    double total(sum.total());
    addtotal(total);
    m_lost += sum.lost();
}



void LivetimeCube::load_table(const tip::Table * scData, 
                    const Gti& gti, 
                    bool verbose) {
   
   tip::Table::ConstIterator it = scData->begin();
   const tip::ConstTableRecord & row = *it;
   long nrows = scData->getNumRecords();

   for (long irow = 0; it != scData->end(); ++it, ++irow) {
      if( verbose && (nrows>=20) && (irow % (nrows/20)) == 0 ) std::cerr << ".";
      processEntry( row, gti) ;
   }
   if (verbose) std::cerr << "!\t"<< total() << std::endl;
   m_gti |= gti;
}


bool LivetimeCube::processEntry(const tip::ConstTableRecord & row, const skymaps::Gti& gti)
{
    double  start, stop, livetime; 
    row["livetime"].get(livetime);  // assume this takes care of any entries during SAA
    if(livetime==0 ) return false;
    row["start"].get(start);
    row["stop"].get(stop);
    double deltat = livetime; 
    double fraction(1); 
    bool  done(false);
    if( gti.getNumIntervals()>0 ) {
        fraction = 0;
        Gti::ConstIterator it  = gti.begin();
        fraction = gti.getFraction(start,stop, it);
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
        
            fill_zenith(scz, zenith, deltat* fraction);
        }
    }
    return done; 
}
double LivetimeCube::value(const astro::SkyDir& dir, double costh)
{
    if( costh==0) costh=1e-9; // avoid error

    return bins(dir)[costh];
}

void LivetimeCube::load(std::string scfile, const skymaps::Gti & gti, std::string tablename)
{
     tip::Table * scData = tip::IFileSvc::instance().editTable(scfile, tablename);
     this->load_table(scData, gti, !m_quiet); // should allow a GTI?
     delete scData; // closes the file, I hope
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                  create FITS file
void LivetimeCube::write(const std::string& outputfile, const std::string& tablename)const
{
    std::string dataPath = 
        facilities::commonUtilities::getDataPath("Skymaps");
    std::string templateFile = 
        facilities::commonUtilities::joinPath(dataPath, "LivetimeCubeTemplate");
    tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
    fileSvc.createFile(outputfile, templateFile);
    writeFilename(outputfile);

    healpix::HealpixArrayIO::instance().write(data(), outputfile, tablename);
#if 0 // doesn't work?
    writeCosbins(outputfile);
#endif

    m_gti.writeExtension(outputfile);

}

//===================================================================
// copied from LivetimeCube.cxx
void LivetimeCube::setCosbinsFieldFormat(const std::string & outfile) const {
   int status(0);

   fitsfile * fptr(0);
   std::string extfilename(outfile + "[EXPOSURE]");
   fits_open_file(&fptr, extfilename.c_str(), READWRITE, &status);
   fitsReportError(status, "LivetimeCube::setCosbinsFieldFormat");
   
   int colnum(1); // by assumption
   fits_modify_vector_len(fptr, colnum, data().at(0).size(), &status);
   fitsReportError(status, "LivetimeCube::setCosbinsFieldFormat");

   fits_close_file(fptr, &status);
   fitsReportError(status, "LivetimeCube::setCosbinsFieldFormat");
}

void LivetimeCube::
computeCosbins(std::vector<double> & mubounds) const {
   bool sqrtbins(healpix::CosineBinner::thetaBinning() == "SQRT(1-COSTHETA)");
   double cosmin(healpix::CosineBinner::cosmin());
   size_t nbins(healpix::CosineBinner::nbins());
   mubounds.clear();
   for (int i(nbins); i >= 0; i--) {
      double factor(static_cast<double>(i)/nbins);
      if (sqrtbins) {
         factor *= factor;
      }
      mubounds.push_back(1. - factor*(1. - cosmin));
   }
}
void LivetimeCube::writeCosbins(const std::string & outfile) const {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   tip::Table * table = fileSvc.editTable(outfile, "CTHETABOUNDS");
   table->setNumRecords(healpix::CosineBinner::nbins());

   tip::Table::Iterator it(table->begin());
   tip::TableRecord & row(*it);
   
   std::vector<double> mubounds;
   computeCosbins(mubounds);

   for (size_t i(0); i < mubounds.size() -1; i++, ++it) {
      row["CTHETA_MIN"].set(mubounds.at(i));
      row["CTHETA_MAX"].set(mubounds.at(i+1));
   }
   delete table;
}

void LivetimeCube::
fitsReportError(int status, const std::string & routine) const {
   if (status == 0) {
      return;
   }
   fits_report_error(stderr, status);
   std::ostringstream message;
   message << routine << ": CFITSIO error " << status;
   throw std::runtime_error(message.str());
}
void LivetimeCube::writeFilename(const std::string & outfile) const {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   tip::Image * phdu(fileSvc.editImage(outfile, ""));
   phdu->getHeader()["FILENAME"].set(facilities::Util::basename(outfile));
   delete phdu;
}

