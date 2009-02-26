/** @file LivetimeCube.h
    @brief definition of the class LivetimeCube

    @author T.Burnett
    $Header$
    copied from:
    /nfs/slac/g/glast/ground/cvs/map_tools/map_tools/Exposure.h,v 1.23 2007/12/11 05:06:57 burnett Exp $
*/
#ifndef MAP_TOOLS_LivetimeCube_H
#define MAP_TOOLS_LivetimeCube_H


#include "astro/SkyDir.h"
#include "astro/SkyFunction.h"

#include "skymaps/Gti.h"

#include "healpix/HealpixArray.h"
#include "healpix/CosineBinner.h"

namespace tip { class Table; class ConstTableRecord;}

#include <utility> // for std::pair


/** @class BasicLivetime
@brief template for differential LivetimeCube

@param S Pixelization class, must implement dir(), is a list of C objects
@param C angular binner class, must implement operator()(const F&)const
@param F a function of one parameter
*/


template< class S, class C>
class BasicLivetime {
public:
    BasicLivetime(S sky):m_sky(sky), m_total(0){}

    virtual ~BasicLivetime(){}

    template<class F>
        double operator()(const astro::SkyDir& dir, const F& fun)const
    {
        const C& binner = m_sky[dir];
        return binner(fun);
    }
    const S& data()const{return m_sky;}
    S& data(){return m_sky;}
    double total()const{return m_total;}
    void addtotal(double t){ m_total+=t;}

    void setData(const S& data){m_sky=data;}

    const C& bins(const astro::SkyDir& dir)const{ return m_sky[dir]; }
private:
    S m_sky;
    double m_total;
};

// define LivetimeCube as specific instantiation of the above
typedef healpix::HealpixArray<healpix::CosineBinner> SkyBinner;
typedef BasicLivetime<SkyBinner, healpix::CosineBinner> SkyLivetimeCube;


namespace skymaps {
    class SkyImage; // forward declaration

/**
@class LivetimeCube
@brief Manage a differential LivetimeCube database.

It is a pixelated using Healpix binning, and the CosineBinner class

*/

    class LivetimeCube : public SkyLivetimeCube{
public:
    //! @brief ctor create object with specified binning
    //! @param inputfile [""] if set, load from file, ignore rest 
    //! @param dir direction to select (ignore if cone_angle is 180)
    //! @param cone_angle [180] half angle of cone (180 degrees default: all sky) 
    //! @param zcut [-1]  zenith angle cut angle
    //! @param pixelsize (deg) [1] Approximate pixel size, in degrees
    //! @param cosbinsize bin size in the cos(theta) binner
    //! Note that all parameters have defaults to allow keyword args in Python interface
    
    LivetimeCube(
          const std::string& inputfile=""
         ,const astro::SkyDir& dir=astro::SkyDir()
         ,double cone_angle=180
         ,double zcut=-1.0
         ,double pixelsize=1. 
         ,double cosbinsize=1./healpix::CosineBinner::nbins()
         ,double quiet=false
          );
  
    //! add a time interval at the given position
    virtual void fill(const astro::SkyDir& dirz, double deltat);

    //! write out to a file.
    void write(const std::string& outputfile, const std::string& tablename="EXPOSURE")const;

    //! load a set of history intervals from a FITS S/C file (FT2)
    void load(std::string scfile, const skymaps::Gti & gti= skymaps::Gti(), std::string tablename="SC_DATA");

    //! load a set of history intervals from a table, qualified by a set of "good-time" intervals 
    void load_table(const tip::Table * scData, 
        const skymaps::Gti & gti= skymaps::Gti(), 
                    bool verbose=true);

    double lost()const{return m_lost;}

    void useZenith(){m_zenith_frame=true;} ///< to convert to Earth coordinates

    //! @brief evaluate the livetime for the given direction,and the theta, phi
    //! @param dir
    //! @param costh
    //! @param phi (optional: if negative or not present, return total)
    double value(const astro::SkyDir& dir, double costh, double phi=-1);

    /// access time
    double total()const{return SkyLivetimeCube::total();} // make accessible to SWIG

    const skymaps::Gti& gti()const{return m_gti;}

private:
    bool processEntry(const tip::ConstTableRecord & row, const skymaps::Gti& gti);

   /** @brief  allow horizon cut, possible if FOV includes horizon
        @param dirz direction of z-axis of instrument
        @param dirzenith direction of local zenith
        @param deltat time interval
    */
    virtual void fill_zenith(const astro::SkyDir& dirz,const astro::SkyDir& dirx, const astro::SkyDir& dirzenith, double deltat);

    /** @brief set up the cache of vectors associated with cosine histograms

    */
    void create_cache();

    // for setting an output FITS LivetimeCube (copied from Likelihood::LikeExposure
    void computeCosbins(std::vector<double> & mubounds) const;
    void writeCosbins(const std::string & outfile) const;
    void fitsReportError(int status, const std::string & routine) const;
    void setCosbinsFieldFormat(const std::string & outfile) const;
    void writeFilename(const std::string & outfile) const;

    /** @class Simple3Vector 
    @brief replacement for Hep3Vector for speed of dot product

    */
    class Simple3Vector {
    public: 
        Simple3Vector(const CLHEP::Hep3Vector& v=CLHEP::Hep3Vector())
            : x(v.x()),y(v.y()),z(v.z()){};
        Simple3Vector(double a, double b, double c):x(a),y(b), z(c){}
        double dot(const Simple3Vector& u)const{return x*u.x+y*u.y+z*u.z;}
        CLHEP::Hep3Vector transform(const CLHEP::HepRotation& R)const{
            return CLHEP::Hep3Vector(
                R.xx()*x+R.xy()*y+R.xz()*z,
                R.yx()*x+R.yy()*y+R.yz()*z,
                R.zx()*x+R.zy()*y+R.zz()*z);
        }
        double x,y,z;
    };
    std::vector< std::pair<healpix::CosineBinner* ,  Simple3Vector> > m_dir_cache;
    class Filler ; ///< class used to fill a CosineBinner object with a value

    double m_zcut; ///< value for zenith angle cut
    double m_lost; ///< keep track of lost

    bool m_zenith_frame; ///< set true to make livetime map in Earth, or zenith frame.
    double m_cone_angle; ///< half-angle of cone used
    astro::SkyDir m_dir; ///< direction to sample. (ignore if 180 deg)

    skymaps::Gti m_gti; ///< over-all Gti
    bool m_quiet;       ///< suppress output if set
};


} // namespace map_tools
#endif
