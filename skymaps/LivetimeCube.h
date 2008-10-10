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
@param C angualar binner class, must implement operator()(const F&)
@param F a function of one parameter
*/


template< class S, class C>
class BasicLivetime {
public:
    BasicLivetime(S sky):m_sky(sky), m_total(0){}

    virtual ~BasicLivetime(){}

    virtual void fill(const astro::SkyDir& dirz, double deltat)=0;
    virtual void fill(const astro::SkyDir& dirz, const astro::SkyDir& dirzenith, double deltat)=0;

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

    const C& bins(const astro::SkyDir& dir){ return m_sky[dir]; }
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

    class LivetimeCube : public SkyLivetimeCube, public astro::SkyFunction {
public:
    //! create object with specified binning
    //! @param pixelsize (deg) Approximate pixel size, in degrees
    //! @param cosbinsize bin size in the cos(theta) binner
    LivetimeCube(double pixelsize=1., double cosbinsize=1./healpix::CosineBinner::nbins(), double zcut=-1.0);

    //! create object from the data file (FITS for now)
    LivetimeCube(const std::string& inputfile, const std::string& tablename="Exposure");

    //! add a time interval at the given position
    virtual void fill(const astro::SkyDir& dirz, double deltat);

    /** @brief this was added by Julie, to allow horizon cut, possible if FOV includes horizon
        @param dirz direction of z-axis of instrument
        @param dirzenith direction of local zenith
        @param deltat time interval
    */
    virtual void fill(const astro::SkyDir& dirz, const astro::SkyDir& dirzenith, double deltat);

    //! write out to a file.
    void write(const std::string& outputfile, const std::string& tablename="Exposure")const;

    //typedef std::vector<std::pair<double, double> > GTIvector;

    //! load a set of history intervals from a FITS S/C file (FT2)
    void load(std::string scfile, const skymaps::Gti & gti= skymaps::Gti(), std::string tablename="SC_DATA");

    //! load a set of history intervals from a table, qualified by a set of "good-time" intervals 
    void load(const tip::Table * scData, 
        const skymaps::Gti & gti= skymaps::Gti(), 
                    bool verbose=true);

    //! implement the SkyFunction interface
    double operator()(const astro::SkyDir& sdir)const; 

    double lost()const{return m_lost;}

    void useZenith(){m_zenith_frame=true;} ///< to convert to Earth coordinates

    double value(const astro::SkyDir& dir, double costh);

    //! create a map in a given sky image
    void createMap(skymaps::SkyImage & image);

    /// access time
    double total()const{return SkyLivetimeCube::total();} // make accessible to SWIG


private:
    bool processEntry(const tip::ConstTableRecord & row, const skymaps::Gti& gti);

    /** @brief set up the cache of vectors associated with cosine histograms

    */
    void create_cache();

    /** @class Simple3Vector 
    @brief replacement for Hep3Vector for speed of dot product

    */
    class Simple3Vector {
    public: 
        Simple3Vector(const CLHEP::Hep3Vector& v=CLHEP::Hep3Vector())
            : x(v.x()),y(v.y()),z(v.z()){};
        double dot(const Simple3Vector& u)const{return x*u.x+y*u.y+z*u.z;}
        double x,y,z;
    };
    std::vector< std::pair<healpix::CosineBinner* ,  Simple3Vector> > m_dir_cache;
    class Filler ; ///< class used to fill a CosineBinner object with a value

    double m_zcut; ///< value for zenith angle cut
    double m_lost; ///< keep track of lost

    bool m_zenith_frame; ///< set true to make livetime map in Earth, or zenith frame.
};


} // namespace map_tools
#endif
