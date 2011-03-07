/** @file PythonUtilities.h
    @author M. Kerr

    This file collects routines for sending/receiving data to/from the Python
    interpreter.
*/

#ifndef skymaps_PythonUtilities_h
#define skymaps_PythonUtilities_h

#include <vector>
#include "skymaps/WeightedSkyDir.h"

#include "astro/SkyDir.h"
#include "astro/SkyFunction.h"

namespace skymaps {
  /// Class to group utilities in a single namespace.
class PythonUtilities
{
public:
    
    // This will take the flat grid at the equator, rotate it to the
    // ROI frame, and evaluate a SkyFunction over it.
    static void val_grid(double* rvals, int rvals_size,
                         const std::vector<double>& lons, const std::vector<double>&lats, 
                         const astro::SkyDir& roi_center, const astro::SkyFunction& myfunc);


    // This will rotate a grid of lons/lats from the ROI frame to the
    // flat grid frame on the equator.  (Presumably the lons/lats
    // are in the ROI.)  They are in the Galactic coordinate system.
    // @brief get the rotated lons and lats for a list of positions
    // @param rvals pre-allocated memory, 2 float for each member of sdirs
    // @param sdirs vector of grid points to rotate
    static void rot_grid (double *rvals, int rvals_size,
                          const std::vector<astro::SkyDir> sdirs,
                          const astro::SkyDir& roi_center);
    
    // Return the arclength (in radians) between the SkyDirs in sdirs
    // and the SkyDir specified by roi_center
    static void arclength(double *rvals, int rvals_size,
                          const std::vector<astro::SkyDir> sdirs,
                          const astro::SkyDir& roi_center);

    // Get a column of floats from a FITS file
    static void get_float_col(double *rvals, int rvals_size,
                              const std::string file_name,
                              const std::string table_name,
                              const std::string col_name);

    // convert barycentered times to topocentric (at LAT) MET
    // input should be TDB in MJD units
    // tolerance is in (s)
    static void tdb2met(double *rvals, int rvals_size,
                        double ra, double dec,
                        const std::string sc_file,
                        double tol);

    // convert MET (topocentric TT) to barycentered times
    // input should be MET in s; output is in TDB in s referenced to MET origin
    static void met2tdb(double *rvals, int rvals_size,
                        double ra, double dec,
                        const std::string sc_file);

    // convert MET (topocentric TT) to geocentric TT
    // input should be MET in s; output is in TT@geocenter referenced to MET origin
    static void met2geo(double *rvals, int rvals_size,
                        double ra, double dec,
                        const std::string sc_file);

    static void get_wsdl_weights(double *rvals, int rvals_size,
                                 const skymaps::BaseWeightedSkyDirList& wsdl);

    static void set_wsdl_weights(const std::vector<double>& weights,
                                 skymaps::BaseWeightedSkyDirList& wsdl);

};

}


#endif
