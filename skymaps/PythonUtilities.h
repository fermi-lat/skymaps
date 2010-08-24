/** @file PythonUtilities.h
    @author M. Kerr

    This file collects routines for sending/receiving data to/from the Python
    interpreter.
*/

#ifndef skymaps_PythonUtilities_h
#define skymaps_PythonUtilities_h

#include <vector>
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
};

}


#endif
