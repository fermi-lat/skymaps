/** @file PsfSkyFunction.h
@brief declare PsfSkyFunction.

$Header$
*/
#ifndef skymaps_PsfSkyFunction_h
#define skymaps_PsfSkyFunction_h

#include "astro/SkyDir.h"
#include "astro/SkyFunction.h"
#include "skymaps/PsfFunction.h"

namespace skymaps{

    class PsfFunction;

    /**
    @class PsfSkyFunction
    @brief Represents the distribution of the PSF function 

    */

    class PsfSkyFunction : public astro::SkyFunction{
    public:
        /** @brief ctor
        @param roi_dir center of PSF
        @param gamma value of the power law index (around 2)
        @param sigma value for the sigma parameter (radians)
        */

        PsfSkyFunction(const astro::SkyDir & roi_dir, double gamma, double sigma);

        /// @brief implement interface by returning the value
        double operator () (const astro::SkyDir & r)const;

        /// @brief the average value integrated for the region
        double average(const astro::SkyDir& dir, double angle, double tolerance)const;

        double level_ave(const astro::SkyDir& dir, double angle, int level) const;
   


    private:
        astro::SkyDir m_dir;
        skymaps::PsfFunction m_psf;
        double m_sigma;
    };

}

#endif
