/** @file PsfSkySpectrum.h
@brief declare PsfSkySpectrum.

$Header$
*/
#ifndef skymaps_PsfSkySpectrum_h
#define skymaps_PsfSkySpectrum_h

#include "astro/SkyDir.h"
#include "skymaps/SkySpectrum.h"

namespace skymaps{

    class PsfFunction;

    /**
    @class PsfSkyFunction
    @brief Represents the distribution of the PSF function 

    */

    class PsfSkySpectrum : public skymaps::SkySpectrum{
    public:
        /** @brief ctor
        @param roi_dir center of PSF
        */

        PsfSkySpectrum(const astro::SkyDir & roi_dir);

        /// @brief the PSF value at the specified direction and energy
        double value(const astro::SkyDir &r, double energy) const;

        double integral(const astro::SkyDir &dir, double a, double b) const;

        /// @brief update the PSF direction
        void set_skydir(const astro::SkyDir &nd);
        /// @brief set the event class
        void set_event_class(int ec);

        /// @brief turn on/off jacobean -- default is to give dN/dOmega
        void set_jacobean(bool b);


    private:
        astro::SkyDir m_dir;
        int m_ec;
        bool m_jac;
    };

}

#endif
