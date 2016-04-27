/** @file Photon.h
@brief Declaration of Photon class, derived from astro::Photon

$Header$
*/

#ifndef skymaps_Photon_h
#define skymaps_Photon_h

#include "astro/Photon.h"
#include "astro/SkyDir.h"


namespace skymaps{

    /** @class Photon
    @brief derive from astro::Photon to add event type
    */
    class Photon : public astro::Photon {
    public:
        /// ctor sets photon, also rotation info.
        Photon(const astro::SkyDir& dir, double energy, 
				double time, int event_class=-1, int source=-1,
				int event_type=-1)
			   : astro::Photon(dir,energy,time,event_class,source)
			   // If event type not specified, set it to event_class if given
			   , m_event_type(((event_class>=0)&&(event_type<0))?event_class:event_type)
		{}

        int event_type()const{ return m_event_type;}
    private:
        int m_event_type;
	};
    
}


#endif
