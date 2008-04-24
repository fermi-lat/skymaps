/** @file BinnedPhoton.h
@brief declare class BinnedPhoton

$Header$

*/

#include "astro/SkyDir.h"

#ifndef skymaps_BinnedPhoton_h
#define skymaps_BinnedPhoton_h
#include <vector>

namespace skymaps {

    class Band;

    class BinnedPhoton {
    public:
        BinnedPhoton(const skymaps::Band* band=0, unsigned int index=0)
            : m_band(band)
            , m_index(index)
        {}
         
        BinnedPhoton(const skymaps::Band* band, const astro::SkyDir dir);

        const Band* band() const {return m_band;}

        unsigned int index()const { return m_index;}
        astro::SkyDir dir()const;  

        bool invalid()const {return (m_band==0);}

        void query_disk( double radius, std::vector<int>& v)const;

    private:
        const Band* m_band;
        unsigned int m_index; ///< pixel index
    };
}
#endif
