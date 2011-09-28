/** @file WeightedSkyDir.h
@brief declare class WeightedSkyDir, BaseWeightedSkyDirList,
and WeightedSkyDirList. WeightedSkyDir is a subclass of astro::SkyDir
also holds a weight. BaseWeightedSkyDir implmeents a list of SkyDirs
and 

Probably a better naming convention would be WeightedSkyDirList, and
BandWeightedSkyDirList, but is defined this way to preserve backwards
compatability.

$Header: 

*/

#include "skymaps/Band.h"
#include "astro/SkyFunction.h"
#include "astro/SkyDir.h"

#ifndef skymaps_WeightedSkyDir_h
#define skymaps_WeightedSkyDir_h

#include <vector>

namespace skymaps {

   /** @class WeightedSkyDir
         @brief a weighted SkyDir, used to describe pixels
    */
    class WeightedSkyDir : public astro::SkyDir {
    public:
        WeightedSkyDir(const astro::SkyDir& sdir=astro::SkyDir(), double weight=1): SkyDir(sdir), m_weight(weight){}
        double weight()const{return m_weight;}
        double set_weight(double weight) {m_weight=weight; return m_weight;}
    private:
        double m_weight;
    };

    /** @class WeightedSkyDirList
        @brief a vector of WeightedSKyDir objects

        */
    class BaseWeightedSkyDirList : public std::vector<WeightedSkyDir> {
    public:

        /// @brief calculate the arclength between a SkyDir and the elements of the list
        void arclength(const astro::SkyDir& sdir,std::vector<double>& output)const;
    };

    class WeightedSkyDirList : public astro::SkyFunction, public BaseWeightedSkyDirList {
    public:
        /** @brief ctor creates an object from a band
        @param sdir direction for cone
        @param radius cone half-angle (radians)
        @param include_empty [false] include empty pixels
        */
        WeightedSkyDirList( const Band& band, const astro::SkyDir& sdir, double radius, bool include_empty=false);

        /// @brief implement SkyFunction: return weight, or zero
        double operator()(const astro::SkyDir& sdir)const;

        /// @brief total pixels in list, including empties
        int total_pix()const {return m_pix;}
        int counts()const {return m_counts;}
    private:
        const skymaps::Band& m_band; ///< the band used to make
        int m_pix;
        int m_counts;
    };
}

#endif

