/**
* @file EffectiveArea.h
* @brief declare class EffectiveArea. 
*
* Author T. Burnett
* $Header$
*/

#ifndef skymaps_EffectiveArea_h
#define skymaps_EffectiveArea_h

#include <string>
#include <map>



namespace skymaps {

    /**
    * @class EffectiveArea
    * Interpolate a table of effective area vs. energy and angle with respect to the boresight
    *
    * Must set the path to the CALDB files, for example by
    *    EffectiveArea::set_CALDB("f:\\glast\\packages\\ScienceTools-v9r7p1\\irfs\\caldb\\v0r7p1\\CALDB\\data\\glast\\lat");
    *
    *  or by setting the environment variable CALDB.
    *
    */

    class EffectiveArea  {

    public:
        ///@brief ctor
        ///@param irfname name of irf, e.g. P6_v1_diff_front. Special value "simple" returns linear functionl Otherwise
        ///      the name is env('CALDB')/bcf/ea/aeff_<irfname>.fits 
        ///@param filename name of explicit file, if irfname not specified
        ///@param table_name Name of the HDU. Default EFFECTIVE_AREA for old IRF format, 
        ///
        EffectiveArea(std::string irfname="", std::string filename="", std::string table_name="EFFECTIVE AREA");
        ~EffectiveArea();

        ///@brief return a value 

        ///@param e energy in MeV
        ///@param costh cosine of theta
        double value(double energy=1000, double costh=1.0)const;

        ///@brief return a value 
        ///@param e energy in MeV
        ///@param costh cosine of theta
        double operator()(double energy=1000, double costh=1.0)const{return value(energy,costh);}
        
        ///@brief return livetime fraction dependent efficiency factors -  added 5/19/2010 EEW
        ///@param e energy in MeV
        std::pair<double,double> getLivetimeFactors(double energy)const;
        static void set_CALDB(std::string CALDB);

        static bool enable_cache(bool value=true);
    private:
        EffectiveArea(const EffectiveArea& other){} // hide copy ctor
        class FitsTable; // forward declaration
        bool m_simple;
        FitsTable * m_aeffTable;
        static std::string s_CALDB;
        static bool s_cache_enabled;
        class EfficiencyParameter;
        
        EfficiencyParameter * m_p0;
        EfficiencyParameter * m_p1; 
        bool m_haveEfficiencyPars;

        
        /// @class private class to hash a key for lookup
        class CacheKey{
        public:
            CacheKey(double logenergy, double costheta);
            operator unsigned int()const{return m_key;}
        private:
            unsigned int m_key;
        };

        typedef std::map<CacheKey, float> Cache;
        mutable Cache m_cache; ///< cache values for speed
    };
} // namespace skymaps

#endif // astro_EffectiveArea_h
