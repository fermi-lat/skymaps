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
        ///
        EffectiveArea(std::string irfname="", std::string filename="");
        ~EffectiveArea();

        ///@brief return a value 

        ///@param e energy in MeV
        ///@param costh cosine of theta
        double value(double energy=1000, double costh=1.0)const;

        ///@brief return a value 
        ///@param e energy in MeV
        ///@param costh cosine of theta
        double operator()(double energy=1000, double costh=1.0)const{return value(energy,costh);}

        static void set_CALDB(std::string CALDB);
    private:
        class FitsTable; // forward declaration
        bool m_simple;
        FitsTable * m_aeffTable;
        static std::string s_CALDB;
    };
} // namespace skymaps

#endif // astro_EffectiveArea_h
