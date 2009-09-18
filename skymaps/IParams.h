/** @file IParams.h
    @brief declare class IParams

$Header$

*/
#ifndef skymaps_IParams_h
#define skymaps_IParams_h

#include <string>
#include <vector>

#include "skymaps/Exposure.h"

namespace skymaps {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class IParams
    @brief a class that manages the parameters of the PSF from a database
*/

    class IParams {
public:
    
    /** @brief sigma returns PSF parameter sigma
        @param energy 
        @param event_class 0=front, 1=back
        */
    static double sigma(double energy, int event_class);

    /** @brief gamma returns PSF parameter gamma
        @param energy 
        @param event_class 0=front, 1=back
        */
    static double gamma(double energy, int event_class);

    /** @brief effaWeight returns effective area weight for event class
        @param energy 
        @param event_class 0=front, 1=back
        */
    static double effaWeight(double energy, int event_class);

    //sets front and back parameterization for sigma
    // s = (a[0]^2+a[2]*(E/100)^(a[1])+a[3]*(E/100)^(2*a[1]))^0.5
    static void set_fp(std::vector<double> params);
    static void set_bp(std::vector<double> params);
    static void set_fsig(std::vector<double> sigs);
    static void set_bsig(std::vector<double> sigs);

    //sets gamma values for front and back set for 4/decade from 100MeV->100GeV
    //with a bin for <100 MeV and one for >100GeV (15 in all)
    static void set_fgam(std::vector<double> gams);
    static void set_bgam(std::vector<double> gams);

    //set front and back effective area weights
    static void set_fwght(std::vector<double> wght);
    static void set_bwght(std::vector<double> wght);

    //set energy bins
    static void set_elist(std::vector<double> elist);

    //run the initialization of psf parameters
    //through CALDB files
    static void init(const std::string& name="P6_v3",
        const std::string& clevel="diff",const std::string& file="");

    //set the CALDB directory
    static void set_CALDB(const std::string& dir);

    static void set_livetimefile(const std::string& ltfile);
    static void set_skydir(const astro::SkyDir& dir);

    //returns 
    static std::vector<double> params(int event_class);
    static double scale(double energy, int event_class);

    // output
    static void print_parameters();

private:

    /** @brief ctor not used
    */
    IParams();

    static std::vector<double> s_elist; //energy bin list 4/decade
    static std::vector<double> s_fparams; //front width parameterization
    static std::vector<double> s_bparams; //back width parameterization
    static std::vector<double> s_fsig;
    static std::vector<double> s_bsig;
    static std::vector<double> s_fgam; //front tail parameterization
    static std::vector<double> s_bgam; //back tail parameterization
    static std::vector<double> s_fwght; //front effective area weight
    static std::vector<double> s_bwght; //back effective area weight

    static std::string s_CALDB;
    static std::string s_livetimefile;

    static bool s_init; //checks to see whether initialization has been performed

    static astro::SkyDir dir;
    
};

class TopHat {
public:
    /**
    @param c_lo lower cosine limit
    @param c_hi upper cosine limit
    */
    TopHat( double c_lo, double c_hi ) : m_clo(c_lo) , m_chi(c_hi) {}

    double operator()(double costh) const
    {
        return ( (costh > m_clo) && (costh <= m_chi) ) ? 1 : 0;
    }
    double m_clo, m_chi;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class ExposureWeighter
@helper class to provide proper weighting of PSF over effective area and/or livetime

*/

class ExposureWeighter {
public:

    ExposureWeighter(const std::string& faeff_str, const std::string& baeff_str, const std::string& livetimefile) {
        
        m_faeff = new EffectiveArea("",faeff_str);
        m_baeff = new EffectiveArea("",baeff_str);

        if (livetimefile.empty()) {
            m_uselt = false;
        }
        else {
            m_uselt = true;
            m_lt = new LivetimeCube(livetimefile);
        }
    }

    ~ExposureWeighter() {
        delete m_faeff;
        delete m_baeff;
        if (m_uselt) { delete m_lt; }
    }


    double operator()(double c_lo, double c_hi, double e_lo, double e_hi, int event_class, astro::SkyDir& dir) {

        EffectiveArea* aeff = (event_class == 0 ? m_faeff : m_baeff);
        double estep ( log(e_hi/e_lo) / 8 ),
               cstep ( (c_hi - c_lo)  / 4 ),
               eratio( exp(estep) ),
               ec(2.),
               e(e_lo),
               ae(0);

        double loop_result,c,cc;

        if (m_uselt) {
            for( int i = 0; i < 9; ++i) {
                c  = c_lo;
                cc = 2.;
                loop_result = aeff->value(e,c)* m_lt->value(dir,c);
                for ( int j = 1; j < 4; ++j) { 
                    c += cstep;
                    cc = 6 - cc;
                    loop_result += cc * aeff->value(e,c) * m_lt->value(dir,c);
                }
                loop_result += aeff->value(e,c_hi)* m_lt->value(dir,c);
                if (i==0 || i==8) {
                    ae += e * loop_result;
                }
                else {
                    ae += e * ec * loop_result;
                    ec = 6 - ec;
                }
                e *= eratio;
            }
            ae *= (cstep*estep/9);
        }
        else {
            for( int i = 0; i < 9; ++i) {
                c  = c_lo;
                cc = 2.;
                loop_result = aeff->value(e,c);
                for ( int j = 1; j < 4; ++j) { 
                    c += cstep;
                    cc = 6 - cc;
                    loop_result += cc * aeff->value(e,c);
                }
                loop_result += aeff->value(e,c_hi);
                if (i==0 || i==8) {
                    ae += e * loop_result;
                }
                else {
                    ae += e * ec * loop_result;
                    ec = 6 - ec;
                }
                e *= eratio;
            }
            ae *= (cstep*estep/9);
        }
        return ae;
        //double ae(event_class == 0 ? m_faeff->value(e, (c_hi+c_lo)/2.) : m_baeff->value(e, (c_hi+c_lo)/2.));
        
        if (m_uselt) {
            TopHat fun(c_lo,c_hi);
            return ae * (m_lt->bins(dir))(fun);
        }
        return ae;
    }

private:

    EffectiveArea* m_faeff;
    EffectiveArea* m_baeff;
    LivetimeCube* m_lt;
    bool m_uselt;
};


}
#endif

