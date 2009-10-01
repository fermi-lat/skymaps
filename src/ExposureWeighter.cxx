/** @file ExposureWeighter.cxx
@brief implement class ExposureWeighter

$Header$

*/

#include "skymaps/ExposureWeighter.h"
#include "skymaps/EffectiveArea.h"
#include "skymaps/LivetimeCube.h"

using namespace skymaps;

ExposureWeighter::ExposureWeighter(const std::string& faeff_str, const std::string& baeff_str, const std::string& livetimefile) :
m_faeff(new EffectiveArea("",faeff_str)),
m_baeff(new EffectiveArea("",baeff_str))
{
    if (livetimefile.empty()) {
        m_uselt = false;
    }
    else {
        m_uselt = true;
        m_lt = new LivetimeCube(livetimefile);
    }

}

ExposureWeighter::~ExposureWeighter() {
    delete m_faeff;
    delete m_baeff;
    if (m_uselt) { delete m_lt; }
}

double ExposureWeighter::operator()(double c_lo, double c_hi, double e_lo, double e_hi, int event_class, astro::SkyDir& dir) {

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
