#include <cstdlib>
#include "skymaps/PythonPsf.h"
#include "skymaps/WeightedSkyDir.h"
#include "healpix/Healpix.h"

namespace {
    
typedef skymaps::svd svd;

const svd* vec_product(const svd& arg1, const svd& arg2) {
    svd* myvec = new svd(arg1);
    svd::const_iterator a2_it = arg2.begin();
    for (svd::iterator mv_it(myvec->begin()); mv_it != myvec->end(); ++mv_it) {
        *mv_it *= *a2_it++;
    }
    return (const svd*)(myvec);
}

class CSimpsQuad
{
public:
    CSimpsQuad() {}
    virtual double integrand(double x)=0;
    double quad(double a, double b) {
        m_ns0 += (m_ns0 & 1);
        svd fv; fv.reserve(m_ns0);
        double w( (b-a)/m_ns0 ),r(integrand(a)),c(2),x(a);
        fv.push_back(r);
        for (int i(1); i < m_ns0; ++i) {
            c = 6 - c;
            x += w;
            fv.push_back(integrand(x));
            r += c*fv.back();
        }
        fv.push_back(integrand(b));
        r += fv.back();
        r *= (w/3);
        double r2(0);
        for (int i(0); i < 10; ++i) {
            r2 = fv.back()-fv.front();
            int nsimps(m_ns0*pow(2.,i+1));
            double w( (b-a)/nsimps ), w2(2*w),x(a+w);
            fv.reserve(nsimps);
            for (int j(0); j < nsimps/2; ++j) {
                fv.push_back(integrand(x));
                x += w2;
                r2 += 2*fv[j] + 4*fv.back();
            }
            r2 *= (w/3);
            if ((abs(r-r2) < m_atol) || (abs(r/r2-1) < m_rtol)) return r2;
            r = r2;
        }
        return r2;
    }
private:
    static double m_atol,m_rtol;
    static int m_ns0;
};

double CSimpsQuad::m_atol = 1e-4;
double CSimpsQuad::m_rtol = 1e-4;
int CSimpsQuad::m_ns0  = 8;



class Interior : public CSimpsQuad
{
public:
    Interior(double offset, double roi_radius, skymaps::PythonPsf* psf)
    : m_off(offset)
    , m_rr_sq(roi_radius*roi_radius)
    , m_psf(psf) {}
    
    virtual double integrand(double x) {
        double c(cos(x));
        double eff_rad( sqrt(m_rr_sq + m_off*m_off*(c*c - 1)) - m_off*c);
        return m_psf->integral_from_zero(eff_rad);
    }
    double calculate() { return quad(0,M_PI)/M_PI; }

private:
    double m_off,m_rr_sq;
    skymaps::PythonPsf* m_psf;
};

class Exterior : public CSimpsQuad
{
public:
    Exterior(double offset, double roi_radius, skymaps::PythonPsf* psf)
    : m_off(offset)
    , m_rr_sq(roi_radius*roi_radius)
    , m_ulimit(asin(roi_radius/offset))
    , m_psf(psf) {}
    
    virtual double integrand(double x) {
        if (x==m_ulimit) return 0;
        double c(cos(x));
        double r2( sqrt(m_rr_sq + m_off*m_off*(c*c - 1)) );
        double de( m_off*c - r2 );
        return m_psf->integral(de+2*r2,de);
    }
    double calculate() { return quad(0,m_ulimit)/M_PI; }

private:
    double m_off,m_rr_sq,m_ulimit;
    skymaps::PythonPsf* m_psf;
};

class XSQ { public: double operator()(double x){return x*x;} };
class FX  { public: double operator()(double x){return cos(100*x)*exp(x);} };

}

namespace skymaps {

bool PythonPsf::density(true);
void PythonPsf::set_density(bool dens) {density=dens;}
bool PythonPsf::get_density() {return density;}

PythonPsf::PythonPsf(const svd& sigma1_in, const svd& sigma2_in,
                     const svd& gamma1_in, const svd& gamma2_in,
                     const svd& nc_in, const svd& nt_in,
                     const svd& weights_in)
: sigma1(vec_product(sigma1_in,sigma1_in)) // NB pre-multiplication
, sigma2(vec_product(sigma2_in,sigma2_in))
, gamma1(new const svd(gamma1_in))
, gamma2(new const svd(gamma2_in))
, nc(vec_product(nc_in,weights_in)) // NB pre-multiplication
, nt(vec_product(nt_in,weights_in))
, weights(0)
, newstyle(true)
, m_size(sigma1_in.size())
, m_result1(new double[m_size])
, m_result2(new double[m_size])
{}

PythonPsf::PythonPsf(const svd& sigma1_in,
                     const svd& gamma1_in,
                     const svd& weights_in) 
: sigma1(vec_product(sigma1_in,sigma1_in))
, sigma2(0)
, gamma1(new const svd(gamma1_in))
, gamma2(0)
, nc(0)
, nt(0)
, weights(new const svd(weights_in))
, newstyle(false)
, m_size(sigma1_in.size())
, m_result1(new double[m_size])
, m_result2(0)
{}

PythonPsf::~PythonPsf() {
    delete sigma1;
    delete gamma1;
    delete[] m_result1;
    if (newstyle) {
        delete sigma2;
        delete gamma2;
        delete nc;
        delete nt;
        delete[] m_result2;
    }
    else{
        delete weights;
    }
}

void PythonPsf::psf_base(double delta, double* result, svd_cit sit, svd_cit git, svd_cit wit) {
    double d(0.5*delta*delta);
    for (double* ptr (result); ptr< result + m_size; ++ptr,++sit,++wit,++git) {
        double u(d/(*sit)),g(*git);
        double f( (1.-1./g)*pow(1.+u/g,-g) );
        *ptr = (*wit)*f ;
    }
}

void PythonPsf::int_base(double delta, double* result, svd_cit sit, svd_cit git, svd_cit wit) {
    double d(0.5*delta*delta);
    for (double* ptr (result); ptr< (result + m_size); ++ptr,++sit,++wit,++git) {
        double u( d/(*sit) ),g(*git);
        double f( pow(1.+u/g,1.-g) );
        *ptr = (*sit)*(*wit)*f;
    }
}

double PythonPsf::operator() (const astro::SkyDir& sd1, const astro::SkyDir& sd2) {
    return operator()(sd1.difference(sd2));
}

double PythonPsf::operator() (double delta) {
    double result(0);
    if (newstyle) {
        psf_base(delta,m_result1,sigma1->begin(),gamma1->begin(),nc->begin());
        psf_base(delta,m_result2,sigma2->begin(),gamma2->begin(),nt->begin());
        for (double* ptr (m_result2); ptr < (m_result2 + m_size); ++ptr) result += *ptr;
    }
    else {
        psf_base(delta,m_result1,sigma1->begin(),gamma1->begin(),weights->begin());
    }
    for (double* ptr (m_result1); ptr < (m_result1 + m_size); ++ptr) result += *ptr;
    return result;
}

double PythonPsf::integral_from_zero(double delta) {
    double result(0);
    if (newstyle) {
        int_base(delta,m_result1,sigma1->begin(),gamma1->begin(),nc->begin());
        int_base(delta,m_result2,sigma2->begin(),gamma2->begin(),nt->begin());
        for (double* ptr (m_result2); ptr < (m_result2 + m_size); ++ptr) result += *ptr;
    }
    else {
        int_base(delta,m_result1,sigma1->begin(),gamma1->begin(),weights->begin());
    }
    for (double* ptr (m_result1); ptr < (m_result1 + m_size); ++ptr) result += *ptr;
    return 1 - (2*M_PI)*result; // 2pi from normalization
}

double PythonPsf::integral(double dmax, double dmin) {
    if (dmin==0.) {return integral_from_zero(dmax);}
    return integral_from_zero(dmax) - integral_from_zero(dmin);
}

double PythonPsf::inverse_integral_on_axis(double frac) {
    double g1((*gamma1)[0]),s1((*sigma1)[0]);
    double u1( g1*(pow(1-frac,1./(1-g1)) -1));
    if (newstyle) {
        double g2((*gamma2)[0]),s2((*sigma2)[0]);
        double w1((*nc)[0]),w2((*nt)[0]);
        double u2( g1*(pow(1-frac,1./(1-g1)) -1));
        return (sqrt(2*u1)*s1*w1 + sqrt(2*u2)*s2*w2)/(w1+w2);
    }
    return sqrt(2*u1)*s1;
}

void PythonPsf::array_val(double* rvals, int rvals_size) {
    for (double* ptr(rvals); ptr < rvals + rvals_size; ++ptr) {
        *ptr = operator()(*ptr);
    }
}

void PythonPsf::wsdl_val(double* rvals, int rvals_size,
              astro::SkyDir& src_pos, skymaps::BaseWeightedSkyDirList& data_pos) {
    skymaps::BaseWeightedSkyDirList::iterator wsdl_it(data_pos.begin());
    for (double* ptr(rvals); ptr < rvals + rvals_size; ++ptr) {
        *ptr = operator()(src_pos.difference(*wsdl_it++));
    }
}
double PythonPsf::overlap_circle(const astro::SkyDir& center, double radius, const astro::SkyDir& src_pos) {
    double delta(center.difference(src_pos));
    if (delta <= radius) {
        if (delta < 1e-5) return integral_from_zero(radius);
        Interior integrator(delta,radius,this);
        return integrator.calculate();
    }
    else {
        Exterior integrator(delta,radius,this);
        return integrator.calculate();
    }
}

double PythonPsf::overlap_healpix(int nside_base, long index, int nside_sub, const astro::SkyDir& src_pos) {
    //if ( (nside_base&(nside_base-1))==0 ) ...
    healpix::Healpix hp_base(nside_base,healpix::Healpix::RING,astro::SkyDir::GALACTIC);
    healpix::Healpix hp_sub(nside_sub,healpix::Healpix::RING,astro::SkyDir::GALACTIC);
    double accum(0),radius(sqrt(hp_base.pixelArea())); // conservative
    std::vector<int> iv;
    hp_sub.query_disc(healpix::Healpix::Pixel(index,hp_base),radius,iv);
    int tpix(0);
    return 0;
    for (std::vector<int>::iterator iv_it = iv.begin(); iv_it != iv.end(); ++iv_it) {
        healpix::Healpix::Pixel p_sub(*iv_it,hp_sub);
        healpix::Healpix::Pixel p_base(p_sub,hp_base);
        if (p_base.index() != index) continue;
        ++tpix;
        accum += operator()(src_pos.difference(p_sub));
    }    
    std::cout << tpix << std::endl;
    return accum*hp_sub.pixelArea();
}


void PythonPsf::print_parameters() {
    if (!newstyle) {
        svd::const_iterator sit = sigma1->begin();
        svd::const_iterator git = gamma1->begin();
        svd::const_iterator wit = weights->begin();
        std::cout << "Sigma   " << "Gamma   " << "Weight" << std::endl;
        for (int i(0); i < m_size; ++i) {
            std::cout << *sit++ << "  " << *git++ << "  " << *wit++ << std::endl;
        }
    }
}

//double PythonPsf::c_simps_a(double (*func)(double),double a, double b,int nsimps_0,double atol_in,double rtol_in) {
//    c_simps_atol = atol_in; c_simps_rtol = rtol_in;
//    return c_simps(func,a,b,nsimps_0);
//}

}

/*
double c_simps_atol(1e-4);
double c_simps_rtol(1e-4);
double c_simps(double (*func)(double),double a, double b,int nsimps_0) {
    nsimps_0 += (nsimps_0 & 1);
    svd fv; fv.reserve(nsimps_0);
    double w( (b-a)/nsimps_0 ),r(func(a)),c(2),x(a);
    fv.push_back(r);
    for (int i(1); i < nsimps_0; ++i) {
        c = 6 - c;
        x += w;
        fv.push_back(func(x));
        r += c*fv.back();
    }
    fv.push_back(func(b));
    r += fv.back();
    r *= (w/3);
    double r2(0);
    for (int i(0); i < 10; ++i) {
        r2 = fv.back()-fv.front();
        int nsimps(nsimps_0*pow(2,i+1));
        double w( (b-a)/nsimps ), w2(2*w),x(a+w);
        fv.reserve(nsimps);
        for (int j(0); j < nsimps/2; ++j) {
            fv.push_back(func(x));
            x += w2;
            r2 += 2*fv[j] + 4*fv.back();
        }
        r2 *= (w/3);
        if ((abs(r-r2) < c_simps_atol) || (abs(r/r2-1) < c_simps_rtol)) return r2;
        r = r2;
    }
    return r2;        
}
*/
