/** @file PythonPsf.h
    @author M. Kerr

    This class takes setup parameters from Python but implements the evaluation
    of the PSF in C++ for speed.
*/

#ifndef skymaps_PythonPsf_h
#define skymaps_PythonPsf_h

#include <vector>
#include "astro/SkyDir.h"
#include "skymaps/Band.h"

namespace skymaps {

    typedef  std::vector<double> svd;
    typedef  std::vector<double>::const_iterator svd_cit;
    //typedef  void (*base_func_ptr) (double, double*, const svd*, const svd*);

class PythonPsf
{
public:
    PythonPsf(const svd& sigma1, const svd& sigma2,
              const svd& gamma1, const svd& gamma2,
              const svd& nc, const svd& nt,
              const svd& weights);

    PythonPsf(const svd& sigma1_in,
              const svd& gamma1_in,
              const svd& weights_in);

    ~PythonPsf();

    double operator() (const astro::SkyDir& sd1, const astro::SkyDir& sd2);
    double operator() (double delta);
    double integral(double dmin, double dmax);
    double integral_from_zero(double delta);
    double inverse_integral_on_axis(double frac);
    double overlap_circle(const astro::SkyDir& center, double radius, const astro::SkyDir& src_pos);
    double overlap_healpix(int nside_base, long index, int nside_sub, const astro::SkyDir& src_pos);
    
    void array_val(double* rvals, int rvals_size);
    void wsdl_val(double* rvals, int rvals_size,
                  astro::SkyDir& src_pos, skymaps::WeightedSkyDirList& data_pos);

    static void set_density(bool dens);
    static bool get_density();

    void print_parameters();
    //double c_simps_a(double (*func)(double),double a, double b,int nsimps_0,double atol_in,double rtol_in);

private:
    static bool density;
    
    const svd* sigma1;
    const svd* sigma2;
    const svd* gamma1;
    const svd* gamma2;
    const svd* nc;
    const svd* nt;
    const svd* weights;
    
    bool newstyle;
    
    const int mysize;

    double* result1;
    double* result2;

    void psf_base(double delta, double* result, svd_cit sit, svd_cit git, svd_cit wit);
    void int_base(double delta, double* result, svd_cit sit, svd_cit git, svd_cit wit);
};
}

#endif