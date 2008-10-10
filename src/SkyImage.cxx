/** @file SkyImage.cxx

@brief implement the class SkyImage
$Header$
*/

#include "skymaps/SkyImage.h"
#include "astro/SkyProj.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"

#include <cctype>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <stdexcept>
#include <errno.h> // to test result of std::remove()

namespace {
    static unsigned long lnan[2]={0xffffffff, 0x7fffffff};
    static double& dnan = *( double* )lnan;
}
using namespace skymaps;

double SkyImage::setNaN(double nan)
{ 
    double oldnan = dnan;
    dnan=nan;
    return oldnan;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SkyImage::SkyImage(const astro::SkyDir& center,  
                   const std::string& outputFile, 
                   double pixel_size, double fov, int layers, 
                   const std::string& ptype,
                   bool galactic, bool earth)
: m_naxis3(layers)  
, m_image(0)
, m_save(false)
, m_layer(0)
, m_interpolate(false)
{

    if( fov>90) {
        m_naxis1=0;
        std::string types[]={"" ,"CAR","AIT","ZEA", "SIN", "CEA"};
        int xsize[] =       {360, 360,  325,  230,  115,    360}; 
        int ysize[] =       {180, 180,  162,  230,  115,    115}; 
        for( unsigned int i = 0; i< sizeof(types)/sizeof(std::string); ++i){
            if( ptype == types[i]) {
                m_naxis1 = static_cast<int>(xsize[i]/pixel_size);
                m_naxis2 = static_cast<int>(ysize[i]/pixel_size);
                break;
            }
        }

        if( m_naxis1==0) {
            throw std::invalid_argument("SkyImage::SkyImage -- projection type " 
                +ptype +" does not have default image size");
        }
    }else{

        m_naxis1=m_naxis2 = static_cast<int>(fov/pixel_size + 0.5);
    }

    double crval[2] = { galactic?center.l():center.ra(),galactic? center.b(): center.dec()};
    double cdelt[2] = { earth? pixel_size: -pixel_size, pixel_size };
    double crpix[2] = { (m_naxis1+1)/2.0, (m_naxis2+1)/2.0};

    m_wcs = new astro::SkyProj(ptype, crpix, crval, cdelt, 0., galactic);
    m_pixelCount = m_naxis1*m_naxis2*m_naxis3;
    m_imageData.resize(m_pixelCount);

    if( ! outputFile.empty() ){
        this->setupImage(outputFile);
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SkyImage::setupImage(const std::string& outputFile,  bool clobber)
{
    std::string extension("skyimage"); // maybe a parameter?

    if( clobber ){
        int rc = std::remove(outputFile.c_str());
        if( rc==-1 && errno ==EACCES ) throw std::runtime_error(
            std::string("SkyImage: cannot remove file "+outputFile)
            );
    }

    // setup the image: it needs an axis dimension array and the file name to write to
    std::vector<long> naxes(3);
    naxes[0]=m_naxis1;
    naxes[1]=m_naxis2;
    naxes[2]=m_naxis3;

    // now add an image to the file
    tip::IFileSvc::instance().appendImage(outputFile, extension, naxes);
    // create a float image
    m_image = tip::IFileSvc::instance().editImageFlt(outputFile, extension);
    m_save=true;

    m_pixelCount = m_naxis1*m_naxis2*m_naxis3;
    m_imageData.resize(m_pixelCount);

    // fill the boundaries with NaN
    //if( pars.projType()!="CAR") clear();

    m_wcs->setKeywords(m_image->getHeader());
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SkyImage::SkyImage(const std::string& fits_file, const std::string& extension, bool interpolate)
: m_save(false)
, m_layer(0)
, m_wcs(0)
, m_interpolate(interpolate)
{
    // note expect the image to be float
    m_image = tip::IFileSvc::instance().editImageFlt(fits_file, extension);
    tip::Header& header = m_image->getHeader();

    int naxis;
    header["NAXIS"].get(naxis);
    header["NAXIS1"].get(m_naxis1);
    header["NAXIS2"].get(m_naxis2);
    if( naxis>2 ){
       header["NAXIS3"].get(m_naxis3);
    }else{ m_naxis3=1;}
    m_pixelCount = m_naxis1*m_naxis2*m_naxis3;

#if 0 // this seems not to interpret CAR properly
    m_wcs = new astro::SkyProj(fits_file,1);
#else
    m_wcs = new astro::SkyProj(fits_file, "");
#endif
    // finally, read in the image: assume it is float
    dynamic_cast<tip::TypedImage<float>*>(m_image)->get(m_imageData);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
unsigned int SkyImage::setLayer(unsigned int newlayer)
{
    checkLayer(newlayer);
    unsigned int t = m_layer;
    m_layer = newlayer;
    return t;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool SkyImage::addPoint(const astro::SkyDir& dir, double delta, unsigned int layer)
{
    std::pair<double,double> p= dir.project(*m_wcs);
    // ignore if not in the image.
    if( p.first<0 || p.first >= m_naxis1 || p.second<0 || p.second>=m_naxis2) return false;
    unsigned int 
        i = static_cast<unsigned int>(p.first),
        j = static_cast<unsigned int>(p.second),
        k = i+m_naxis1*(j + layer*m_naxis2);
    
    if(  k< m_pixelCount){
        m_imageData[k] += delta;
        m_total += delta;
    }
    return true;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SkyImage::checkLayer(unsigned int layer)const
{
    if( layer >= (unsigned int)m_naxis3){
        std::stringstream errmsg;
        errmsg << "SkyImage: requested layer " << layer 
            << " not compatible with axis3: " << m_naxis3 << std::endl;
        throw std::out_of_range(errmsg.str());
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SkyImage::fill(const astro::SkyFunction& req, unsigned int layer)
{
    checkLayer(layer);
    m_total=m_count=m_sumsq=0;
    m_min=1e10;m_max=-1e10;
    int offset = m_naxis1* m_naxis2 * layer;
    for( size_t k = 0; k< (unsigned int)(m_naxis1)*(m_naxis2); ++k){
        double 
            x = static_cast<int>(k%m_naxis1)+1.0, 
            y = static_cast<int>(k/m_naxis1)+1.0;
        if( m_wcs->testpix2sph(x,y)==0) {
            astro::SkyDir dir(x,y, *m_wcs);
            double t= req(dir);
            m_imageData[k+offset] = t;
            m_total += t;
            ++m_count;
            m_sumsq += t*t;
            m_min = t<m_min? t:m_min;
            m_max = t>m_max? t:m_max;
        }else{
            // not valid (off the edge, perhaps)
            m_imageData[k+offset]=dnan; 
        }
    }
    return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SkyImage::fill(const skymaps::SkySpectrum& req, unsigned int layer)
{
    fill(dynamic_cast<const astro::SkyFunction&>(req),layer);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SkyImage::clear()
{
    size_t s = m_imageData.size();
    for( size_t k = 0; k< s; ++k){
        m_imageData[k]=0; 
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SkyImage::~SkyImage()
{
    if( m_save) {
        dynamic_cast<tip::TypedImage<float>*>(m_image)->set(m_imageData);
    }
    delete m_image; 
    delete m_wcs;
}
void SkyImage::save()
{
    dynamic_cast<tip::TypedImage<float>*>(m_image)->set(m_imageData);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double SkyImage::pixelValue(const astro::SkyDir& pos,unsigned  int layer)const
{
    checkLayer(layer); 
    double v(0) ;
    
    if( m_interpolate ){
        // interpolating
        // project using wcslib interface
        std::pair<double,double> p= pos.project(*m_wcs);
        double x( floor(p.first-0.5) ),
               y( floor(p.second-0.5) ),
               dx( p.first-x -1 ),
               dy( p.second-y-1 );
        unsigned int k = static_cast<unsigned int>(x + m_naxis1*(y+layer*m_naxis2));
        v = m_imageData.at(k);
        if( dx>0. && x < m_naxis1){ 
            v = v*(1.-dx) + m_imageData[k+1]*dx;
        }else if( x >1) {
            v = v*(1+dx) - m_imageData[k-1]*dx;
        }
        if( dy>0. && y < m_naxis2){ 
            v = v*(1.-dy) + m_imageData[k+m_naxis1]*dy;
        }else if( y > 1) {
            v = v*(1.+dy) - m_imageData[k-m_naxis1]*dy;
        }
    }else{
        unsigned int k = pixel_index(pos,layer);
        v = m_imageData.at(k);
    }
   
    return v;        
    
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
float &  SkyImage::operator[](const astro::SkyDir&  pixel)
{
    unsigned int k = pixel_index(pixel);
    return m_imageData[k];        

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
const float &  SkyImage::operator[](const astro::SkyDir&  pixel)const
{
    unsigned int k = pixel_index(pixel);
    return m_imageData[k];        

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double SkyImage::operator()(const astro::SkyDir& s)const
{
    return pixelValue(s, m_layer);        
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SkyImage::getNeighbors(const astro::SkyDir& pos, std::vector<double>&neighbors)const
{
    int layer = 0; ///@todo: get neighbors on a different layer
    std::pair<double,double> p= pos.project(*m_wcs);
    if( p.first<0) p.first += m_naxis1;
    if(p.second<0) p.second += m_naxis2;
    unsigned int 
        i = static_cast<unsigned int>(p.first-0.5),
        j = static_cast<unsigned int>(p.second-0.5),
        k = i+m_naxis1*(j + layer*m_naxis2);
    if(i+1<(unsigned int)m_naxis1)neighbors.push_back(m_imageData[k+1]); 
    if(i>0) neighbors.push_back(m_imageData[k-1]);
    if(j+1<(unsigned int)m_naxis2)neighbors.push_back(m_imageData[k+m_naxis1]);
    if(j>0)neighbors.push_back(m_imageData[k-m_naxis1]);

}

// internal routine to convert a SkyDir to a pixel index
unsigned int SkyImage::pixel_index(const astro::SkyDir& pos, int layer) const
{
    // if not specified, use the data member
    if( layer<0 ) layer = m_layer;

    // project using wcslib interface, then adjust to be positive
    std::pair<double,double> p= pos.project(*m_wcs);
    if( p.first<0) p.first += m_naxis1;
    if(p.second<0) p.second += m_naxis2;
    unsigned int 
        i = static_cast<unsigned int>(p.first-0.5),
        j = static_cast<unsigned int>(p.second-0.5),
        k = i+m_naxis1*(j + layer*m_naxis2);
     if( k > m_pixelCount ) {
        throw std::range_error("SkyImage::pixel_index -- outside image hyper cube");
    }
    return k;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SkyImage::reimage( const astro::SkyDir& center,
                        const std::string& outputFile, 
                   double pixel_size, double fov 
                   ,const std::string& ptype
                   ,bool galactic)
{
     SkyImage image(center,outputFile, pixel_size, fov, m_naxis3, ptype, galactic);
     for(int layer(0); layer!= m_naxis3; ++layer){
         image.fill(*this, layer);
     }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void SkyImage::setEnergies(const std::vector<double>& energies)
{
    if (energies.size()!=m_naxis3) {
        throw std::invalid_argument("SkyImage::setEnergies: wrong size for energy array, must be same as layers");
    }
    m_energy.resize(energies.size());
    std::copy(energies.begin(), energies.end(), m_energy.begin());
}

