/** @file SkyImage.cxx

@brief implement the class SkyImage
$Header$
*/

#include "skymaps/SkyImage.h"
#include "astro/SkyProj.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"
#include "st_facilities/Util.h"


#include <cctype>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <stdexcept>
#include <errno.h> // to test result of std::remove()
#include "fitsio.h"


namespace {
    static unsigned long lnan[2]={0xffffffff, 0x7fffffff};
    static double& dnan = *( double* )lnan;

       void fitsReportError(int status) {
      fits_report_error(stderr, status);
      if (status != 0) {
          throw std::string("skymaps::SkyImage " +
                           std::string("cfitsio error."));
      }
   }

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
, m_outfile(outputFile)
, m_ax1_offset(1)
, m_ax2_offset(1)
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
    double crpix[2] = { (m_naxis1+1)/2.0, (m_naxis2+1)/2.0}; // center pixel; WCS convention is that center of a pixel is a half-integer

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

    m_outfile = outputFile; // Kerr added; this prevent writeEnergies from failing if setupImage called independently

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

    header["CRPIX1"].get(m_ax1_offset);
    header["CRPIX2"].get(m_ax2_offset);

    m_ax1_offset = m_ax1_offset - int(m_ax1_offset) == 0.5 ? 1.0 : 1.0;
    m_ax2_offset = m_ax2_offset - int(m_ax2_offset) == 0.5 ? 1.0 : 1.0;

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
            //x = static_cast<int>(k%m_naxis1)+0.5, // center of pixel is half/integer! 
            //y = static_cast<int>(k/m_naxis1)+0.5;
            x = static_cast<int>(k%m_naxis1)+0 + m_ax1_offset, // wcs is 1-indexed
            y = static_cast<int>(k/m_naxis1)+0 + m_ax2_offset;

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
        if( m_energy.size()>0) writeEnergies();
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

    /*
    if( m_interpolate ){
        if( layer<0 ) layer = m_layer;
        // project using wcslib interface
        std::pair<double,double> p= pos.project(*m_wcs);
        if (p.first < 0) {p.first += m_naxis1; }
        if (p.second < 0) {p.second += m_naxis2; }

        double x( floor(p.first-0.5 ) ),
               y( floor(p.second-0.5) ),
               dx( p.first-x-1 ),
               dy( p.second-y-1 );

        unsigned int k = static_cast<unsigned int>(x + m_naxis1*(y+layer*m_naxis2));
        v = m_imageData.at(k);
        if( dx>0. && x < m_naxis1){ 
            v = v*(1.-dx) + m_imageData[k+1]*dx;
        }else if( x >1 ) {
            v = v*(1+dx) - m_imageData[k-1]*dx;
        }
        if( dy>0. && y < m_naxis2){ 
            v = v*(1.-dy) + m_imageData[k+m_naxis1]*dy;
        }else if( y > 1 ) {
            v = v*(1.+dy) - m_imageData[k-m_naxis1]*dy;
        }
    }
    */
    
    if (m_interpolate) {
        // using a bilinear interpolation scheme
        // project using wcslib interface
        if( layer<0 ) layer = m_layer;
        std::pair<double,double> p= pos.project(*m_wcs);

        p.first  -= m_ax1_offset; // convert from wcs base-1 indexing
        p.second -= m_ax2_offset;

        int x1(p.first),y1(p.second),x2(x1+1),y2(y1+1);
        double dx(p.first - x1),dy(p.second - y1);

        unsigned int offset = layer * m_naxis1 * m_naxis2;

        //if (x1 < 0 || y1 < 0 || x2 > m_naxis1 - 1 || y2 > m_naxis2 -1) {
        //    std::cout << x1 << "\t" << x2 << "\t" << y1 << "\t" << y2 << std::endl;
        //}

        // protect against going over edges -- not sure why this should happen, but it does.
        if (x1 < 0) { x1 = 0; x2 = 0; }
        if (y1 < 0) { y1 = 0; y2 = 0; }
        if (x2 > m_naxis1 -1) { x2 = m_naxis1 - 1; x1 = m_naxis1 - 1; }
        if (y2 > m_naxis2 -1) { y2 = m_naxis2 - 1; y1 = m_naxis2 - 1; }
              
        double v11( m_imageData[offset + x1 + y1*m_naxis1] );
        double v12( m_imageData[offset + x1 + y2*m_naxis1] );
        double v21( m_imageData[offset + x2 + y1*m_naxis1] );
        double v22( m_imageData[offset + x2 + y2*m_naxis1] );

        return v11*(1-dx)*(1-dy) + v12*(1-dx)*dy + v21*dx*(1-dy) + v22*dx*dy;
    }
    else{
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
    p.first  -= m_ax1_offset;
    p.second -= m_ax2_offset;
    
    // round pixel value and subtract 1 (wcs indexing starts at 1)
    int i = static_cast<int>(p.first);
    i = p.first  - i >= 0.5 ? i : i - 1;
    int j = static_cast<int>(p.second);
    j = p.second - j >= 0.5 ? j : j - 1;
    i += 1;
    j += 1;

    //if ((!m_isPeriodic && (i < 0 || i >= m_nlon)) 
    //    || j < 0 || j >= m_nlat) {
    //   return 0;
    //}

    //if ( (i < 0 || i > m_naxis1) || j < 0 || j >= m_naxis2 ) {
    //    return 0;
    //}
    if (i < 0 && i >= -1) {
        i = 0;
    }
    return (layer*m_naxis2 + j)*m_naxis1 + i;
    /*
    // if not specified, use the data member
    if( layer<0 ) layer = m_layer;

    // project using wcslib interface, then adjust to be positive
    std::pair<double,double> p= pos.project(*m_wcs);
    if( p.first<0) p.first  += m_naxis1;
    if(p.second<0) p.second += m_naxis2;
    unsigned int 
        i = static_cast<unsigned int>(p.first-0.5),
        j = static_cast<unsigned int>(p.second-0.5),
        k = i+m_naxis1*(j + layer*m_naxis2);
    if( k > m_pixelCount ) {
        throw std::range_error("SkyImage::pixel_index -- outside image hyper cube");
    }
    return k;
    */
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

void SkyImage::writeEnergies(){
    if( m_outfile.empty() ){
        throw std::invalid_argument("SkyImage::writeEnergies: no file");
    }
    std::string tablename("ENERGIES");
    const std::string& filename(m_outfile);

    // Check if the extension exists already. If not, add it.
    try{
        std::auto_ptr<const tip::Table> 
            gtiTable(tip::IFileSvc::instance().readTable(filename, tablename));
    } catch (tip::TipException & eObj) {
        if (!st_facilities::Util::expectedException(eObj, "Could not open FITS extension"))
        { throw;
        }
        int status(0);
        fitsfile * fptr;
        fits_open_file(&fptr, filename.c_str(), READWRITE, &status);
        ::fitsReportError(status);
        char *ttype[] = {"ENERGY"};
        char *tform[] = {"D"};
        char *tunit[] = {"MeV"};
        int cols(1);

        fits_create_tbl(fptr, BINARY_TBL, 0, cols, ttype, tform, tunit,
            (char*)tablename.c_str(), &status);
        ::fitsReportError(status);

        fits_close_file(fptr, &status);
        ::fitsReportError(status);
    }
   std::auto_ptr<tip::Table> 
      gtiTable(tip::IFileSvc::instance().editTable(filename, tablename));


   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   tip::Table * table = fileSvc.editTable(filename,tablename);

   table->setNumRecords(m_energy.size());

   tip::Table::Iterator it(table->begin());
   tip::TableRecord & row(*it);
   
   for (size_t i(0); i < m_energy.size(); i++, ++it) {
      row["ENERGY"].set(m_energy.at(i));
   }
   delete table;
}

