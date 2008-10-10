/** @file SkyImage.h

    @brief declare  the class SkyImage

    @author Toby Burnett <tburnett@u.washington.edu>
    $Header$

*/

#ifndef skymaps_SKYIMAGE_H
#define skymaps_SKYIMAGE_H

#include "skymaps/SkySpectrum.h"
#include "astro/SkyDir.h"

#include <string>
#include <vector>

// forward declarations of classes involved in implementaion
namespace tip   { class ImageBase; }
namespace astro { class SkyProj; }

namespace skymaps {
/**
    @class SkyImage
    @brief manage a FITS image containing float values

    Implements the class SkyFunction, to define a function over the sky derived
    from the image.

*/
class SkyImage : public astro::SkyFunction
{
public:

    /** @brief load an image from a file.
        @param filename name of the file, only FITS for now
        @param extension Name of an extension: if blank, assume primary
        @param interpolate [false] return interpolated values

    */
   
    SkyImage(const std::string& filename, const std::string& extension="", bool interpolate=false);

    /** @brief create an image, using the projection
        @param center coords of image center
        @param outputFile FITS file to write the image to. If null, do not setup
        @param pixel_size [0.5] degree size of indivitual pixel
        @param fov [20] (degrees) size of field of view, square if <90, full sky if>90
        @param layers [1] number of layers to allocate
        @param ptype ["ZEA"] projection type.
        @param galactic [false] use galactic or equatorial coords
        @param earth [false] looking at Earth: reverse x axis

        Note that if the outputFile is empty, which is useful to get the the transformation, most methods are not valid
        until setupImage() is called with a filename.
    */
    SkyImage(const astro::SkyDir& center,  
                   const std::string& outputFile, 
                   double pixel_size=0.5, double fov=20, int layers=1
                   ,const std::string& ptype="ZEA"
                   ,bool galactic=false
                   ,bool earth=false);
    
    /**
        @brief add a count to the map, using current SkyDir projection
        @param dir A SkyDir object
        @param delta incremental value (default 1 if not present)
        @param layer for multi-layer app. 0 (default) means the first layer 
        @return true if in the image, and added to it; false otherwise (not added)
    */
    bool addPoint(const astro::SkyDir& dir, double delta=1.0, unsigned int layer=0);

 
     /** @brief direct access to the pixel at the given direction and current layer
    */
    float & operator[](const astro::SkyDir&  pixel);

     const float & operator[](const astro::SkyDir&  dir)const;

    void save();

    ~SkyImage();

    //! set default layer, return previous 
    unsigned int setLayer(unsigned int newlayer);
 
    //! set, or return a reference to, the energies
    void getEnergies(std::vector<double> & energy) const { energy = m_energy; }
    const std::vector<double> & energies()const{return m_energy;}
    void setEnergies(const std::vector<double>& energies);

    /**
    @brief loop over all internal bins, request the intensity from a functor derived
    from SkyFunction
    @param req a functor that returns a double for a SkyDir
    @param layer layer number to fill [default 0]
    */
    void fill( const astro::SkyFunction& req, unsigned int layer=0);

    /// @brief needed since SWIG does not check inheritance?
    void fill( const skymaps::SkySpectrum& req, unsigned int layer=0);

    /** brief clear the image, putting nulls around a AIT map
    */
    void clear();

    //! @brief return the sum of all pixel values in the image
    double total()const{return m_total;}
    double minimum()const{return m_min;}
    double maximum()const{return m_max;}
    double count()const{return m_count;}

    /** @brief get value of the pixel at given skydir location
        @param pos position in the sky
        @param layer number
        @return value of the pixel corresponding to the given direction
    */
    double pixelValue(const astro::SkyDir& pos, unsigned int layer=0)const;
    
    /** @brief  set a list of the neighbor values
    @param pos position in the sky
    @param nlist list of neighbor values to set
    */
    void getNeighbors(const astro::SkyDir& pos, std::vector<double>& neighbors)const ;
    
    /// @brief implement SkyFunction interface by returning value at the selected pixel
    /// @param s the direction
    /// note that if there are multiple layers, it will choose the selected layer, see setLayer.
    double operator()(const astro::SkyDir& s)const;

    /// @brief access to number of layers
    int layers()const{return m_naxis3;}

    /// @brief access to the projection object
    const astro::SkyProj* projector()const { return m_wcs; }

    /// @brief setup a FITS image 
    void setupImage(const std::string& outputFile,  bool clobber=true);

    int naxis1()const{return m_naxis1;}
    int naxis2()const{return m_naxis2;}

    const std::vector<float>& image()const{return m_imageData;}

    /** @brief create new FITS file with different projection, or portion
        @param center coords of image center
        @param outputFile FITS file to write the image to. 
        @param pixel_size [0.5] degree size of indivitual pixel
        @param fov [20] (degrees) size of field of view
        @param ptype ["ZEA"] projection type.
        @param galactic [false] use galactic or equatorial coords
        
        If multiple layers, will do all.
    */
    void reimage( const astro::SkyDir& center,
                   const std::string& outputFile, 
                   double pixel_size=0.5, double fov=20 
                   ,const std::string& ptype="ZEA"
                   ,bool galactic=false);



    /// @brief change the value to use for invalid
    static double setNaN(double nan);

private:
    //! @brief internal routine to convert SkyDir to pixel index
    unsigned int pixel_index(const astro::SkyDir& pos, int layer=-1) const;

    /// @brief internal routine to check layer, or perhaps extend
    void checkLayer(unsigned int layer)const;

    //! sizes of the respective axes.
    int   m_naxis1, m_naxis2, m_naxis3;

    //! for statistics of a fill
    double m_total, m_sumsq, m_count, m_min, m_max;

    //! pointer to the associated tip Image class. The type is the base class, with 
    //! dynamic cast when necessary.
    tip::ImageBase* m_image;
    //! the actual image data
    std::vector<float>m_imageData;
    //! energy bounds for layers
    std::vector<double>m_energy;

    unsigned int m_pixelCount;
    bool m_save; 
    unsigned int m_layer;

    /// associated projection object, initialized from a par file, or a FITS header
    astro::SkyProj* m_wcs; 

    bool m_interpolate; ///< flag to determine if interpolate
};
} //namesace skymaps

#endif
