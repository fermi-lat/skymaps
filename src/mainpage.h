// $Header$
// Mainpage for doxygen

/*! \mainpage package skymaps

   \author  Toby Burnett, Marshall Roth, Bruce Lesnick 


 Class hierarchy:

- astro::SkyFunction - abstract base class defining a real function on the sphere (a function of a astro::SkyDir)
- skymaps::SkyImage - implement SkyFunction with FITS image; also create from a SkyFunction
- skymaps::SkySpectrum - abstract, allow specification of energy spectrum at any point
- skymaps::DiffuseFunction - interpolate a FITS cube. Used for the background for point source fits.
- skymaps::PhotonMap - pixelized photon data 
- skymaps::Convolution - convolution of a SkySpectrum object with another SkySpectrum, perhaps a PSF.
- skymaps::CompositeSkySpectrum - linear combination of SkySpectrum objects. Used to combine the galactic diffuse with nearby (< 1deg) strong sources
- skymaps::Exposure - Integrate an expouse cube over the acceptance to define the exposure at any point.


\section notes release notes
  release.notes
\section requirements requirements
\include requirements

*/

