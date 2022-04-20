include(UseSWIG)

add_library(
  skymaps SHARED
  src/Background.cxx
  src/Band.cxx
  src/BandBackground.cxx
  src/BinnedPhotonData.cxx
  src/CompositeSkyFunction.cxx
  src/CompositeSkySpectrum.cxx
  src/Convolution.cxx
  src/DiffuseFunction.cxx
  src/EffectiveArea.cxx
  src/Exposure.cxx
  src/ExposureWeighter.cxx
  src/Gti.cxx
  src/GtiBase.cxx
  src/HealpixDiffuseFunc.cxx
  src/IParams.cxx
  src/IsotropicConstant.cxx
  src/IsotropicPowerLaw.cxx
  src/IsotropicSpectrum.cxx
  src/LivetimeCube.cxx
  src/PhotonBinner.cxx
  src/PhotonMap.cxx
  src/PsfFunction.cxx
  src/PsfSkyFunction.cxx
  src/PsfSkySpectrum.cxx
  src/PySkyFunction.cxx
  src/PySkySpectrum.cxx
  src/PythonPsf.cxx
  src/PythonUtilities.cxx
  src/SkyImage.cxx
  src/SkyImageSpectrum.cxx
  src/SkyIntegrator.cxx
  src/SkySpectrum.cxx
  src/SmoothedSkySpectrum.cxx
  src/SpectralFunction.cxx
  src/WeightedSkyDir.cxx
)
target_include_directories(
  skymaps PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

target_link_libraries(skymaps PUBLIC facilities tip st_facilities healpix embed_python timeSystem)

set(UseSWIG_TARGET_NAME_PREFERENCE STANDARD)

### SWIG generates the python wrapper library
set_property(SOURCE src/swig_setup.i PROPERTY CPLUSPLUS ON)

swig_add_library(
  swig_skymaps
  TYPE SHARED
  LANGUAGE python
  SOURCES src/swig_setup.i
)

set_property(
  TARGET swig_skymaps
  PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE
)
set_property(TARGET swig_skymaps PROPERTY SWIG_MODULE_NAME _skymaps)

target_include_directories(
  swig_skymaps PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

target_link_libraries(
  swig_skymaps
  PRIVATE skymaps
  astro
  CLHEP::CLHEP
  Python3::Python
)

###############################################################
# Installation
###############################################################
install(DIRECTORY skymaps DESTINATION ${FERMI_INSTALL_INCLUDEDIR})
install(DIRECTORY data/ DESTINATION ${FERMI_INSTALL_DATADIR})

install(
  TARGETS skymaps swig_skymaps
  # EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
