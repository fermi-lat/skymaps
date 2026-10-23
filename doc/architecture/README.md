---
package: skymaps
summary: NOT BUILT. Legacy C++/SWIG library of sky functions, spectra, binned photon data, PSF/exposure helpers (pointlike-era); present as a submodule but excluded from the Fermitools build.
layer: instrument
depends_on: [facilities, tip, st_facilities, healpix, embed_python, timeSystem, astro, CLHEP, Python3]
used_by: []
last_reviewed_commit: ac2444b
---

# skymaps Architecture

> **Status: not built.** `src/skymaps` is a git submodule of ScienceTools,
> but it is **not** in `FERMI_COMPONENTS` in the top-level `CMakeLists.txt`,
> so nothing here is compiled, installed or tested in a Fermitools build.
> It was dropped from that list in ScienceTools commit `a4e9834d`
> (2025-06-10, "Vendor WCSlib & XercesC"); the commit gives no reason. An
> earlier comment (ScienceTools `17ea9c4c`, 2022) said "Small ROOT
> dependence must be optionally removed from the swig module". The ROOT
> dependence was later excluded (skymaps `f22c158`).

## Purpose and scope

This is a library of functions defined on the sky, written by T. Burnett,
M. Roth and others for the UW `pointlike` analysis. That connection is
inferred from the authors, the commit messages and comments in
`BinnedPhotonData.h` and `SkyIntegrator.h`. The library provides:

- `SkySpectrum` (derived from `astro::SkyFunction`), which adds energy
  dependence. Its implementations are `DiffuseFunction` (interpolates a
  FITS cube), `IsotropicSpectrum`, `IsotropicPowerLaw`,
  `IsotropicConstant`, `CompositeSkySpectrum`, `Convolution` (HEALPix
  a_lm convolution with a PSF), `SmoothedSkySpectrum`, `Background`,
  `BandBackground`, `PsfSkySpectrum`, `SkyImageSpectrum`, and the Python
  callback adapters `PySkyFunction` and `PySkySpectrum`.
- `SkyImage`, a FITS image that implements `SkyFunction`.
- Binned photon data: `Photon`, `Band`, `BinnedPhotonData`,
  `PhotonBinner`, `PhotonMap` and `WeightedSkyDir`, all on HEALPix.
- Exposure and IRF helpers: `LivetimeCube`, `Exposure`, `EffectiveArea`,
  `IParams`, `PsfFunction` and `PythonPsf`. These read CALDB `aeff_*` and
  `psf_*` files **directly** from `$CALDB/bcf/...`, not through `irfs`.
- Local copies of `Gti`/`GtiBase`. `src/mainpage.h` says these were
  "copied … to avoid considerable unrelated dependencies".
- A SWIG Python module, `skymaps` (`src/swig_setup.i`), which also wraps
  parts of `astro` and `healpix`.

It does not provide any `gt*` tools or par files. `doc/release.notes`
(last tagged entry: `skymaps-02-17-00`, 2017) and `src/mainpage.h` are
historical. Commits since then only change the build.

## Position in the system

| Direction | Package | Why |
|-----------|---------|-----|
| Depends on | `astro` | `SkyDir`, `SkyFunction`, `SkyProj`, `Photon`, `PointingTransform`, `EarthCoordinate`. It gets these only through `st_facilities`'s PUBLIC link to `astro` (the library has no direct link). |
| Depends on | `healpix` | `Healpix`, `HealPixel`, `HealpixArray(IO)`, `Map`, `AlmOp`, `CosineBinner`, `alm_filter_tools`. |
| Depends on | `tip`, `facilities`, `st_facilities` | FITS I/O, `getDataPath`, `FileSys`/`Util`. |
| Depends on | `timeSystem` | `AbsoluteTime`, `BaryTimeComputer`, and related classes. |
| Depends on | `embed_python`, Python3 | `embed_python/Module.h` and the Python C API (`PySkyFunction`, `PythonUtilities`). |
| Depends on | CLHEP | `Hep3Vector` and related classes. Linked to the SWIG module as `CLHEP::VectorS`, and reaches the library transitively. |
| Used by | — | **Nothing in `src/` uses it.** There is no `#include "skymaps/..."` and no Python `import skymaps` outside this package. |

The tier is `instrument`. skymaps models the sky and the LAT response
(diffuse models, exposure, PSF) on top of tier-1 `io` packages, but its
dependency on `timeSystem` (tier 3) breaks the tier rule. If skymaps were
ever built again, it would belong in `analysis`.

## Key abstractions

See *Purpose and scope*. The central relationship is
`astro::SkyFunction` ← `SkyImage`, and `astro::SkyFunction` ←
`SkySpectrum` ← (diffuse, isotropic, composite, convolution, background
and PSF spectra). `BinnedPhotonData` is a `std::list<Band>` and is also a
`SkySpectrum`.

## Data flow and control flow

Not applicable in the Fermitools. Nothing invokes this code. In its
original use, Python scripts (pointlike) drove the SWIG module.

## Public interfaces

- **C++ headers:** `skymaps/*.h` (37 headers).
- **Python:** the SWIG module `skymaps` (`_skymaps`), from
  `src/swig_setup.i`.
- **Data:** `data/LivetimeCubeTemplate` and
  `data/aeff_P6_v1_diff_front.fits` (used by the test). `CMakeLists.txt`
  would install them to `share/fermitools/data/skymaps`, but
  `LivetimeCube.cxx` and the test look them up with
  `getDataPath("skymaps")`, which is `refdata/fermi/skymaps`, so they
  would not be found at run time.
- **Applications, par files:** none.

### Environment variables

| Variable | Read in | Effect |
|----------|---------|--------|
| `CALDB` | `IParams.cxx`, `EffectiveArea.cxx` (cached in a static on first use) | Root for hand-built `$CALDB/bcf/psf/...` and `$CALDB/bcf/ea/aeff_<irf>.fits` paths when no file name is given. Unset throws `std::invalid_argument`. |

`FERMI_INST_DIR` is used indirectly through `getDataPath`. Variables set
by conda activation are described under "Runtime environment" in
ScienceTools' `docs/architecture/README.md`.

## Build and test

`CMakeLists.txt` would be valid if the package were added back to
`FERMI_COMPONENTS`. It cannot be configured on its own, because it relies
on `FERMI_INSTALL_*` variables, `UseSWIG` settings and imported targets
from the top-level build. It defines:

- `skymaps` (SHARED) from 35 of the 36 `src/*.cxx`. `ComplexSkySpectrum.cxx`
  is omitted because its header includes ROOT's `TF1.h`.
- `swig_skymaps` (the Python module `_skymaps`).
- `test_skymaps` (`src/test/test_main.cxx`). It is not registered with
  CTest, but it is installed to `bin/`. It uses
  `getDataPath("skymaps")`, writes `test_philivetimecube.fits` to the
  current directory, and probably needs `CALDB` for the IRF classes.
  `recipe/tests/ST-unit-test` does not list it (see "Testing" in
  ScienceTools' `docs/architecture/README.md`), and no skymaps reference
  data exists anywhere.

**Would it build?** On 2026-10-22, against the current `RelWithDebInfo`
tree and the heads of all submodules, a `-fsyntax-only` compile of all 35
library sources and the test passed with no errors. Running `swig4.0` on
`swig_setup.i` and syntax-checking the generated wrapper also passed.
Linking, installing and running were **not** checked.

## Invariants and gotchas

- **Do not assume it is maintained.** Since 2017 the only changes are to
  the build: the CMake port, the ROOT exclusion, numpy removal (`fa6249e`)
  and the CLHEP target rename (`307d4cb`). Its IRF readers build CALDB
  paths by hand (`$CALDB/bcf/psf/psf_<name>_<class>_{front,back}.fits`,
  in `IParams.cxx`). Whether they match current CALDB contents has not
  been checked.
- **Duplicate GTI implementation.** `skymaps::Gti` is a copy of the
  ScienceTools GTI code. It is not `evtbin::Gti`, and the two have not
  been compared.
- **Bypasses `irfs`.** `EffectiveArea`/`IParams` open CALDB FITS files by
  path, and `IParams` throws if `CALDB` is unset. They do not use
  `irfInterface`/`irfLoader`.
- **Tier violation.** It depends on `timeSystem` (analysis) while its
  other dependencies are at tier 2 or lower. See *Position in the
  system*.

## Extension points

Not applicable. The package is not built. To revive it: add `skymaps` to
`FERMI_COMPONENTS`, decide what to do with `test_skymaps` (register it
with CTest, or stop installing it), and reconcile its tier.

## Further reading

- `src/mainpage.h` and `doc/release.notes` (historical).
- `docs/architecture/dependency-layers.md` ("Not built") in ScienceTools.
