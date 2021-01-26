# spacepy
Python module for spacecraft simulation and orbit determination

# Author Background
Hello. I am James Scales, Chief Science Officer at Apogee Industries and student of aerospace engineering at the University of Texas at Austin. As the administrator of Apogee Industries' GitHub account, I am using it to store personal projects such as these. When Apogee Industries is no more, I may move this repository to another account.
# Introduction
To my knowledge, no mainstream package for spacecraft simulation is present in common scientific python installs such as Anaconda. This package does not seek to fill that void, as doing so is far beyond my capabilities as a programmer! My goal with spacepy is twofold: I wish to incorporate features that I will find useful in academic and personal spaceflight projects, and I wish to consolidate my knowledge of spacecraft dynamics into a single, easily accessible Python package. Several of my goals for spacepy's features are as follows:

1. Object-oriented modeling of solar system bodies, complete with their physical data, orbit characteristics, and relation to each other.

2. Numerical propagation of spaceraft orbits in two-body, restricted three-body, n-body, and patched conic scenarios, with modeling of perturbations included where possible.

3. Analysis of spacecraft subsystems including electrical, propulsion, thermal, and attitude control over a given trajectory.

# Installation and Use
I have not yet implemented the `setup.py` file for spacepy. Spacepy is designed to work inside an Anaconda environment, as it requires many of the scientific packages contained including `numpy`, `scipy`, and `matplotlib`. In addition, I am using the latest available versions of these packages during development. Spacepy also requires the package `spiceypy`, which can be installed on conda as follows:

`conda config --add conda-forge`

`conda install spiceypy`

Spiceypy is a wrapper for the C implementation of JPL's SPICE Toolkit, written by Andrew Annex and released under the MIT license. (Annex et al, (2020). SpiceyPy: a Pythonic Wrapper for the SPICE Toolkit. Journal of Open Source Software, 5(46), 2050, https://doi.org/10.21105/joss.02050). The source code for spiceypy can be found here:

https://github.com/AndrewAnnex/SpiceyPy

# Current Features
## Physical Data
Physical parameters for a number of bodies can be accessed via the `spacepy.objects.Sol`, `spacepy.objects.Planet`, `spacepy.objects.SmallBody`, and `spacepy.objects.Moon` classes. This data is taken from JPL and NASA. Currently, the following objects are supported:

1. The Sun
2. All 8 recognized major planets
3. Earth's Moon
4. Pluto
5. The asteroids Ceres, Vesta, Pallas, Psyche, Lutetia, Kleopatra, and Eros

Note that you will need to download the appropriate ephemeris files in order to use these bodies in solar system simulations. See the section on "Ephemeris Data" for more details.
## Numerical Propagation
Spacecraft trajectories can be integrated using a wrapper for `scipy.integrate.solve_ivp()` with method `DOP853`, an 8th-order variable-step Runge-Kutta method. Currently, nonsphericity (J2 effect) and thrust perturbations can be simulated.

## Ephemeris Data
Using integration with spiceypy, spacepy supports simulation of all solar system bodies, provided necessary SPK kernels are installed. These files tend to be quite large, so are not included in this repository. The required SPK kernels can be found on the NAIF website here: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/
Spacepy requires the following ephemeris files:

`de440.bsp` (general-purpose solar system ephemeris, contains data for the barycenters of all planets, the Sun, Earth, and the Moon)

In order to use bodies other than the barycenters of the 8 planets, the Earth, the Moon, and the Sun, you will need additional ephemerides:

`mar097.bsp` (ephemeris for the Martian system)

`jup310.bsp` (Jovian system)

`sat427.bsp` (Saturnian system)

`ura111.bsp` (Uranian system)

`nep095.bsp` (Neptunian system)

I have included, for convenience, an ephemeris file which contains data for asteroids Ceres, Vesta, Pallas, and Eros. This ephemeris file, along with the others listed above, are required by a metakernel file which spacepy will attempt to load using the function `spacepy.data.bodydata.load()` which in turn calls `spiceypy.furnsh()`. If you do not wish to download the extra ephemeris files (other than the general-purpose), you will need to remove the filenames from `spacepy/data/metakr.tm` before calling `spacepy.data.bodydata.load()`.

Place the SPK files in `spacepy/data/kernels/spk` before using spacepy. 

# TODO
## Short-term
1. ~~Add attrubute under Planet class for the integer NAIF object code~~
2. Add `spacepy.objects.System` class to handle creation and simulation of heliocentric and Earth-Moon systems
3. Add wrappers for spiceypy time conversion functions 
4. Add physical and ellipsoid data for the Moon
5. Document existing functions and add import statements in package `__init.py__`

## Mid-term
1. Implement restricted three-body equations of motion
2. Implement SOI boundary detection for patched-conic modeling
3. ~~Add data for asteroid 1 Ceres~~

## Long-term
1. Add three-axis spacecraft attitude to EOM state vector along with perturbing torques
2. Add attitude control law simulation
3. Add thermal and electrical subsystem simulation
