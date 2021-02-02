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
I have not yet implemented the `setup.py` file for spacepy. Spacepy is designed to work inside an Anaconda environment, as it requires many of the scientific packages contained including `numpy`, `scipy`, `matplotlib`, and `numba`. In addition, I am using the latest available versions of these packages during development. Spacepy also requires the package `spiceypy`, which can be installed on conda as follows:

`conda config --add conda-forge`

`conda install spiceypy`

Spiceypy is a wrapper for the C implementation of JPL's SPICE Toolkit, written by Andrew Annex and released under the MIT license. (Annex et al, (2020). SpiceyPy: a Pythonic Wrapper for the SPICE Toolkit. Journal of Open Source Software, 5(46), 2050, https://doi.org/10.21105/joss.02050). The source code for spiceypy can be found here:

https://github.com/AndrewAnnex/SpiceyPy

# Current Features
## Physical Data
Physical parameters for a number of bodies can be accessed via the `spacepy.objects.Sol`, `spacepy.objects.Planet`, `spacepy.objects.SmallBody`, and `spacepy.objects.Moon` classes. This data is taken from JPL and NASA, and for all bodies includes the gravitational parameter and radius. Some objects have more data, including flattening/zonal harmonics, atmospheric pressure, surface gravity and escape velocity, blackbody temperature, albedo, and more. See the docstring for each body's class for a full list. Currently, the following objects are supported:

1. The Sun
2. All 8 recognized major planets
3. Earth's Moon
4. Pluto
5. The asteroids Ceres, Vesta, Pallas, Psyche, Lutetia, Kleopatra, and Eros

Note that you will need to download the appropriate ephemeris files in order to use these bodies in solar system simulations. See the section on "Ephemeris Data" for more details.
## Numerical Propagation
Spacecraft trajectories can be integrated using a wrapper for `scipy.integrate.solve_ivp()` with method `DOP853`, an 8th-order variable-step Runge-Kutta method. Currently, nonsphericity (J2 effect) and thrust perturbations can be simulated.

A Kepler orbit solver using f and g functions may be implemented later. Spacepy currently has a class for using classical orbit elements, but this will likely change to a more robust implementation in future.

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

# Coming Soon

Here are some features that I aim to implement at some point, depending on how much time I have and when in particular I need their functionality:

1. Lambert problem solver for elliptic and potentially hyperbolic orbits (will likely use an algorithm from Vallado)
2. Using the above, a "porkchop" plot generator to analyze ballistic transfers between any planet or asteroid for which you have ephemerides
3. Also using the Lambert solver: a patched-conic trajectory solver, with potential use of an n-body integrator near SOI boundaries or other points
4. A `spacepy.objects.SpacecraftPart` class for creating thrusters, fuel tanks, solar panels, radiators, and more in a flexible way
5. Module for electrical subsystem analysis, with functions for solar panel and battery sizing based on cell specifications and power requirements
6. Module for thermal subsystem analysis, with functions for heat balancing, solar/albedo heating, and radiator sizing
