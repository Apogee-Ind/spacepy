# spacepy
Python module for spacecraft simulation and orbit determination

# Author Background
Hello. I am James Scales, Chief Science Officer at Apogee Industries and student of aerospace engineering at the University of Texas at Austin. As the administrator of Apogee Industries' GitHub account, I am using it to store personal projects such as these. When Apogee Industries is no more, I may move this repository to another account.
# Introduction
To my knowledge, no mainstream package for spacecraft simulation is present in common scientific python installs such as Anaconda. This package does not seek to fill that void, as doing so is far beyond my capabilities as a programmer! My goal with spacepy is twofold: I wish to incorporate features that I will find useful in academic and personal spaceflight projects, and I wish to consolidate my knowledge of spacecraft dynamics into a single, easily accessible Python package. Several of my (admittedly aspirational) goals for spacepy's features are as follows:

1. Provide an object-oriented framework for retrieving and using physical and orbit data of planets, moons, and other solar system bodies
2. Interface with the JPL Horizons system to acquire recent or possibly realtime solar system ephemerides
3. Provide an object-oriented framework for creating modifiable spacecraft with a range of attributes
4. Include routines for numerical integration of n-body dynamics and perturbations for which analytic formulations exist
5. Include functions for analysis of electrical, thermal, and environmental subsystems of spacecraft objects
6. Include routines for spacecraft attitude control, including simulation of perturbations

It should be noted that although I hope to include powerful numerical integration routines such as the multi-step Gauss-Jackson method, I lack the knowlegde to create C implementations of these algorithms (as would be done in numpy or scipy, for instance). Given this reality, it stands to reason that this package will not be as capable as the Aerospace Toolbox of MATLAB, or comparable programs. Do also note, however, that scipy contains no implementation of Gauss-Jackson, and whenever possible I will not reinvent the wheel and will use the undoubtedly faster algorithms of numpy and scipy.