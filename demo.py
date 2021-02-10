import numpy as np
from spacepy.tests import lowthrust, ephemeris, random
import spiceypy as spice

#spice.furnsh('spacepy/data/metakr.tm')
#ephemeris.main()
#lowthrust.main()

random.ceres_landing()
