# scientific library imports
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

# internal imports
import spacepy.data.constants as const
from .data.planetdata import planets_orb, planets_phys
from .helpers import to_deg, to_rad, MA_to_nu, set_3daxes_equal

# root space object class
class SpaceObject():

    def __init__(self):
        # position, velocity, and acceleration in inertial frame fixed to parent body
        self.rvec = np.zeros((1,3))
        self.vvec = np.zeros((1,3))
        self.avec = np.zeros((1,3))
        # position, velocity, and acceleration in orbit perifocal frame
        self.r_pf = np.zeros((1,3))
        self.v_pf = np.zeros((1,3))
        self.a_pf = np.zeros((1,3))
        # spacecraft time
        self.t = np.zeros(1)
        self.parent = None

    def set_orbit(self, oe_vec, is_deg=False, contains_MA=True):
        if hasattr(self, 'parent'):
            self.oe = OrbitElements(oe_vec, self.parent, is_deg, contains_MA)
            r_new, v_new, r_pf_new, v_pf_new = self.oe.to_rv()
            self.rvec = np.append(self.rvec, r_new, axis=0)
            self.vvec = np.append(self.vvec, v_new, axis=0)
            self.r_pf = np.append(self.r_pf, r_pf_new, axis=0)
            self.v_pf = np.append(self.v_pf, v_pf_new, axis=0)
        else:
            raise AttributeError('Parent body is not defined.')

# derivative class for the Sun
class Sol(SpaceObject):

    def __init__(self):
        super().__init__()
        self.name = 'Sol'
        self.system = 'Sol'
        self.parent = None
        self.type = 'star'
        self.gm = 1.3271244e11 # km^3 s^-2
        self.r = 695700 # km
        self.m = 1.9885e30 # kg
        self.rho = 1408 # kg/m**3
        self.sday = 25.05*86400
        self.vesc = 617.7 # km/s
        self.g = 274 # m/s**2
        self.j2 = 0.0
        self.f = 9.0e-6
        self.p = 86.8 # Pa
        self.alb = 0.0
        self.Tsurf = 5777 # K
        self.Teff = 5777 # K

# derivative class for planets
class Planet(Sol):
    """
        Create a new Planet object.
        Args:
        name    Type: string. Must match one of: 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', or 'Neptune'. If not specified, defaults to 'Earth'.

        Attributes:
        gm      Gravitational Parameter. Units: km^3 * s^-2
        r       Equatorial radius. Unit: km
        rmean   Mean radius. Unit: km
        m       Mass. Unit: kg
        rho     Bulk density. Unit: kg * m^-3
        sday    Sidereal rotation period. Unit: s
        syr     Sidereal revolution period. Unit: yr (Earth year)
        vesc    Escape velocity. Unit: km * s^-1
        g       Equatorial surface gravity. Unit: m * s^-2
        j2      J2 zonal harmonic term. Unit: dimensionless
        f       Ellipsoid flattening term. Unit: dimensionless
        p       Equatorial sea-level/datum atmospheric pressure. Unit: Pa
        alb     Bond albedo. Unit: dimensionless
        Tsurf   Surface temperature. Unit: K
        Teff    Effective blackbody temperature. Unit: K

        Note: Not all bodies have data for all attributes. In this case, the attribute is not defined.
        References: JPL Planetary and Lunar Ephemerides DE430 & DE431
    """
    bodytype = "major_planet"

    def __init__(self, name="Earth"):
        super().__init__()
        if type(name) != str:
            raise TypeError
        else:
            self.name = name
            self.type = 'major_planet'
            self.system = 'Sol'
            self.parent = Sol()
        if name in planets_phys:
            for key in planets_phys[name]:
                setattr(self, key, planets_phys[name][key])
        else:
            for key in planets_phys['Earth']:
                setattr(self, key, planets_phys['Earth'][key])
        if name in planets_orb:
            oe_vec = np.array(list(planets_orb[name].values()))
            self.oe = OrbitElements(oe_vec, self.parent, is_deg=True, contains_MA=True)

    def gen_ellipsoid(self):
        # equatorial radius a, polar radius c
        a = self.r
        c = a*(1 - self.f)

        # spherical angles
        u = np.linspace(0, 2*np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        # conversion to cartesian
        x = a*np.outer(np.cos(u), np.sin(v))
        y = a*np.outer(np.sin(u), np.sin(v))
        z = c*np.outer(np.ones_like(u), np.cos(v))

        self.ellipsoid_coords = np.array([x, y, z])

# class to contain orbit elements
class OrbitElements():

    def __init__(self, oe_vec, body: SpaceObject, is_deg=False, contains_MA=False):

        if is_deg:
            oe_vec[2:] = to_rad(oe_vec[2:])
        if contains_MA:
            oe_vec[5] = MA_to_nu(oe_vec[5], oe_vec[1])

        self.a = oe_vec[0]
        self.e = oe_vec[1]
        self.i = oe_vec[2]
        self.w = oe_vec[3]
        self.lan = oe_vec[4]
        self.nu = oe_vec[5]

        self.gm = body.gm

    def R313(self):
        R_w = np.array([[np.cos(self.w),   -np.sin(self.w),    0],
                        [np.sin(self.w),   np.cos(self.w),     0],
                        [0,            0,              1]
                        ])
        R_i = np.array([[1,     0,          0         ],
                        [0,    np.cos(self.i),  -np.sin(self.i)],
                        [0,    np.sin(self.i),  np.cos(self.i) ]
                        ])
        R_lan = np.array([[np.cos(self.lan),     -np.sin(self.lan),   0],
                        [np.sin(self.lan),     np.cos(self.lan),    0],
                        [0,               0,              1]
                        ])

        return R_lan@R_i@R_w
    
    def to_rv(self):
        p = self.a*(1 - self.e**2)
        r = p/(1 + self.e*np.cos(self.nu))

        r_pqw = np.array([r*np.cos(self.nu), r*np.sin(self.nu), 0])
        v_pqw = np.array([-np.sqrt(self.gm/p)*np.sin(self.nu), np.sqrt(self.gm/p)*(self.e + np.cos(self.nu)), 0])

        r_ijk = np.array([self.R313()@r_pqw])
        v_ijk = np.array([self.R313()@v_pqw])
        r_pqw = np.array([r_pqw])
        v_pqw = np.array([v_pqw])

        return r_ijk, v_ijk, r_pqw, v_pqw