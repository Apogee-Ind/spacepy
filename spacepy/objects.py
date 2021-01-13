# scientific library imports
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

# internal imports
import spacepy.data.constants as const
from .data.planetdata import planets_orb, planets_phys
from .helpers import to_deg, to_rad, MA_to_nu, set_3daxes_equal, unpack_geom
from .frames import pqw2ijk

# root space object class
class SpaceObject():

    bodytype = 'unknown'
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

    bodytype = 'star'
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
    
    def to_rv(self):
        p = self.a*(1 - self.e**2)
        r = p/(1 + self.e*np.cos(self.nu))

        r_pqw = np.array([r*np.cos(self.nu), r*np.sin(self.nu), 0])
        v_pqw = np.array([-np.sqrt(self.gm/p)*np.sin(self.nu), np.sqrt(self.gm/p)*(self.e + np.cos(self.nu)), 0])

        r_ijk = np.array([pqw2ijk(r_pqw, self)])
        v_ijk = np.array([pqw2ijk(v_pqw, self)])
        r_pqw = np.array([r_pqw])
        v_pqw = np.array([v_pqw])

        return r_ijk, v_ijk, r_pqw, v_pqw
    
    def update_oe(self, r_ijk, v_ijk):
        # magnitudes
        r_mag = np.linalg.norm(r_ijk)
        v_mag = np.linalg.norm(v_ijk)
        # angular momentum vector (points normal to orbit plane)
        h_vec = np.cross(r_ijk, v_ijk)
        h_mag = np.linalg.norm(h_vec)
        z = np.array([0, 0, 1])
        # line of nodes (points in direction of ascending node)
        n = np.cross(z, h_vec)
        n_hat = n/np.linalg.norm(n)
        # eccentricity vector (points toward periapsis)
        e_vec = np.cross(v_ijk, h_vec)/self.gm - r_ijk/r_mag
        # update orbit elements
        # semi-major axis
        self.a = -self.gm*r_mag / (r_mag*v_mag**2 - 2*self.gm) # vis-viva eqn
        # eccentricity
        self.e = np.linalg.norm(e_vec)
        # inclination
        self.i = np.arccos(np.dot(h_vec, z) / h_mag)
        # argument of periapsis
        self.w = np.arccos(np.dot(n_hat, e_vec) / (np.linalg.norm(n_hat)*self.e))
        if e_vec[2] < 0.0:
            self.w = -self.w
        # longitude of ascending node
        self.lan = np.arctan2(n_hat[1], n_hat[0])
        # true anomaly
        self.nu = np.arccos(np.dot(r_ijk, e_vec) / (r_mag*self.e))
        if np.dot(r_ijk, v_ijk) < 0.0:
            self.nu = -self.nu


class SpaceCraft(SpaceObject):
    bodytype = 'spacecraft'

    def __init__(self, name, mass=1.0, shape='point', dims=(0.0)):
        super().__init__()
        self.name = name
        # spacecraft components
        self.parts = {'m': [], 'dims': [], 'shapes': [], 'V': [], 'A': []} # dictionary containing mass, dimensions, shape, volume, and area of each component
        self.m = np.sum(self.parts['m'])
        # add initial part
        self.add_part(mass, shape, dims)
        self.events = {}

    def add_part(self, mass, shape='sphere', dims=(1.0), name='generic_part'):
        V, A = unpack_geom(dims, shape)
        to_add = [mass, dims, shape, V, A]
        for key, prop in zip(self.parts, to_add):
            self.parts[key].append(prop)
        self.m = np.sum(self.parts['m'])

    def add_thruster(self, thrust, Isp, orientation=np.array([0, 1, 0])):
        """
        Orientation: unit vector defining direction of thrust when spacecraft is 
            pointed in the direction of the velocity vector (prograde), such that
            the thruster unit vector is expressed in the NTW coordinate frame.
        """
        self.thruster = {'max_thrust': thrust, 'Isp': Isp, 'orientation': orientation}
        self.m_fuel = 0.0

    def add_fuel(self, m_fuel):
        self.m_fuel = m_fuel
        self.m = self.m + self.m_fuel

    def place_in_orbit(self, around=Planet(), h_p=400.0, h_a=400.0, i=0.0, w=0.0, lan=0.0, nu=0.0):
        self.parent = around
        a = (2.0*self.parent.r + h_p + h_a)/2.0
        e = 1 - (self.parent.r + h_p)/a
        oe_vec = np.array([a, e, i, w, lan, nu])

        self.set_orbit(oe_vec, is_deg=True)

def create_LEO(h_p=400.0, h_a=400.0, i=0.0, w=0.0, lan=0.0, nu=0.0):
    Earth = Planet()
    spacecraft = SpaceObject()
    spacecraft.parent = Earth
    a = (2.0*Earth.r + h_p + h_a)/2.0
    e = 1 - (Earth.r + h_p)/a
    oe_vec = np.array([a, e, i, w, lan, nu])

    spacecraft.set_orbit(oe_vec, is_deg=True)
    return spacecraft