# scientific library imports
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import spiceypy as spice

# internal imports
import spacepy.data.constants as const
from .data.bodydata import planet_data, planet_systems, smallbody_data, moon_data
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
        self.id = 10
        self.system = 'Sol'
        self.parent = None
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
        name            dtype: str  | Must match one of: 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', or 'Neptune'. If not specified, defaults to 'Earth'.
        barycenter      dtype: bool | Specifies whether this planet should have the NAIF integer code corresponding to its system's barycenter.

        Initial physical attributes are generated from values stored in spacepy.data.bodydata.planet_data. This dictionary contains the following fields for all planets:
        id      dtype: int      | NAIF integer code for this body. Unit: dimensionless
        gm      dtype: float    | Gravitational parameter. Unit: km**3 / s**2
        r       dtype: float    | Equatorial radius. Unit: km
        rmean   dtype: float    | Mean radius. Unit: km
        m       dtype: float    | Mass. Unit: kg
        rho     dtype: float    | Bulk density. Unit: kg / m**3
        sday    dtype: float    | Sidereal rotation period. Unit: s
        syr     dtype: float    | Sidereal revolution period. Unit: yr
        mag     dtype: float    | V-band magnitude. Unit: dimensionless
        alb_g   dtype: float    | Geometric albedo. Unit: dimensionless
        g       dtype: float    | Equatorial surface gravity. Unit: m / s**2
        vesc    dtype: float    | Escape velocity. Unit: km / s**2
        j2      dtype: float    | 2nd zonal harmonic coefficient. Unit: dimensionless

        In addition, all inner planets have the following physical parameters defined:
        f       dtype: float    | Flattening. Unit: dimensionless
        p       dtype: float    | Surface atmospheric pressure. Unit: Pa
        Tsurf   dtype: float    | Surface temperature. Unit: K
        Tbb     dtype: float    | Blackbody temperature. Unit: K
        alb_b   dtype: float    | Bond albedo. Unit: dimensionless

        In addition, all outer planets have the following physical parameters defined:
        j4      dtype: float    | 4th zonal harmonic coefficient. Unit: dimensionless

    """
    bodytype = "major_planet"

    def __init__(self, name='Earth', barycenter=False):
        super().__init__()
        if type(name) != str:
            raise TypeError
        else:
            self.name = name
            self.system = 'Sol'
            self.parent = Sol()

            if name in planet_data:
                for key in planet_data[name]:
                    setattr(self, key, planet_data[name][key])
            else:
                for key in planet_data['Earth']:
                    setattr(self, key, planet_data['Earth'][key])
            if barycenter:
                if name in planet_systems:
                    setattr(self, 'id', planet_systems[name])

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

class SmallBody(Sol):
    """
    Create a new SmallBody object. These objects include dwarf planets, asteroids, and comets.

    Args:
    name    dtype: str  | Must match one of: 'Pluto', 'Ceres', 'Vesta', 'Pallas', 'Psyche', 'Lutetia', 'Kleopatra', 'Eros'. Defaults to 'Ceres' if unspecified.

    Initial physical attributes are generated from values stored in spacepy.data.bodydata.smallbody_data. The following are defined for all small bodies:
    id      dtype: int      | NAIF integer code for this body.
    gm      dtype: float    | Gravitational parameter. Unit: km**3 / s**2
    r       dtype: float    | Equatorial radius. (Note: for asteroids, this value is the radius of an equivalent spherical body.) Unit: km
    sday    dtype: float    | Sidereal rotaion period. Unit: s
    alb_g   dtype: float    | Geometric albedo. Unit: dimensionless

    In addition, the dwarf planet Pluto has the following physical parameters defined:
    rmean   dtype: float    | Mean radius. Unit: km
    rho     dtype: float    | Bulk density. Unit: kg / m**3
    mag     dtype: float    | V-band magnitude: Unit: dimensionless
    g       dtype: float    | Equatorial surface gravity. Unit: m / s**2
    vesc    dtype: flaot    | Escape velocity. Unit: km / s

    In addition, some asteroids have one or more of the following physical parameters defined:
    extent  dtype: ndarray  | Dimensions of an equivalent triaxial ellipsoid. Unit: km
    rho     dtype: float    | Bulk density. Unit: kg / m**3  
    """
    bodytype = 'small_body'

    def __init__(self, name='Ceres'):
        super().__init__()
        if type(name) != str:
            raise TypeError
        else:
            self.name = name
            self.system = 'Sol'
            self.parent = Sol()

            if name in smallbody_data:
                for key in smallbody_data[name]:
                    setattr(self, key, smallbody_data[name][key])
            else:
                for key in smallbody_data['Ceres']:
                    setattr(self, key, smallbody_data['Ceres'][key])

class Moon(Sol):
    """
    Create a new Moon object. Moons are natural satellites of planets or minor bodies.

    Args:
    name    dtype: str  | Must match one of: 'Moon'. Defaults to 'Moon' if not specified.

    Initial physical attributes are generated from values stored in spacepy.data.bodydata.moon_data. For Earth's Moon, the following parameters are defined:
    id      dtype: int      | NAIF integer code for this moon.
    parent  dtype: str      | Body around which this moon orbits.
    gm      dtype: float    | Gravitational parameter. Unit: km**3 / s**2
    r       dtype: float    | Equatorial radius. Unit: km
    rmean   dtype: float    | Mean radius. Unit: km
    rho     dtype: float    | Bulk density. Unit: kg / m**3
    sday    dtype: float    | Sidereal rotation period. Unit: s
    """
    bodytype = 'moon'

    def __init__(self, name='Moon'):
        super().__init__()
        assert (type(name) == str), 'Name must be a string.'
        self.name = name

        if name in moon_data:
            for key in moon_data[name]:
                setattr(self, key, moon_data[name][key])
            else:
                for key in moon_data['Moon']:
                    setattr(self, key, moon_data['Moon'][key])
        
        self.parent = Planet(self.system)


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
    """
    Create a new SpaceCraft object.
    Args:
    name    dtype: str      | Name for this spacecraft. Can be whatever you like.
    id      dtype: int      | User-defined id value for this spacecraft. Can be anything, but conventionally is a negative integer. Defaults to -1 if not specified.
    mass    dtype: float    | Initial spacecraft dry mass. Defaults to 1.0 if not specified. Unit: kg
    shape   dtype: str      | Must be one of: 'point', 'sphere', 'cylinder', 'cuboid', 'torus', or 'cone'. Defaults to 'point' if not specified.
    dims    dtype: tuple    | Contains critical dimensions as required by spacepy.helpers.unpack_geom(). Defaults to (0.0) if not specified. Unit: m

    Attributes:
    bodytype    dtype: str  | Has the value 'spacecraft' for all SpaceCraft objects.
    parts       dtype: dict | Dictionary with keys: 'm', 'dims', 'shapes', 'V', 'A'. Each field contains a list of values for each part, and is updated as parts are added.
    m           dtype: float| Total mass of spacecraft, including fuel if fuel is added. Updates dynamically during simulation.
    events      dtype: dict | Contains events generated by the integrator or other user-defined processes. Events are timestamped and include state vector information.
    """
    bodytype = 'spacecraft'

    def __init__(self, name, id=-1, mass=1.0, shape='point', dims=(0.0)):
        super().__init__()
        self.name = name
        self.id = id
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

class System:
    bodytype = 'system'

    def __init__(self, epoch, *bodies):
        self.contents = {
            'star':{},
            'major_planet':{},
            'small_body':{},
            'moon':{},
            'spacecraft':{}
        }
        for body in bodies:
            self.contents[body.bodytype][body.id] = body
        self.epoch = epoch
        self.et_start = spice.str2et(epoch)
    
    def add_body(self, body: SpaceObject):
        self.contents[body.bodytype][body.id] = body

    def _gen_ICRF_vectors(self, stop, step, include_moons=True):
        assert (10 in self.contents['star']), 'System must contain the Sun in order to use ICRF coordinates.'
        self.epoch_end = stop
        self.et_end = spice.str2et(stop)
        self.step_time = step

        t_iterable = np.arange(self.et_start, self.et_end, step)
        for t in t_iterable:
            state = spice.spkssb(10, t, 'J2000')
            sun = self.contents['star'][10]
            sun.rvec = np.append(sun.rvec, np.array([state[0:3]]), axis=0)
            sun.vvec = np.append(sun.vvec, np.array([state[3:6]]), axis=0)
            if self.contents['major_planet']:
                for pid in self.contents['major_planet']:
                    state = spice.spkssb(pid, t, 'J2000')
                    planet = self.contents['major_planet'][pid]
                    planet.rvec = np.append(planet.rvec, np.array([state[0:3]]), axis=0)
                    planet.vvec = np.append(planet.vvec, np.array([state[3:6]]), axis=0)
            if self.contents['small_body']:
                for sbid in self.contents['small_body']:
                    state = spice.spkssb(sbid, t, 'J2000')
                    sbody = self.contents['small_body'][sbid]
                    sbody.rvec = np.append(sbody.rvec, np.array([state[0:3]]), axis=0)
                    sbody.vvec = np.append(sbody.vvec, np.array([state[3:6]]), axis=0)
            if (include_moons & bool(self.contents['moon'])):
                for mid in self.contents['moon']:
                    state = spice.spkssb(mid, t, 'J2000')
                    mbody = self.contents['moon'][mid]
                    mbody.rvec = np.append(mbody.rvec, np.array([state[0:3]]), axis=0)
                    mbody.vvec = np.append(mbody.vvec, np.array([state[3:6]]), axis=0)