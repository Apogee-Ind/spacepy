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
        ## will be deprecated soon
        self.rvec = np.zeros((1,3))
        self.vvec = np.zeros((1,3))
        self.avec = np.zeros((1,3))

        # position, velocity, and acceleration in orbit perifocal frame
        ## will be deprecated soon
        self.r_pf = np.zeros((1,3))
        self.v_pf = np.zeros((1,3))
        self.a_pf = np.zeros((1,3))

        # generic state vector container, for storing position and velocity in multiple frames
        self.state = {}
        # structure is as follows:
        #self.state = {
        #    # keys must be the NAIF id code of the body/location that coordinates are expressed relative to. 0 is the solar system barycenter.
        #    0:{
        #        'epoch':'epoch_str', # human-readable initial epoch of the object, with the format 'YYYY MMM DD HH:MM:SS'
        #        'step':86400, # step size for integration and/or plotting, in seconds
        #        'frame':'ECLIPJ2000' # SPICE-compatible string code for the reference frame in which to express the state vector
        #        't':np.zeros(n) # ndarray of time values, represents time since the J2000 epoch, in seconds
        #        'state':np.zeros((n,6)) # ndarray of state vector. Columns 0-2 are position, columns 3-5 are velocity, in km and km/s
        #    },
        #    # 399 is Earth's center of mass
        #    399:{
        #        'epoch':'epoch_str',
        #        'step':60,
        #        'frame':'IAU_EARTH'
        #        't':np.zeros(n)
        #        'state':np.zeros((1,6))
        #    }
        #}
        # for n-body simulations, each body should contain vectors relative to the other bodies, all expressed in a frame centered on the system barycenter.
        # for example: an n-body heliocentric simulation might include the following bodies: Sun, Earth, Jupiter; with spacecraft coordinates expressed in ECLIPJ2000 frame relative to solar system barycenter.
        
        # time
        ## will be deprecated eventually
        self.epoch = None
        self.t = np.zeros(1)
        self.parent = None

    def _add_state(self, epoch, step, obs_id=0, frame='ECLIPJ2000', state=np.zeros((1,6)), t=np.zeros(1)):
        state_keys = ('epoch', 'step', 'frame', 'state', 't')
        state_vals = (epoch, step, frame, state, t)

        self.state[obs_id] = {k:v for k, v in zip(state_keys, state_vals)}

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
    Create a new SpaceCraft instance.
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
    flags       dtype: dict | Dictionary of flags, each flag is an ndarray of 0's or 1's, tracking some property over the time span of integration.

    A SpaceCraft instance has certain simulation routines not available to other SpaceObject classes. These include the ability to define parts, surfaces, thrusters, and unique simulation events. Additional attributes and methods for attitude dynamics, the electrical subsystem, and the thermal subsystem will be added in a later update. SpaceCraft instances will contain a unique state vector, with observer ID equal to that of the SpaceCraft instance, for purposes of storing the three attitude angles and corresponding rates. 
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
        self.flags = {}

    def add_part(self, mass, shape='sphere', dims=(1.0), name='generic_part'):
        """
        Add a new part to a SpaceCraft instance.

        Args:
        mass    dtype: float    | Total mass of the part. Unit: kg
        shape   dtype: str      | Shape of the part, must match one of: 'point', 'sphere', 'cuboid', 'cylinder', 'cone', 'torus'. Defaults to 'sphere' if not specified. Note that all part shapes are assumed solid; thin-walled structures will be added later.
        dims    dtype: tuple    | Contains critical dimensions as required by spacepy.helpers.unpack_geom(). Defaults to (0.0) if not specified. Unit: m
        name    dtype: str      | User-given name for the part. May be anything, although if another part of the same name exists, the new part will overwrite the old one.

        Calling spacepy.objects.SpaceCraft.add_part() updates the SpaceCraft instance's mass automatically. Note that at the moment, all parts of a SpaceCraft instance are assumed to exist at the instance's center of mass for purposes of simulation, and are rigid bodies. This will change in a later update, and add_part() will require a user to specify the attachment point(s) and orientation of the new part. Potentially, later updates could add non-rigid body dynamics to joints between parts of a SpaceCraft instance. In that case, calling add_part() will also require a user to specify spring and damping constants for the joint(s) between the new part and other parts.
        """
        V, A = unpack_geom(dims, shape)
        to_add = [mass, dims, shape, V, A]
        for key, prop in zip(self.parts, to_add):
            self.parts[key].append(prop)
        self.m = np.sum(self.parts['m'])

    def add_thruster(self, thrust, Isp, orientation=np.array([0, 1, 0])):
        """
        Add a thruster to a SpaceCraft instance.

        Args:
        thrust          dtype: float    | Maximum thrust of the thruster. Unit: N
        Isp             dtype: float    | Specific impulse of the thruster, measured using Earth's standard gravity. Unit: s
        orientation     dtype: ndarray  | ndarray of shape (3,) containing a vector describing the direction of thrust in the inertial, non-rotating NTW frame. Defaults to [0, 1, 0] (a thruster pointing in the direction of the velocity vector). Note that this orientation is expressed with respect to the spacecraft's velocity, not its attitude (roll/pitch/yaw). This is subject to change when attitude dynamics modeling is added to SpaceCraft instances, at which point calling add_thruster() will require specification of an orientation vector in a rotating, non-inertial body fixed frame, optionally with a position vector in the same frame specifying the offset of the thruster from the SpaceCraft instance's center of mass.

        In a later update, thrusters will be merged with other parts in a new class, spacepy.objects.SpacecraftPart, which will not be a subclass of SpaceObject as parts will not have uniquely defined state vectors (in other words, parts cannot detach from their parent SpaceCraft instance). This will also allow for the creation of multiple thrusters, something currently not supported by add_thruster() (a new thruster will overwrite an existing thruster).
        """
        self.thruster = {'max_thrust': thrust, 'Isp': Isp, 'orientation': orientation}
        self.m_fuel = 0.0

    def get_max_dV(self):
        assert hasattr(self, 'thruster'), 'SpaceCraft instance must have a thruster to calculate maximum delta-V.'
        self.dV = self.thruster['Isp'] * const.g_0 * np.log(self.m/(self.m - self.m_fuel))
        return self.dV

    def add_fuel(self, m_fuel):
        """
        Add fuel to a SpaceCraft instance.

        Args:
        m_fuel  dtype: float    | Mass of fuel to add. Unit: kg

        Calling add_fuel() automatically updates the SpaceCraft instance's total mass. When a simulation including thrust is run, the spacecraft's thruster will burn fuel at a rate determined by the thrust equation (thrust and mass are assumed constant over an integration time step). When the SpaceCraft instance's fuel mass reaches 0.0, an event with key 'out_of_fuel' will be created. This event can be used in plotting routines from spacepy.simulation.

        In a later update, fuel mass will be an attribute of a SpacecraftPart of type 'fuel_tank'. This will allow for separate fuel tanks, as may be required by multiple thrusters using different fuel types.
        """
        self.m_fuel = m_fuel
        self.m = self.m + self.m_fuel
        self.get_max_dV()

    def place_in_orbit(self, around=Planet(), h_p=400.0, h_a=400.0, i=0.0, w=0.0, lan=0.0, nu=0.0):
        self.parent = around
        a = (2.0*self.parent.r + h_p + h_a)/2.0
        e = 1 - (self.parent.r + h_p)/a
        oe_vec = np.array([a, e, i, w, lan, nu])

        self.set_orbit(oe_vec, is_deg=True)
    

def create_LEO(h_p=400.0, h_a=401.0, i=0.0, w=0.0, lan=0.0, nu=0.0):
    """
    Create a SpaceCraft instance and place it in Low Earth Orbit. WARNING: Will be deprecated in a future update.

    Args:
    h_p     dtype: float    | Geometric height (altitude) of orbit perigee above Earth's surface. Defaults to 400.0 if not specified. Unit: km
    h_a     dtype: float    | Geometric height (altitude) of orbit apogee above Earth's surface. Defaults to 401.0 if not specified. Unit: km
    i       dtype: float    | Inclination of orbit plane with respect to Earth's equator. Defaults to 0.0 if not specified. Unit: degree; Range: 0 to 180
    w       dtype: float    | Argument of perigee. Defaults to 0.0 if not specified. Unit: degree; Range: 0 to 360
    lan     dtype: float    | Longitude of the ascending node. Defaults to 0.0 if not specified. Unit: degree; Range: 0 to 360
    nu      dtype: float    | True anomaly (angle between spacecraft's position and perigee, measured positive in the direction of motion). Defaults to 0.0 if not specified. Unit: degree; Range: -90 to +90

    Note that very low eccentricities and inclinations (< 0.001) create near-singular orbit elements, which can cause undesirable effects during numerical integration as certain angles will change rapidly. create_LEO() and the spacepy.objects.OrbitElements class do not currently support alternative orbit elements for circular or equatorial orbits.

    This function will be modified in a future update, as two-body dynamics will use the spacepy.objects.System() class along with three- and n-body dynamics. In future, this function will create a System instance containing a SpaceCraft instance and a Planet instance, with classical orbit elements passed to the SpaceCraft instance and greater support for near-singular elements.
    """
    Earth = Planet()
    spacecraft = SpaceObject()
    spacecraft.parent = Earth
    a = (2.0*Earth.r + h_p + h_a)/2.0
    e = 1 - (Earth.r + h_p)/a
    oe_vec = np.array([a, e, i, w, lan, nu])

    spacecraft.set_orbit(oe_vec, is_deg=True)
    return spacecraft


from .wrappers import append_ssb_state
class System:
    bodytype = 'system'

    def __init__(self, epoch, *bodies: SpaceObject):
        """
        Create a new System instance.

        Args:
        epoch       | dtype: str                            | Start epoch, formatted as 'YYYY MMM DD HH:MM:SS'
        *bodies     | dtype: spacepy.objects.SpaceObject    | Bodies to include in the simulation. Must be an instance of subclasses Sol, Planet, Moon, SmallBody, or SpaceCraft
        """
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
    
    def remove_body(self, body: SpaceObject):
        try:
            removed = self.contents[body.type].pop(body.id)
        except KeyError:
            print(f'System already does not contain {body.name}')

    def _add_state_to_all(self, obs_id, frame, epoch, step):
        for bodytype, bodies in self.contents.items():
            for bid, body in bodies.items():
                body._add_state(epoch, step, obs_id, frame)

    def _gen_ssb_vectors(self, stop, step, include_moons=True, frame='ECLIPJ2000'):
        """
        Generate state vectors for all bodies in System, relative to the solar system barycenter and between the defined times.

        Args:
        stop            dtype: str      | Simulation stop time, formatted as 'YYYY MMM DD HH:MM:SS'
        step            dtype: float    | Desired interval between subsequent values in the state vector, in seconds  
        include_moons   dtype: bool     | Whether to include moons contained in System in the simulation. Defaults to True if not specified.
        frame           dtype: str      | Reference frame in which to express output state vectors. Defaults to 'ECLIPJ2000' if not specified.
        """
        assert (10 in self.contents['star']), 'System must contain the Sun in order to use coordinates centered on the Solar System Barycenter.'
        self.epoch_end = stop
        self.et_end = spice.str2et(stop)
        self.step_time = step

        self._add_state_to_all(0, frame, self.epoch, step)

        t_iterable = np.arange(self.et_start, self.et_end, step)
        
        for bodytype, bodies in self.contents.items():
            if bodies:
                for bid, body in bodies.items():
                    append_ssb_state(body, t_iterable, frame)