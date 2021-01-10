import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from scipy import integrate

# planet physical and orbit data
from ..data.planetdata import planets_phys, planets_orb

# physical constants
import spacepy.data.constants as const

# helper functions
def to_deg(theta):
    return theta*(180/np.pi)

def to_rad(theta):
    return theta*(np.pi/180)

def MA_to_nu(ma, e, order=3, is_deg=False):
    if is_deg:
        ma = to_rad(ma)
    terms = np.array([ma,
                      (2*e - 0.25*e**3)*np.sin(ma),
                      1.25*e**2*np.sin(2*ma),
                      (13/12)*e**3*np.sin(3*ma)])
    nu = np.sum(terms[:order + 1])
    return nu

def set_3daxes_equal(ax: plt.Axes):
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d()
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5*np.max(np.abs(limits[:,1] - limits[:,0]))

    ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
    ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
    ax.set_zlim3d([origin[2] - radius, origin[2] + radius])

# space object master class
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

    def twobody_sim(self, tmax, step, tol=1e-6):
        if hasattr(self, 'parent'):
            if tmax <= self.t[-1]:
                raise ValueError('Simulation end time cannot be in the past.')
            r0 = self.rvec[-1]
            v0 = self.vvec[-1]
            x0 = np.concatenate((r0, v0))
            def eqn(t, x):
                r_mag = np.linalg.norm(x[0:3])
                if hasattr(self.parent, 'j2'):
                    j2_term_xy = self.parent.j2 * 1.5*(self.parent.r/r_mag)**2 * (5*(x[2]/r_mag)**2 - 1)
                    j2_term_z = self.parent.j2 * 1.5*(self.parent.r/r_mag)**3 * (3 - 5*(x[2]/r_mag)**2)
                else:
                    j2_term_xy = 0.0
                    j2_term_z = 0.0
                x_dot = np.array([x[3],
                                x[4],
                                x[5],
                                ((-self.parent.gm*x[0])/r_mag**3)*(1 - j2_term_xy),
                                ((-self.parent.gm*x[1])/r_mag**3)*(1 - j2_term_xy),
                                ((-self.parent.gm*x[2])/r_mag**3)*(1 + j2_term_z)
                                ])
                return x_dot
            
            t_span = (self.t[-1], tmax)
            t_eval = np.arange(self.t[-1], tmax, step=step)
            output = integrate.solve_ivp(eqn, t_span, x0, method='RK45', t_eval=t_eval, rtol=tol, atol=tol)
            self.t = np.append(self.t, output.t[1:])
            self.rvec = np.append(self.rvec, output.y[0:3,1:].T, axis=0)
            self.vvec = np.append(self.vvec, output.y[3:,1:].T, axis=0)
        else:
            raise AttributeError('Parent body is not defined.')

    def encke(self, tmax, step, tol=1e-4):
        if hasattr(self, 'parent'):
            if tmax <= self.t[-1]:
                raise ValueError('Simulation end time cannot be in the past.')
        

    def plot_trajectory(self, is_3d=True):
        if hasattr(self, 'parent'):
            if is_3d:
                self.parent.gen_ellipsoid()
                x_planet = self.parent.ellipsoid_coords[0,:]
                y_planet = self.parent.ellipsoid_coords[1,:]
                z_planet = self.parent.ellipsoid_coords[2,:]
                self.plot = plt.figure()
                ax = self.plot.add_subplot(111, projection='3d')
                ax.plot_wireframe(x_planet, y_planet, z_planet, rcount=24, ccount=13)
                ax.plot(self.rvec[1:,0], self.rvec[1:,1], self.rvec[1:,2], 'r-')
                ax.set_box_aspect([1,1,1])
                set_3daxes_equal(ax)
                ax.set_xlabel('X, km')
                ax.set_ylabel('Y, km')
                ax.set_zlabel('X, km')
                ax.set_title('Trajectory Plot, ' + self.parent.name + '-Centered Inertial Frame')
                plt.show()
            else:
                raise NotImplementedError
        else:
            raise AttributeError('Parent body is not defined.')

    def dV(self, *maneuvers, interval=1):
        pass

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
        self.Tsurf = 5777 # K
        self.Teff = 5777 # K

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
        c = a
        if hasattr(self, 'f'):
            c = a*(1 - self.f)

        # spherical angles
        u = np.linspace(0, 2*np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        # conversion to cartesian
        x = a*np.outer(np.cos(u), np.sin(v))
        y = a*np.outer(np.sin(u), np.sin(v))
        z = c*np.outer(np.ones_like(u), np.cos(v))

        self.ellipsoid_coords = np.array([x, y, z])

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

def create_LEO(h_p=400.0, h_a=400.0, i=0.0, w=0.0, lan=0.0):
    Earth = Planet()
    spacecraft = SpaceObject()
    spacecraft.parent = Earth
    a = (2.0*Earth.r + h_p + h_a)/2.0
    e = 1 - (Earth.r + h_p)/a
    oe_vec = np.array([a, e, i, w, lan, 0.0])

    spacecraft.set_orbit(oe_vec, is_deg=True)
    return spacecraft

def create_GEO(i=0.0, w=0.0, lan=0.0):
    Earth = Planet()
    spacecraft = SpaceObject()
    spacecraft.parent = Earth
    a, e = 42164, 0.0
    oe_vec = np.array([a, e, i, w, lan, 0.0])

    spacecraft.set_orbit(oe_vec, is_deg=True)
    return spacecraft