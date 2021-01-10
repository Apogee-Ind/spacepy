# scientific library imports
import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import animation

# internal imports
from .objects import SpaceObject, Sol, Planet, OrbitElements
from .helpers import to_deg, to_rad, MA_to_nu, set_3daxes_equal
import spacepy.data.constants as const

def _get_rhs(x: np.ndarray, sc: SpaceObject, j2=True):
    r_mag = np.linalg.norm(x[0:3])
    xd_twobody = x[3:6]
    xd_perturbing = []

    xdd_twobody = (-sc.parent.gm*x[0:3])/r_mag**3
    xdd_perturbing = []

    if (j2 & hasattr(sc.parent, 'j2')):
        xdd_j2 = np.array([
            -sc.parent.j2 * 1.5*(sc.parent.r/r_mag)**2 * (5*(x[2]/r_mag)**2 - 1),
            -sc.parent.j2 * 1.5*(sc.parent.r/r_mag)**2 * (5*(x[2]/r_mag)**2 - 1),
            sc.parent.j2 * 1.5*(sc.parent.r/r_mag)**3 * (3 - 5*(x[2]/r_mag)**2)
        ])
        xdd_j2 = xdd_twobody * xdd_j2
        xdd_perturbing.append(xdd_j2)
    
    xd_total = xd_twobody + np.sum(np.array(xd_perturbing), axis=0)
    xdd_total = xdd_twobody + np.sum(np.array(xdd_perturbing), axis=0)
    return np.concatenate((xd_total, xdd_total))

def rk8(sc: SpaceObject, tmax, tol=1e-6):
    if hasattr(sc, 'parent'):
        r0 = sc.rvec[-1]
        v0 = sc.vvec[-1]
        x0 = np.concatenate((r0, v0))
        def eqn(t, x):
            x_dot = _get_rhs(x, sc)
            return x_dot
        
        t_span = (sc.t[-1], tmax)
        output = integrate.solve_ivp(eqn, t_span, x0, method='DOP853', rtol=tol, atol=tol, max_step=120.0)
        sc.t = np.append(sc.t, output.t[1:])
        sc.rvec = np.append(sc.rvec, output.y[0:3,1:].T, axis=0)
        sc.vvec = np.append(sc.vvec, output.y[3:,1:].T, axis=0)
    else:
        raise AttributeError('Parent body is not defined.')

def plot_twobody(sc: SpaceObject):
    if hasattr(sc, 'parent'):
        assert (np.shape(sc.rvec)[0] > 1), 'No trajectory found. Please run an orbit simulation before plotting.'
        sc.parent.gen_ellipsoid()
        x_planet = sc.parent.ellipsoid_coords[0,:]
        y_planet = sc.parent.ellipsoid_coords[1,:]
        z_planet = sc.parent.ellipsoid_coords[2,:]

        sc.plot = plt.figure()
        ax = sc.plot.add_subplot(111, projection='3d')
        ax.plot_wireframe(x_planet, y_planet, z_planet, rcount=24, ccount=13)
        ax.plot(sc.rvec[1:,0], sc.rvec[1:,1], sc.rvec[1:,2], 'r-', alpha=0.5)
        ax.set_box_aspect([1,1,1]) # does not work here
        set_3daxes_equal(ax)
        ax.set_xlabel('X, km')
        ax.set_ylabel('Y, km')
        ax.set_zlabel('Z, km')
        ax.set_title('Trajectory Plot, ' + sc.parent.name + '-Centered Inertial Frame')
        plt.show()
    else:
        raise AttributeError('Parent body is not defined.')

def animate_twobody(sc: SpaceObject):
    if hasattr(sc, 'parent'):
        assert (np.shape(sc.rvec)[0] > 1), 'No trajectory found. Please run an orbit simulation before plotting.'
        sc.parent.gen_ellipsoid()
        x_planet = sc.parent.ellipsoid_coords[0,:]
        y_planet = sc.parent.ellipsoid_coords[1,:]
        z_planet = sc.parent.ellipsoid_coords[2,:]

        sc.plot = plt.figure()
        ax = sc.plot.add_subplot(111, projection='3d')
        ln, = ax.plot(sc.rvec[1,0], sc.rvec[1,1], sc.rvec[1,2], 'r-')
        frames = np.arange(2, np.shape(sc.t)[0], step=1)
        xdata = []
        ydata = []
        zdata = []

        def init():
            ax.plot_wireframe(x_planet, y_planet, z_planet, rcount=24, ccount=12)
            ax.set_box_aspect([1,1,1])
            set_3daxes_equal(ax)
            ax.set_xlabel('X, km')
            ax.set_ylabel('Y, km')
            ax.set_zlabel('Z, km')
            ax.set_title('Trajectory Plot, ' + sc.parent.name + '-Centered Inertial Frame')

        def update(i):
            xdata.append(sc.rvec[i,0])
            ydata.append(sc.rvec[i,1])
            zdata.append(sc.rvec[i,2])
            ln.set_data_3d(xdata, ydata, zdata)
            return ln
        
        fig_anim = animation.FuncAnimation(sc.plot, update, frames, init_func=init, interval=10, repeat_delay=1000)
        plt.show()
    else:
        raise AttributeError('Parent body is not defined.')
        


# class for gauss-jackson multi-step integrator
class _GJ8:
    # 8th-order Gauss-Jackson difference coefficients alpha(j,i)
    coef = np.array([
        []
    ])

    j_indices = np.array([-4, -3, -2, -1, 0, 1, 2, 3, 4, 5])

    def __init__(self, x0, tmax, step, tol):
        pass

# class to include f and g functions
class _Kepler:

    def __init__(self, x0):
        pass