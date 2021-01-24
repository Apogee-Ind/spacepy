# scientific library imports
import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import animation

# internal imports
from .objects import SpaceObject, Sol, Planet, SmallBody, Moon, OrbitElements, SpaceCraft, System
from .helpers import to_deg, to_rad, MA_to_nu, set_3daxes_equal
from .frames import ntw2ijk
import spacepy.data.constants as const

def _update_mass(t, x, thrust, thrust_spec, sc: SpaceCraft):
    g0 = const.g_0
    dt = t - sc.t_old
    dm = (thrust*1000*dt)/(sc.thruster['Isp']*g0)
    if (sc.m_fuel - dm) <= 0:
        thrust_spec[:,1] = 0
        sc.events['out of fuel'] = {'t': t, 'rvec': x[0:3], 'vvec': x[3:6]}
        return False
    else:
        sc.m_fuel = sc.m_fuel - dm
        sc.m = sc.m - dm
        return True

def _thrust_accel(t, x: np.ndarray, sc: SpaceCraft, thrust_spec, do_mass_update):
    throttle = thrust_spec[thrust_spec[:,0] <= t, 1][-1]
    thrust_mag = throttle*sc.thruster['max_thrust']/1000 # convert to kN
    xdd_thrust_ntw = (thrust_mag*sc.thruster['orientation'])/sc.m
    xdd_thrust = ntw2ijk(xdd_thrust_ntw, sc.oe)
    if do_mass_update:
        do_mass_update = _update_mass(t, x, thrust_mag, thrust_spec, sc)

    return xdd_thrust

def _twobody_accel(t, x: np.ndarray, sc: SpaceCraft, j2):
    xdd_twobody = (-sc.parent.gm*x[0:3])/r_mag**3
    if (j2 & hasattr(sc.parent, 'j2')):
        xdd_j2 = np.array([
            -sc.parent.j2 * 1.5*(sc.parent.r/r_mag)**2 * (5*(x[2]/r_mag)**2 - 1),
            -sc.parent.j2 * 1.5*(sc.parent.r/r_mag)**2 * (5*(x[2]/r_mag)**2 - 1),
            sc.parent.j2 * 1.5*(sc.parent.r/r_mag)**3 * (3 - 5*(x[2]/r_mag)**2)
        ])
        xdd_j2 = xdd_twobody * xdd_j2
    
    return xdd_twobody + xdd_j2

def _nbody_accel(t, x: np.ndarray, sc: SpaceCraft, sys: System):
    pass   

def rk8(sc: SpaceObject, tmax, thrust_spec=None, j2=False, drag=False, **integrator_options):
    if hasattr(sc, 'parent'):
        r0 = sc.rvec[-1]
        v0 = sc.vvec[-1]
        x0 = np.concatenate((r0, v0))
        sc.x_old = x0
        sc.t_old = sc.t[-1]
        do_mass_update = True

        def _get_rhs(t, x: np.ndarray, sc: SpaceObject, thrust_spec, j2, drag, do_mass_update):
            r_mag = np.linalg.norm(x[0:3])
            xd_twobody = x[3:6]
            # update stored orbital elements - necessary for perturbation calculations
            sc.oe.update_oe(x[0:3], xd_twobody)

            xd_perturbing = []

            # define two-body acceleration
            xdd_twobody = (-sc.parent.gm*x[0:3])/r_mag**3
            xdd_perturbing = []

            # spacecraft thrust
            if type(thrust_spec) != type(None):
                throttle = thrust_spec[thrust_spec[:,0] <= t, 1][-1]
                thrust_mag = throttle*sc.thruster['max_thrust']/1000 # convert to kN
                xdd_thrust_ntw = (thrust_mag*sc.thruster['orientation'])/sc.m
                xdd_thrust = ntw2ijk(xdd_thrust_ntw, sc.oe)
                if do_mass_update:
                    do_mass_update = _update_mass(t, x, thrust_mag, thrust_spec, sc)

                xdd_perturbing.append(xdd_thrust)

            # parent body j2 zonal harmonic (nonsphericity)
            if (j2 & hasattr(sc.parent, 'j2')):
                xdd_j2 = np.array([
                    -sc.parent.j2 * 1.5*(sc.parent.r/r_mag)**2 * (5*(x[2]/r_mag)**2 - 1),
                    -sc.parent.j2 * 1.5*(sc.parent.r/r_mag)**2 * (5*(x[2]/r_mag)**2 - 1),
                    sc.parent.j2 * 1.5*(sc.parent.r/r_mag)**3 * (3 - 5*(x[2]/r_mag)**2)
                ])
                xdd_j2 = xdd_twobody * xdd_j2
                xdd_perturbing.append(xdd_j2)
            
            # sum all terms of velocity and acceleration
            xd_total = xd_twobody + np.sum(np.array(xd_perturbing), axis=0)
            xdd_total = xdd_twobody + np.sum(np.array(xdd_perturbing), axis=0)

            # keep current value to use as old value on next iteration
            sc.x_old = x
            sc.t_old = t
            return np.concatenate((xd_total, xdd_total))
            del xd_perturbing
            del xdd_perturbing
        
        t_span = (sc.t[-1], tmax)
        output = integrate.solve_ivp(_get_rhs, t_span, x0, method='DOP853', args=(sc, thrust_spec, j2, drag, do_mass_update), **integrator_options)
        sc.t = np.append(sc.t, output.t[1:])
        sc.rvec = np.append(sc.rvec, output.y[0:3,1:].T, axis=0)
        sc.vvec = np.append(sc.vvec, output.y[3:,1:].T, axis=0)
        print(f'Number of rhs evaluations: {output.nfev}')
        print(output.message)
    else:
        raise AttributeError('Parent body is not defined.')

def plot_twobody(sc: SpaceCraft):
    if hasattr(sc, 'parent'):
        assert (np.shape(sc.rvec)[0] > 1), 'No trajectory found. Please run an orbit simulation before plotting.'
        sc.parent.gen_ellipsoid()
        x_planet = sc.parent.ellipsoid_coords[0,:]
        y_planet = sc.parent.ellipsoid_coords[1,:]
        z_planet = sc.parent.ellipsoid_coords[2,:]

        sc.plot = plt.figure()
        ax = sc.plot.add_subplot(111, projection='3d')
        ax.plot_wireframe(x_planet, y_planet, z_planet, rcount=24, ccount=13)
        ax.plot(sc.rvec[1:,0], sc.rvec[1:,1], sc.rvec[1:,2], 'r-')
        for event in sc.events:
            coords = sc.events[event]['rvec']
            ax.scatter(coords[0], coords[1], coords[2], 'bx')
            ax.text(coords[0], coords[1], coords[2], event)
        ax.set_box_aspect([1,1,1])
        set_3daxes_equal(ax)
        ax.set_xlabel('X, km')
        ax.set_ylabel('Y, km')
        ax.set_zlabel('Z, km')
        plt.legend([sc.name, sc.parent.name])
        ax.set_title('Trajectory Plot, ' + sc.parent.name + '-Centered Inertial Frame')
        plt.show()
    else:
        raise AttributeError('Parent body is not defined.')

def animate_twobody(sc: SpaceCraft):
    if hasattr(sc, 'parent'):
        assert (np.shape(sc.rvec)[0] > 1), 'No trajectory found. Please run an orbit simulation before plotting.'
        sc.parent.gen_ellipsoid()
        x_planet = sc.parent.ellipsoid_coords[0,:]
        y_planet = sc.parent.ellipsoid_coords[1,:]
        z_planet = sc.parent.ellipsoid_coords[2,:]

        sc.plot = plt.figure()
        ax = sc.plot.add_subplot(111, projection='3d')
        ln, = ax.plot(sc.rvec[1,0], sc.rvec[1,1], sc.rvec[1,2], 'r-')
        ln.set_label(sc.name + ' trajectory')
        dot, = ax.plot(sc.rvec[1,0], sc.rvec[1,1], sc.rvec[1,2], 'o')
        dot.set_label(sc.name)
        frames = np.arange(2, np.shape(sc.t)[0], step=1)
        xdata = []
        ydata = []
        zdata = []

        def init():
            wire = ax.plot_wireframe(x_planet, y_planet, z_planet, rcount=24, ccount=12)
            wire.set_label(sc.parent.name)
            ax.set_box_aspect([1,1,1])
            set_3daxes_equal(ax)
            ax.set_xlabel('X, km')
            ax.set_ylabel('Y, km')
            ax.set_zlabel('Z, km')
            ax.set_title('Trajectory Plot, ' + sc.parent.name + '-Centered Inertial Frame')
            plt.legend()

        def update(i):
            xdata.append(sc.rvec[i,0])
            ydata.append(sc.rvec[i,1])
            zdata.append(sc.rvec[i,2])
            ln.set_data_3d(xdata, ydata, zdata)
            dot.set_data_3d(sc.rvec[i,0], sc.rvec[i,1], sc.rvec[i,2])
            return ln
        
        fig_anim = animation.FuncAnimation(sc.plot, update, frames, init_func=init, interval=10, repeat_delay=1000)
        #fig_anim.save('plot.mp4')
        plt.show()
    else:
        raise AttributeError('Parent body is not defined.')
        
def plot_system(sys: System, do_markers=True):
    sys.plot = plt.figure()
    ax = sys.plot.add_subplot(111, projection='3d')

    for bodytype in sys.contents:
        for bid, body in sys.contents[bodytype].items():
            line, = ax.plot(body.rvec[2:,0], body.rvec[2:,1], body.rvec[2:,2])
            line.set_label(body.name)
            if do_markers:
                ax.plot(body.rvec[1,0], body.rvec[1,1], body.rvec[1,2], 'ok')
                ax.plot(body.rvec[-1,0], body.rvec[-1,1], body.rvec[-1,2], 'xk')

            if bodytype == 'spacecraft':
                for event in body.events:
                    coords = body.events[event]['rvec']
                    ax.scatter(coords[0], coords[1], coords[2], 'bx')
                    ax.text(coords[0], coords[1], coords[2], event)

    ax.set_box_aspect([1,1,1])
    set_3daxes_equal(ax)
    ax.set_xlabel('X, km')
    ax.set_ylabel('Y, km')
    ax.set_zlabel('Z, km')
    plt.legend()
    plt.suptitle('Trajectory Plot of System in heliocentric coordinate frame.')
    plt.title(sys.epoch + ' to ' + sys.epoch_end)
    plt.show()

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