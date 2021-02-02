# external imports
import numpy as np 
from numpy.polynomial import Polynomial
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import spiceypy as spice
from numba import jit

# internal imports
from .objects import SpaceObject, Sol, SpaceCraft, OrbitElements
from .helpers import to_deg, to_rad
import spacepy.data.constants as const
import spacepy.data.bodydata as bodydata
    
def v_escape(Body: Sol, altitude):
    r = Body.r + altitude
    return np.sqrt(2*Body.gm / r)

def v_circular(Body: Sol, altitude):
    r = Body.r + altitude
    return np.sqrt(Body.gm/r)

def vis_viva(Body: Sol, r, a):
    return np.sqrt(Body.gm * (2/r - 1/a))





#@jit(nopython=True)
def battin_method(r0: np.ndarray, r1: np.ndarray, t, gm, method=1, tol=1e-6):
    """
    Implements the method for a single-revolution Lambert solver described by Battin (1987:325-342), using the algorithm as presented by Vallado (2013:494-497).
    """
    def kepler_f(a, r, E):
        return 1 - (a/r) * (1 - np.cos(E))

    def kepler_g(a, E, t, gm):
        return t - np.sqrt(a**3 / gm) * (E - np.sin(E))

    def kepler_fh(a, r, H):
        return 1 - (a/r) * (1 - np.cosh(H))

    def kepler_gh(a, H, t, gm):
        return t - np.sqrt(-a**3 / gm) * (np.sinh(H) - H)
    #assert isinstance(method, int), 'Method must be either the integer 1 (for Type I transfers) or 2 (for Type II transfers).'
    if method == 2: # determine transfer method coefficient
        t_m = -1
    else:
        t_m = 1
    # angle between radius vectors (change in true anomaly)
    r0_mag = np.linalg.norm(r0)
    r1_mag = np.linalg.norm(r1)
    cos_delta_nu = np.dot(r0, r1) / (r0_mag*r1_mag)
    sin_delta_nu = t_m * np.sqrt(1 - cos_delta_nu**2)
    delta_nu = np.arctan2(sin_delta_nu, cos_delta_nu)
    #print(delta_nu)

    c = np.sqrt(r0_mag**2 + r1_mag**2 - 2*r0_mag*r1_mag*cos_delta_nu)
    s = 0.5 * (r0_mag + r1_mag + c)
    epsilon = (r1_mag - r0_mag) / r0_mag

    cos_dn2 = np.cos(delta_nu*0.5)
    cos2_dn4 = np.cos(delta_nu*0.25)**2
    sin2_dn4 = np.sin(delta_nu*0.25)**2

    r_ratio = r1_mag/r0_mag
    tan2_2w = (0.25*epsilon**2) / (np.sqrt(r_ratio) + r_ratio*(2 + np.sqrt(r_ratio)))
    r_op = np.sqrt(r1_mag*r0_mag) * (cos2_dn4 + tan2_2w)
    print(r_op)

    if delta_nu < np.pi:
        l = (sin2_dn4 + tan2_2w) / (sin2_dn4 + tan2_2w + cos_dn2)
    elif delta_nu < 2*np.pi:
        l = (cos2_dn4 + tan2_2w - cos_dn2) / (cos2_dn4 + tan2_2w)
    #print(l)
    m = (gm*t**2) / (8*r_op**3)

    f_eta = lambda x: x / (np.sqrt(x + 1) + 1)**2
    xi_coef = np.zeros(50)
    for n in range(50):
        xi_coef[n] = n**2 / ((2*n)**2 - 1)
    f_xi = lambda x, eta: (8*(np.sqrt(1 + x) + 1)) / (3 + 1 / (5 + eta + (9/7)*eta / cont_frac(4, eta)))
    def cont_frac(n, eta):
        if n < 50:
            return 1 + (xi_coef[n]*eta) / cont_frac(n+1, eta)
        else:
            return 1

    u_coef = np.array([1/3, 4/27, 8/27, 208/891, 340/1287])
    f_b = lambda h1, h2: 27*h2 / (4*(1 + h1)**3)
    f_u = lambda B: B / (2*(np.sqrt(1 + B) + 1))
    f_k = lambda U: u_coef[0] / (1 + (u_coef[1]*U) / (1 + (u_coef[2]*U) / (1 + (u_coef[3]*U) / (1 + u_coef[4]*U))))

    # how tf to know if its elliptic or not??
    x = 0.0
    while True:
        eta = f_eta(x)
        xi = f_xi(x, eta)
        if xi==np.nan:
            raise ValueError('something broke')
        if eta==np.nan:
            raise ValueError('something broke')
        #print(x)
        #print(xi)
        #print(eta)

        h1 = ((l + x)**2 * (1 + 3*x + xi)) / ((1 + 2*x + l) * (4*x + xi) * (3 + x))
        h2 = (m * (x - l + xi)) / ((1 + 2*x + l) * (4*x + xi) * (3 + x))
        #p = Polynomial([-h2, 0, (-1 - h1), 1])
        #roots = p.roots()
        #y = roots[0]
        b = f_b(h1, h2)
        u = f_u(b)
        k = f_k(u)
        y = ((1 + h1)/3) * (2 + np.sqrt(1 + b) / (1 + 2*u*k**2))

        x_new = np.sqrt((0.5*(1 - l))**2 + m/y**2) - 0.5*(1 + l)
        if np.abs(x_new - x) < tol:
            break
        x = x_new
    
    a = (gm*t**2) / (16*r_op*x_new*y**2)

    if a > 0.0: # for elliptic orbits
        sin_beta2 = np.sqrt((s - c)/ (2*a))
        beta_e = np.arcsin(sin_beta2)*2
        if delta_nu > np.pi:
            beta_e = -beta_e
        # minimum time and semi-major axis
        a_min = 0.5*s
        t_min = np.sqrt(a_min**3 / gm) * (np.pi - beta_e + np.sin(beta_e))

        sin_alpha2 = np.sqrt(s / (2*a))
        alpha_e = np.arcsin(sin_alpha2)*2
        if t > t_min:
            alpha_e = 2*np.pi - alpha_e
        delta_E = alpha_e - beta_e

        # f and g functions
        f = kepler_f(a, r0_mag, delta_E)
        g = kepler_g(a, delta_E, t, gm)
        gdot = kepler_f(a, r1_mag, delta_E)
    else: # hyperbolic orbits
        sinh_alpha2 = np.sqrt(s / (-2*a))
        alpha_h = np.arcsinh(sinh_alpha2)*2
        sinh_beta2 = np.sqrt((s - c) / (-2*a))
        beta_h = np.arcsinh(sinh_beta2)*2

        delta_H = alpha_h - beta_h

        # hyperbolic f and g functions
        f = kepler_fh(a, r0_mag, delta_H)
        g = kepler_gh(a, delta_H, t, gm)
        gdot = kepler_fh(a, r1_mag, delta_H)
    
    # velocity at departure and arrival
    v0 = (r1 - f*r0) / g
    v1 = (gdot*r1 - r0) / g
    #print(v0)
    #print(v1)
    return v0, v1
