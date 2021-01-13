# external imports
import numpy as np

# internal imports
from .helpers import to_deg, to_rad, MA_to_nu

def rot1(theta):
    # rotation about x-axis
    R = np.array([
        [1,                 0,                  0               ],
        [0,                 np.cos(theta),      -np.sin(theta)  ],
        [0,                 np.sin(theta),      np.cos(theta)   ]
    ])
    return R

def rot2(theta):
    # rotation about y-axis
    R = np.array([
        [np.cos(theta),     0,                  np.sin(theta)   ],
        [0,                 1,                  0,              ],
        [-np.sin(theta),    0,                  np.cos(theta)   ]
    ])
    return R

def rot3(theta):
    # rotation about z-axis
    R = np.array([
        [np.cos(theta),     -np.sin(theta),     0               ],
        [np.sin(theta),     np.cos(theta),      0               ],
        [0,                 0,                  1               ]
    ])
    return R

def R313(a, b, c):
    return rot3(a)@rot1(b)@rot3(c)

def R321(a, b, c):
    return rot3(a)@rot2(b)@rot1(c)

def rsw2ijk(x_rsw, oe):
    # spacecraft-fixed Radial, Transverse, Normal frame to body-centered inertial frame
    u = oe.nu + oe.w
    x_ijk = R313(oe.lan, oe.i, u)@x_rsw
    return x_ijk

def pqw2ijk(x_pqw, oe):
    # perifocal frame to body-centered inertial frame
    x_ijk = R313(oe.lan, oe.i, oe.w)@x_pqw
    return x_ijk

def rsw2pqw(x_rsw, oe):
    return rot3(oe.nu)@x_rsw

def ntw2ijk(x_ntw, oe):
    den = np.sqrt(1 + 2*oe.e*np.cos(oe.nu) + oe.e**2)
    fpa = np.arctan2((oe.e*np.sin(oe.nu))/den, (1 + oe.e*np.cos(oe.nu))/den)
    x_rsw = rot3(-fpa)@x_ntw
    x_ijk = rsw2ijk(x_rsw, oe)
    return x_ijk