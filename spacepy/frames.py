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
    return R3(a)@R1(b)@R3(c)

def R321(a, b, c):
    return R3(a)@R2(b)@R1(c)

