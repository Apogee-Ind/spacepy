import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import spacepy.data.constants as const
import spiceypy as spice

def to_deg(theta):
    """
    Convert an angle or array of angles from radians to degrees.

    Args:
    theta   dtype: float, ndarray   | Angle or ndarray of angles to convert. Units: rad

    Returns an angle of the same dtype as input, in .
    """
    return theta*(180/np.pi)

def to_rad(theta):
    """
    Convert an angle or array of angles from degrees to radians.

    Args:
    theta   dtype: float, ndarray   | Angle or ndarray of angles to convert. Units: deg

    Returns an angle of the same dtype as input, in radians.
    """
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

def unpack_geom(dims, shape):
    """
    Calculate the volume and surface area of a three-dimensional solid of given dimensions.

    Args:
    dims    dtype: tuple    | Tuple of critical dimensions. The number of required dimensions varies depending on the shape, see below. Units: m
    shape   dtype: str      | Must match one of: 'point', 'sphere', 'cylinder', 'cuboid', 'torus', 'cone'. Note that all shapes are assumed solid.

    Returns:
    V   dtype: float    | Volume of 3D solid. Units: m**3
    A   dtype: float    | Surface area of 3D solid. Units: m**2

    Explanation of shapes/dims:
    'point': dims may be empty, None, or contain any values. The returned volume and surface area are always 0.0
    'sphere': dims(0) is the radius
    'cylinder': dims(0) is the radius, dims(1) is the height
    'cuboid': dims(0), dims(1), dims(2) are the length, width, and height of the cuboid
    'torus': dims(0) is the major radius, dims(1) is the minor radius
    'cone': dims(0) is the radius, dims(1) is the height 
    Extra values in dims are ignored. 
    """
    shape_dict = {
        'point': {
            'V': lambda dims: 0.0,
            'A': lambda dims: 0.0
        },
        'sphere': {
            'V': lambda dims: (4/3)*np.pi * dims[0]**3,
            'A': lambda dims: 4*np.pi * dims[0]**2
        },
        'cylinder': {
            'V': lambda dims: np.pi * dims[0]**2 * dims[1],
            'A': lambda dims: 2*np.pi * dims[0] * (dims[0] + dims[1])
        },
        'cuboid': {
            'V': lambda dims: dims[0] * dims[1] * dims[2],
            'A': lambda dims: 2 * (dims[0]*dims[1] + dims[0]*dims[2] + dims[1]*dims[2])
        },
        'torus': {
            'V': lambda dims: 2*np.pi**2 * dims[0] * dims[1]**2,
            'A': lambda dims: 4*np.pi**2 * dims[0] * dims[1]
        },
        'cone': {
            'V': lambda dims: (1/3)*np.pi * dims[0]**2 * dims[1],
            'A': lambda dims: np.pi * dims[0] * (dims[0] + np.sqrt(dims[0]**2 + dims[1]**2))
        }
    }
    V = shape_dict[shape]['V'](dims)
    A = shape_dict[shape]['A'](dims)
    return V, A

