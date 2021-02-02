# external imports
import numpy as np 
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import spiceypy as spice

# internal imports
from .objects import SpaceObject, Sol, SpaceCraft, System, OrbitElements
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

class Lambert:
    """
    Solver for the multiple-revolution Lambert problem.
    """
    def __init__(self, Body1, Body2, epoch_start, epoch_end):
        pass

    def solve_elliptic(r0: np.ndarray, r1: np.ndarray, method=1):
        assert isinstance(method, int), 'Method must be either the integer 1 (for Type I transfers) or 2 (for Type II transfers).'
        if method == 2:
            t_m = -1
        else:
            t_m = 1
        cos_delta_nu = np.dot(r0, r1) / (np.linalg.norm(r0)*np.linalg.norm(r1))
