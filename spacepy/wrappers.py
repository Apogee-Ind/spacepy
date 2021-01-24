# external imports
import numpy as np
import spiceypy as spice

# internal imports
from .objects import SpaceObject

def append_ssb_state(body: SpaceObject, t, frame):
    NAIF_id = body.id
    state = spice.spkssb(NAIF_id, t, frame)
    body.rvec = np.append(body.rvec, np.array([state[0:3]]), axis=0)
    body.vvec = np.append(body.vvec, np.array([state[3:6]]), axis=0)

def append_state(target: SpaceObject, obs: SpaceObject, t, frame, abcorr='NONE'):
    target_id = target.id
    obs_id = obs.id
    state, owlt = spice.spkez(target_id, t, frame, abcorr, obs_ic)
    