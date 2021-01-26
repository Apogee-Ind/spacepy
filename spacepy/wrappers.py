# external imports
import numpy as np
import spiceypy as spice

# internal imports
from .objects import SpaceObject

def append_ssb_state(body: SpaceObject, t_iterable, frame):
    state = np.zeros((np.shape(t_iterable)[0],6))
    NAIF_id = body.id
    for t, nt in zip(t_iterable, range(np.shape(t_iterable)[0])):
        state[nt] = spice.spkssb(NAIF_id, t, frame)
    body.state[0]['state'] = np.append(body.state[0]['state'], np.array(state), axis=0)
    body.state[0]['t'] = np.append(body.state[0]['t'], t_iterable, axis=0)
    

def append_state(target: SpaceObject, obs: SpaceObject, t, frame, abcorr='NONE'):
    target_id = target.id
    obs_id = obs.id
    state, owlt = spice.spkez(target_id, t, frame, abcorr, obs_ic)
    target.state[obs.id]['state'] = np.append(target.state[obs.id]['state'], np.array([state]), axis=0)
    target.state[obs.id]['t'] = np.append(target.state[obs.id]['t'], t, axis=0)
    