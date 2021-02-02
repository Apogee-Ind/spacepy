# external imports
import numpy as np
import spiceypy as spice

# internal imports
from .objects import SpaceObject

def append_ssb_state(body: SpaceObject, t_iterable: np.ndarray, frame):
    """
    Add an array of state vector values to body.state[0]['state'] and time values to body.state[0]['t'], representing the body's position and velocity relative to the solar system barycenter at a range of times.

    Args:
    body        dtype: SpaceObject  | Body must have an integer NAIF ID recognized by the SPICE API, and for which an SPK ephemeris file has been loaded using spacepy.data.bodydata.load().
    t_iterable  dtype: ndarray      | Array of time values, expressed as seconds since the J2000 epoch, at which to generate a state vector.
    frame       dtype: str          | Reference frame recognizeable by the SPICE API, in which to express the state vector.

    Calling append_ssb_state() may return an error from spiceypy if required ephemeris files are missing, or if the reference frame is not recognized. Users should check carefully to ensure they have loaded the necessary ephemeris file for the body, frame, and times at which a state vector is desired.

    append_ssb_state() wraps spiceypy.spkssb(), documentation for which can be found on that module's wiki.
    """
    state = np.zeros((np.shape(t_iterable)[0],6))
    NAIF_id = body.id
    for t, nt in zip(t_iterable, range(np.shape(t_iterable)[0])):
        state[nt] = spice.spkssb(NAIF_id, t, frame)
    body.state[0]['state'] = np.append(body.state[0]['state'], np.array(state), axis=0)
    body.state[0]['t'] = np.append(body.state[0]['t'], t_iterable, axis=0)
    

def append_state(target: SpaceObject, obs_id, t_iterable: np.ndarray, frame, abcorr='NONE'):
    """
    Add an array of state vector and time values to target.state[obs_id], representing the target's position and velocity relative to the body with NAIF ID code obs_id at a range of times.

    Args:
    target      dtype: SpaceObject  | Target body must have an integer NAIF ID recognized by the SPICE API, and for which an SPK ephemeris file has been loaded using spacepy.data.bodydata.load().
    obs_id      dtype: int          | NAIF ID code of the body to which state vectors are to be expressed. If not specified, defaults to 10 (the Sun).
    t_iterable  dtype: ndarray      | Array of time values, expressed as seconds since the J2000 epoch, at which to generate state vectors. Unit: s
    frame       dtype: str          | Reference frame recognizeable by the SPICE API, in which to express the state vector.
    abcorr      dtype: str          | Tells spiceypy.spkez() which type of aberration correction, if any, should be applied to the state vector. If not specified, defaults to 'NONE'.

    append_state() wraps spiceypy.spkez(), documentation for which can be found on that module's wiki. 
    """
    state = np.zeros((np.shape(t_iterable)[0],6))
    owlt = np.zeros_like(t_iterable)
    target_id = target.id
    for t, nt in zip(t_iterable, range(np.shape(t_iterable)[0])):
        state[nt], owlt[nt] = spice.spkez(target_id, t, frame, abcorr, obs_id)
    target.state[obs_id]['state'] = np.append(target.state[obs_id]['state'], np.array(state), axis=0)
    target.state[obs_id]['t'] = np.append(target.state[obs_id]['t'], t_iterable, axis=0)
    