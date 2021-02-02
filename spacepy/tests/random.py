import numpy as np

from ..objects import Planet, SmallBody, SpaceCraft, Moon
import spacepy.data.bodydata as bodydata
import spacepy.data.constants as const

def ss_dv():
    earth = Planet()
    v_circ = np.sqrt(earth.gm/(earth.r + 350)) # km/s

    ss_Isp = 380.0
    m_payload = 100e3
    m_dry = 120e3 + m_payload
    m_fuel = 1200e3
    m_wet = m_dry + m_fuel

    max_dv = ss_Isp * const.g_0 * np.log(m_wet/m_dry) # m/s

    v_p = v_circ + max_dv/1000 # km/s
    vesc = np.sqrt((2*earth.gm)/(earth.r + 350.0))
    c3 = v_p**2 - vesc**2
    max_vinf = np.sqrt(c3)

    target_vinf = 6.6 # km/s
    vp_required = np.sqrt(vesc**2 + target_vinf**2)
    dv_required = np.abs(vp_required - v_circ)
    dv_leftover = max_dv/1000 - dv_required

    m_after = m_wet /  np.exp(dv_required*1000 / (ss_Isp*const.g_0))

    #print(f'Fully loaded & fueled starship has {np.round(max_dv/1000,3)} km/s delta-V')
    print(f'Starship max C3: {c3}, max v_inf: {max_vinf}')
    print(f'Target V_inf of {target_vinf} km/s: {np.round(dv_leftover,3)} km/s delta-V leftover, total mass {np.round(m_after/1000)} mT after departure.')

    ceres = SmallBody()
    vesc_ceres = np.sqrt((2*ceres.gm)/(ceres.r + 200.0))
    vinf_ceres = 6.3
    vp_ceres = np.sqrt(vesc_ceres**2 + vinf_ceres**2)
    v_circ_ceres = np.sqrt(ceres.gm/(ceres.r + 200.0))
    vp_min = vp_ceres - dv_leftover
    print(f'ceres hyperbolic periapsis velocity: {np.round(vp_ceres,3)} and circular orbit velocity at same altitude: {np.round(v_circ_ceres,3)}')
    print(f'minimum v_p after expending all fuel: {np.round(vp_min,3)} km/s')
    a_final = ceres.gm*(ceres.r + 200.0) / (2*ceres.gm - (ceres.r + 200.00)*vp_min**2)
    #print(f'final orbit semi-major axis after expending all fuel decelerating at Ceres: {a_final} km')
    """ ceres = SmallBody('Ceres')
    v_ceres = np.sqrt(ceres.gm/(ceres.r + 110.0))
    omega_c = (2*np.pi)/ceres.sday
    v_rot = ceres.r * omega_c
    print(f'200 km circular low Ceres orbit velocity: {np.round(v_ceres,3)}')
    v_ascent = v_ceres #- v_rot
    v_descent = v_ascent * 1.25
    print(f'Landing delta-V, assuming in direction of Ceres rotation and 25% margin for hovering: {np.round(v_descent,3)}')
    print(f'Ascent delta-V, assuming launching at equator in direction of rotation, to 200 km orbit: {np.round(v_ascent,3)}')

    moon = Moon()
    v_moon = np.sqrt(moon.gm/(moon.r + 110.0))
    v_d_moon = v_moon * 1.25
    print(f'Moon landing delta-V: {np.round(v_d_moon,3)} km/s, ascent delta-V: {np.round(v_moon,3)} km/s')
    """