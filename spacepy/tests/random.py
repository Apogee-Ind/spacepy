import numpy as np

from ..objects import Sol, Planet, SmallBody, SpaceCraft, Moon
import spacepy.data.bodydata as bodydata
import spacepy.data.constants as const
import spacepy.orbits as orb

def ss_dv():
    earth = Planet()
    mars = Planet('Mars')
    ceres = SmallBody()
    sun = Sol()

    r_LEO = earth.r + 350.0
    r_e = const.au
    r_m = 230269198.0
    r_c = 413964389.0

    e_m_dep, e_m_arr = orb.hohmann_dv(sun, r_e, r_m)
    e_c_dep, e_c_arr = orb.hohmann_dv(sun, r_e, r_c)
    m_c_dep, m_c_arr = orb.hohmann_dv(sun, r_m, r_c)
    e_m_dv = orb.capture_dv(earth, e_m_dep, r_LEO, 0.0, capture=False)
    e_c_dv = orb.capture_dv(earth, e_c_dep, r_LEO, 0.0, capture=False)
    m_c_dv = orb.capture_dv(mars, m_c_dep, 350+mars.r, 0.0, capture=False)

    print(f'dv to leave: earth-mars {e_m_dv} km/s, earth-ceres {e_c_dv} km/s, mars-ceres {m_c_dv} km/s')
    print('note: leaving from 350 km altitude circular orbit at earth, and from 350 km circular orbit at mars')

    e_m_cap = orb.capture_dv(mars, e_m_arr, 350+mars.r, 0.0)
    e_c_cap = orb.capture_dv(ceres, e_c_arr, 200+ceres.r, 0.95)
    m_c_cap = orb.capture_dv(ceres, m_c_arr, 200+ceres.r, 0.95)

    print(f'dv to capture: earth-mars {e_m_cap} km/s, earth-ceres {e_c_cap} km/s, mars-ceres {m_c_cap} km/s')
    print('note: Mars capture is to 350 km altitude circular orbit, ceres capture is to 200 km orbit with e=0.95')

def ceres_landing():
    ceres = SmallBody()
    sun = Sol()

    r_parking = ceres.r + 200.0
    v_park = orb.v_circular(ceres, 200.0)
    omega_ceres = (2*np.pi)/ceres.sday
    v_rot = omega_ceres*ceres.r

    v_ascend = v_park
    v_descend = v_ascend * 1.15
    v_hover = v_descend - v_ascend

    v_ascend_min = v_ascend - v_rot
    v_descend_min = v_ascend_min * 1.15
    v_hover_min = v_descend_min - v_ascend_min

    lander = SpaceCraft('lm', mass=3000)
    lander.add_thruster(16000, 310)
    lander.add_fuel(2353.0)
    lander_dv = lander.get_max_dV()

    dv_required = (v_ascend + v_descend) * 1.10
    n_landings = lander_dv/dv_required

    print(f'ceres dV to land: {v_descend} and take off: {v_ascend} to/from 200 km parking orbit\nassuming 15% margin: {v_hover} for hovering')
    print(f'absolute minimum dv to land (using ceres rotation): {v_descend_min} and take off: {v_ascend_min}')
    print(f'apollo ascent module with full fuel has {lander_dv} km/s dV, sufficient for {n_landings} landings')