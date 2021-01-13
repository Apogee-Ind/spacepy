# external imports
import numpy as np

# internal imports
from spacepy.objects import SpaceCraft, Planet
from spacepy.data.planetdata import planets_orb, planets_phys
from spacepy import simulate as sim

def main():
    Earth = Planet()
    sat = SpaceCraft('Sat01', 100.0, shape='cuboid', dims=(2.0, 2.0, 2.0))
    
    sat.place_in_orbit(Earth, 418.0, 419.0, 51.6467, 207.309, 43.0535, 0.0)
    sat.add_thruster(90e-3, 3100.0)
    sat.add_fuel(100.0)

    thrust_spec = np.array([[0.0, 1]])

    sim.rk8(sat, 30*24*3600, thrust_spec=thrust_spec, j2=True, drag=False, rtol=1e-6, atol=1e-6, max_step=120.0)
    print(f'Simulation start time: {sat.t[0]} s\nSimulation end time: {sat.t[-1]} s')
    print(f'Fuel remaining: {sat.m_fuel} kg')
    print(f'Final spacecraft mass: {sat.m} kg')
    print(f'Initial velocity (magnitude): {np.linalg.norm(sat.vvec[1])} km/s')
    print(f'Final velocity (magnitute): {np.linalg.norm(sat.vvec[-1])} km/s')
    #sim.animate_twobody(sat)
    sim.plot_twobody(sat)
    

if __name__ == '__main__':
    main()

