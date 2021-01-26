# external imports
import spiceypy as spice
import numpy as np

# internal imports
from spacepy.objects import Planet, Sol, System, SpaceCraft, SmallBody, Moon
import spacepy.data.bodydata as bodydata
import spacepy.simulate as sim

def main():
    Earth = Planet(barycenter=True)
    Sun = Sol()
    Mercury = Planet('Mercury')
    Mars = Planet('Mars', barycenter=True)
    Venus = Planet('Venus')
    Jupiter = Planet('Jupiter', barycenter=True)
    Saturn = Planet('Saturn', barycenter=True)
    Ceres = SmallBody()
    Uranus = Planet('Uranus', barycenter=True)
    Neptune = Planet('Neptune', barycenter=True)
    Vesta = SmallBody('Vesta')
    Pallas = SmallBody('Pallas')
    Eros = SmallBody('Eros')
    #Psyche = SmallBody('Psyche')

    bodydata.load()
    epoch = '2045 DEC 31 00:00:00'
    epoch_end = '2047 DEC 31 00:00:00'
    step = 2*86400

    sys = System(epoch, Sun, Earth, Mars, Ceres)
    #print(sys.contents)
    sys._gen_ssb_vectors(epoch_end, step)
    #print(np.shape(sys.contents['major_planet'][3].state[0]['state']))
    #print(np.shape(sys.contents['star'][10].state[0]['state']))
    sim.plot_system(sys, 0, 'ECLIPJ2000', do_markers=True)
    #print(sys.nt)
    

if __name__ == "__main__":
    main()