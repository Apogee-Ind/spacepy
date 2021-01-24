# external imports
import spiceypy as spice
import numpy as np

# internal imports
from spacepy.objects import Planet, Sol, System, SpaceCraft, SmallBody, Moon
import spacepy.data.bodydata as planetdata
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

    planetdata.load()
    epoch = '2021 JAN 01 00:00:00'
    epoch_end = '2028 JAN 01 00:00:00'
    step = 2*86400

    sys = System(epoch, Sun, Earth, Venus, Mercury)
    #print(sys.contents)
    sys._gen_ssb_vectors(epoch_end, step)
    sim.plot_system(sys, do_markers=True)
    

if __name__ == "__main__":
    main()