# external imports
import spiceypy as spice
import numpy as np

# internal imports
from spacepy.objects import Planet, Sol, System, SpaceCraft
import spacepy.data.planetdata as planetdata

def main():
    Earth = Planet()
    Sun = Sol()
    Mars = Planet('Mars')
    Jupiter = Planet('Jupiter')


    planetdata.load()
    epoch = '2020 DEC 31 00:00:00'
    et_start = spice.str2et(epoch)
    epoch_end = '2021 JAN 20 00:00:00'
    et_end = spice.str2et(epoch_end)

    sys = System(et_start, Sun, Earth, Mars, Jupiter, sat)
    print(sys.contents)
    

if __name__ == "__main__":
    main()