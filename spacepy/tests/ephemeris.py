# external imports
import spiceypy as spice
import numpy as np

# internal imports
from spacepy.objects import Planet
import spacepy.data.planetdata as planetdata

def main():
    Earth = Planet()
    planetdata.load()
    epoch = '2020 DEC 31 00:00:00'
    et_start = spice.str2et(epoch)
    epoch_end = '2021 JAN 20 00:00:00'
    et_end = spice.str2et(epoch_end)

    et = [et_start, et_end]
    state = spice.spkssb(399, et_start, 'J2000')
    print(state)

if __name__ == "__main__":
    main()