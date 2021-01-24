# external imports
import numpy as np

# internal imports
from spacepy.objects import create_LEO
from spacepy.data.bodydata import planet_data
from spacepy import simulate as sim

def main():
    sc = create_LEO(h_p=418.0, h_a=419.0, i=51.6467, w=207.309, lan=43.0535)
    sim.rk8(sc, 10*24*3600)
    sim.animate_twobody(sc)

if __name__=='__main__':
    main()