import numpy as np
import spiceypy as spice

planets_phys = {'Earth':{
                    'gm':3.98600435136e5,
                    'r':6378.1363,
                    'rmean':6371.0084,
                    'm':5.9724e24,
                    'rho':5.5136e3,
                    'sday':86164,
                    'syr':1.0,
                    'vesc':11.19,
                    'g':9.80,
                    'j2':1.08262545e-3,
                    'p':101325.0,
                    'alb':0.306,
                    'Tsurf':288.0,
                    'Tbb':254.0,
                    'f':0.0033528
                    },
                'Mars':{
                    'gm':4.2828375214e4,
                    'r':3396.2,
                    'rmean':3389.5,
                    'm':6.41712e23,
                    'rho':3.9341e3,
                    'sday':1.02595676*86400,
                    'syr':1.8808476,
                    'vesc':5.03,
                    'g':3.72076,
                    'j2':1.96045e-3,
                    'p':636.0,
                    'alb':0.25,
                    'Tsurf':210.0,
                    'Tbb':209.8,
                    'f':0.00589
                    }
                }

planets_orb = {'Earth':{
                    'a':149598023.0,
                    'e':0.0167086,
                    'i':0.0,
                    'om':-11.26064,
                    'w':114.20783,
                    'ma':358.617
                    },
                'Mars':{
                    'a':227939200.0,
                    'e':0.0934,
                    'i':1.85,
                    'om':49.558,
                    'w':86.502,
                    'ma':19.412
                    }
                }

def load():
    META = 'spacepy/data/metakr.tm'
    spice.furnsh(META)

def unload():
    META = 'spacepy/data/metakr.tm'
    spice.unload(META)