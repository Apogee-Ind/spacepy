import numpy as np
import spiceypy as spice

# 8 recognized planets
planet_data = {
    'Mercury':{
        'id':199,
        #parameters from JPL solar system dynamics
        'gm':2.2032e4,
        'r':2440.53,
        'rmean':2439.4,
        'm':3.30114e23,
        'rho':5.4291e3,
        'sday':58.6462*86400,
        'syr':0.2408467,
        'mag': -0.60,
        'alb_g':0.106,
        'g':3.7,
        'vesc':4.25,
        #parameters from nasa planet fact sheet
        'j2':50.3e-6,
        'f':0.0,
        'p':0.0,
        'Tsurf':439.6,
        'Tbb':439.6,
        'alb_b':0.068
    },
    'Venus':{
        'id':299,
        #parameters from JPL solar system dynamics
        'gm':3.2486e5,
        'r':6051.8,
        'rmean':6051.8,
        'rho':5.243e3,
        'sday': -243.018*86400,
        'syr':0.61519726,
        'mag': -4.47,
        'alb_g':0.65,
        'g':8.87,
        'vesc':10.36,
        #parameters from nasa planet fact sheet
        'j2':4.458e-6,
        'f':0.0,
        'p':92*101325.0,
        'Tsurf':737,
        'Tbb':226.6,
        'alb_b':0.77
    },
    'Earth':{
        'id':399,
        #parameters from JPL solar system dynamics
        'gm':3.98600435136e5,
        'r':6378.1366,
        'rmean':6371.0084,
        'm':5.9724e24,
        'rho':5.5136e3,
        'sday':86164,
        'syr':1.0,
        'mag': -3.86,
        'alb_g':0.367,
        'g':9.80,
        'vesc':11.19,
        #parameters from nasa planet fact sheet
        'j2':1.08262545e-3,
        'f':0.0033528,
        'p':101325.0,
        'Tsurf':288.0,
        'Tbb':254.0,
        'alb_b':0.306
    },
    'Mars':{
        'id':499,
        #parameters from JPL solar system dynamics
        'gm':4.2828375214e4,
        'r':3396.19,
        'rmean':3389.5,
        'm':6.41712e23,
        'rho':3.9341e3,
        'sday':1.02595676*86400,
        'syr':1.8808476,
        'mag': -1.52,
        'alb_g':0.15,
        'g':3.72076,
        'vesc':5.03,
        #parameters from nasa fact sheet
        'j2':1.96045e-3,
        'f':0.00589,
        'p':636.0,
        'Tsurf':210.0,
        'Tbb':209.8,
        'alb_b':0.25
    },
    'Jupiter':{
        'id':599,
        #parameters from JPL solar system dynamics
        'gm':1.26686536e8,
        'r':71492,
        'rmean':69911,
        'm':1.89819e27,
        'rho':1.3262e3,
        'sday':0.44401*86400,
        'syr':11.862615,
        'mag': -9.40,
        'alb_g':0.52,
        'g':24.79,
        'vesc':60.2,
        'j2':14695.6e-6,
        'j4':-591.3e-6
    },
    'Saturn':{
        'id':699,
        #parameters from JPL solar system dynamics
        'gm':3.7931208e7,
        'r':60268,
        'rmean':58232,
        'm':5.683174e26,
        'rho':0.6871e3,
        'sday':0.44401*86400,
        'syr':29.447498,
        'mag':-8.88,
        'alb_g':0.47,
        'g':10.44,
        'vesc':36.09,
        'j2':16290.7e-6,
        'j4':-935.8e-6
    },
    'Uranus':{
        'id':799,
        #parameters from JPL solar system dynamics
        'gm':5.793951e6,
        'r':25559,
        'rmean':25362,
        'm':8.68127e25,
        'rho':1.27e3,
        'sday':-0.71833,
        'syr':84.016846,
        'mag':-7.19,
        'alb_g':0.51,
        'g':8.87,
        'vesc':21.38,
        'j2':3510.7e-6,
        'j4':-34.2e-6
    },
    'Neptune':{
        'id':899,
        #parameters from JPL solar system dynamics
        'gm':6835100e6,
        'r':24764,
        'rmean':24622,
        'm':1.024126e26,
        'rho':1.638e3,
        'sday':0.67125*86400,
        'syr':164.79132,
        'mag':-6.87,
        'alb_g':0.41,
        'g':11.15,
        'vesc':23.56,
        'j2':3408.4e-6,
        'j4':-33.4e-6
    }
}

# NAIF ID codes for planetary system barycenters
planet_systems = {
    'Mercury':1,
    'Venus':2,
    'Earth':3,
    'Mars':4,
    'Jupiter':5,
    'Saturn':6,
    'Uranus':7,
    'Neptune':8,
    'Pluto':9
}

# natural satellites of planets
moon_data = {
    'Moon':{
        'id':301,
        'system':'Earth',
        'gm':4902.801,
        'r':1738.1,
        'rmean':1737.4,
        'm':7.346e22,
        'rho':3344,
        'sday':27.321661*86400,
        'vesc':2.38,
        'g':1.62,
        'j2':202.7e-6,
        'p':1e-7,
        'alb_b':0.11,
        'Tsurf':270.4,
        'Tbb':270.4,
        'f':0.0012
    }
}

# dwarf planets, asteroids, and comets
smallbody_data = {
    'Pluto':{
        'id':999,
        'gm':8.7e2,
        'r':1188.3,
        'rmean':1188.3,
        'm':1.303e22,
        'rho':1.89e3,
        'sday':-6.3872*86400,
        'syr':247.92065,
        'mag':-1.0,
        'alb_g':0.3,
        'g':0.62,
        'vesc':1.21
    },
    'Ceres':{
        'id':2000001,
        'gm':62.809393,
        'r':939.4*0.5,
        'extent':[964.4, 964.2, 891.8],
        'rho':2.162e3,
        'sday':9.07417*3600,
        'alb_g':0.09
    },
    'Pallas':{
        'id':2000002,
        'gm':13.923011,
        'r':545*0.5,
        'extent':[582, 556, 500],
        'sday':7.8132*3600,
        'alb_g':0.101
    },
    'Vesta':{
        'id':2000004,
        'gm':17.288009,
        'r':525.4*0.5,
        'extent':[572.6, 557.2, 446.4],
        'rho':3.456e3,
        'sday':5.34212766*3600,
        'alb_g':0.4228
    },
    'Psyche':{
        'id':2000016,
        'gm':1.530048,
        'r':226*0.5,
        'extent':[279, 232, 189],
        'rho':4.5e3,
        'sday':4.196*3600,
        'alb_g':0.1203
    },
    'Lutetia':{
        'id':2000021,
        'gm':0.113442,
        'r':95.76*0.5,
        'sday':8.1655*3600,
        'alb_g':0.2212
    },
    'Kleopatra':{
        'id':2000216,
        'gm':0.309813,
        'r':122*0.5,
        'extent':[276, 94, 78],
        'sday':5.385*3600,
        'alb_g':0.1164
    },
    'Eros':{
        'id':2000433,
        'gm':4.463e-4,
        'r':16.84*0.5,
        'extent':[34.4, 11.2, 11.2],
        'rho':2.67e3,
        'sday':5.27*3600,
        'alb_g':0.25
    }

}


def load():
    META = 'spacepy/data/metakr.tm'
    spice.furnsh(META)

def unload():
    META = 'spacepy/data/metakr.tm'
    spice.unload(META)