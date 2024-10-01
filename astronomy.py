###############################################################################
## A Relativisitic Space Travel Simulator                                    ##
## Author: Xuesen Na                                                         ##
## Description: Produce sky maps seen by a spaceship during space travel.    ##
## The spaceship travels with constant proper acceleration followed by       ##
## deceleration of 9.8 m/s^2, reaches relativistic speed for any destination ##
## of distance on the order of light-year. The sky map incorporates all      ##
## the stars with valid distance data in the Hipparcos and Tycho catalog     ##  
## and uses equirectangular projection. The sky map seen during travel       ##
## exhibits these interesting visual effects: time dilation, relativistic    ##
## aberration and headlight effect, relativistic Doppler effect              ##
## (blue-shift and red-shift).                                               ##
## File name: astronomy.py                                                   ##
## File content: functions for computation in astronomy                      ##
###############################################################################

# -------------------
# importing libraries
# -------------------

from numpy import log10, cos
from physics import gamma_frm_beta

parsec = 3.26156 # ly
abs_mag_sun = 4.83 # absolute magnitude of the sun

# --------------------
# Functions: astronomy
# --------------------

# converting between apparent and absolute magnitudes, using distance (ly)
def abs_mag(Vmag, D_ly):
    return Vmag - 5*(log10(D_ly/parsec)-1)
def app_mag(Vabs, D_ly):
    return Vabs + 5*(log10(D_ly/parsec)-1)

# parallax (mas) to distance (ly)
def Plx_mas_to_Dist_ly(Plx):
    return (1000/Plx)*parsec

# determine the star size, set a maximum size
def mag_to_size(vmag,star_size_mult,star_size_max,power=2.5):
    return min(star_size_mult*10**(-float(vmag)/power),star_size_max)

# Johnson B-V color to temperature in Kelvin
def BV_to_temp(bv: float) -> float:
    return 4600 * (1 / (0.92 * bv + 1.7) + 1 / (0.92 * bv + 0.62) )
    
# one of the terms of SR correction to effective V-band magnitude
def residual_mV_term(beta, theta):
    return 5*log10(gamma_frm_beta(beta)*(1+beta*cos(theta)))