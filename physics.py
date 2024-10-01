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
## File name: physics.py                                                     ##
## File content: functions from special relativity, Doppler effect           ##
##               and Planck's law of blackbody radiation (and asymptotics)   ##
###############################################################################

# -------------------
# importing libraries
# -------------------
import numpy as np
from numpy import cos, tan, sqrt, arccos, exp, log10
from numpy.linalg import norm

# all these are in SI unit
nm = 1E-9
h_Planck = 6.626E-34
c = 299792458 # duplicate, but makes some formula more familiar-looking
k_Boltzman = 1.38E-23

# l_cent is roughly the center wavelength of V filter
l_cent = 551 # nm

# threshold to use large T asymptotics of Planck's law, to avoid log(negative)
T_large = 2*h_Planck*c/(2*l_cent*nm*k_Boltzman)

# --------------------
# Functions: physics, SR, Doppler effect, Jonhson V filtered blackbody radiation "alpha" value
# --------------------

# Lorentz factor "gamma" in terms of beta ( = v/c)
def gamma_factor(b1,b2,b3): 
    return 1/sqrt(1-norm([b1,b2,b3])**2)

# Lorentz factor from beta
def gamma_frm_beta(beta):
    return 1/sqrt(1-beta**2)

# aberration of angle from Lorentz transform
def Lorentz_aber(beta, theta):
    return arccos((cos(theta)+beta)/(1+beta*cos(theta)))

# relativistic redshift factor, theta: see above in headlight_eff
# usually denoted by 1+z
def redsh(beta, theta):
    return gamma_frm_beta(beta)*(1+beta*cos(theta))

# effective blackbody radiation temperature
def bb_temp_adjust(T_K_emit, beta, theta):
    return T_K_emit*redsh(beta,theta)

# Planck's law of blackbody radiation and its asymptotics
def Planck_law(l_nm, T_K): # unit free
    return (2*h_Planck*c**2/(l_nm*nm)**5)*(1/(exp(h_Planck*c/(k_Boltzman*l_nm*nm*T_K))-1))

def Planck_law_small_asymp(l_nm, T_K):
    return (2*h_Planck*c**2)*exp(-h_Planck*c/(k_Boltzman*l_nm*nm*T_K))

def Planck_law_large_asymp(l_nm, T_K):
    return (2*c/(l_nm*nm)**4)*k_Boltzman*T_K-(h_Planck*c**2)/(l_nm*nm)**5

# asymptotics of the alpha function
def log10_alpha_small(T_K):
    return log10(Planck_law_small_asymp(l_cent, T_K))+37

def log10_alpha_large(T_K):
    return log10(Planck_law_large_asymp(l_cent, T_K))-13.4




