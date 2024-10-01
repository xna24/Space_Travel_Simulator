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
## File name: geometry.py                                                    ##
## File content: functions for computation in geometry, coord. transform     ##
###############################################################################

# -------------------
# importing libraries
# -------------------
import numpy as np

from numpy import sin, cos, sqrt, log, dot, arccos
from numpy.linalg import norm
arctan2 = lambda x,y: np.arctan2(y,x)
arcsinh = lambda x: log(x+sqrt(x**2+1))

# ---------
# Functions
# ---------

# convert btw spherical and Cartesian coordinate
def long_lat_to_UVW(l,b,D): 
    return [D*cos(l)*cos(b), D*sin(l)*cos(b), D*sin(b)]

def UVW_to_long_lat(U,V,W):
    if U==0 and V==0:
        l = 0
    else:
        l = arctan2(U,V)
    b = arctan2(sqrt(U**2+V**2),W)
    return l,b

def UVW_to_theta_phi(U,V,W):
    return arctan2(U,sqrt(V**2+W**2)), arctan2(V,W)

def theta_phi_to_unit_vec(theta,phi):
    return cos(theta), cos(phi)*sin(theta), sin(phi)*sin(theta)

# angle between two vectors in 3-space
def angle_btw(X, Y):
    return arccos(dot(X,Y)/(norm(X)*norm(Y)))

# l = longitude, b = latitude in radian
# x=[U,V,W] vector to be rotated
# return x', vector after first rotating along W-axis by -l and then along V-axis by b
# note: (cos(l)cos(b), sin(l)cos(b), sin(b)) rotates to (1,0,0)
def rotator(l,b,x):
    [U,V,W] = x
    return [cos(b)*cos(l)*U+cos(b)*sin(l)*V+sin(b)*W, -sin(l)*U+cos(l)*V, -cos(l)*sin(b)*U-sin(b)*sin(l)*V+cos(b)*W]
 