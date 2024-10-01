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
## File name: travel.py                                                      ##
## File content: functions for computations of the space travel              ##
###############################################################################

# -------------------
# importing libraries
# -------------------

from numpy import sqrt, log, sinh
from astronomy import Plx_mas_to_Dist_ly

arcsinh = lambda x: log(x+sqrt(x**2+1))

# ----------------------------------------------
# constants from geometry, physics and astronomy
# ----------------------------------------------
speed_of_light = 299792458 # m/s
earth_grav = 9.8 # m/s**2
year_in_sec = 3600*24*365.25
ly_in_m = speed_of_light * year_in_sec

T_gc_SI = speed_of_light/earth_grav # c/g in sec
T_gc_yr = (speed_of_light/earth_grav)/year_in_sec # c/g in year
L_gc_ly = (speed_of_light**2/earth_grav)/ly_in_m # c**2/g in ly

# ------------------------------------
# Functions: relativistic space travel
# ------------------------------------

# maximal gamma factor reached during travel, from L the total distance traveled (in ly)
def gamma_max_frm_dist_ly(dist_ly):
    return dist_ly/(2*L_gc_ly)+1
# computes T from L
def travel_half_time_yr(dist_ly):
    return T_gc_yr*sqrt(gamma_max_frm_dist_ly(dist_ly)**2-1)
# total proper time span (in years) from half-time T (also in years)
def tot_prop_time_frm_T(T_yr):
    return 2*(T_gc_yr)*arcsinh(T_yr/T_gc_yr)

# computes coord time t from proper time tau and half-time T (both in yr)
def coord_time_frm_prop_time(tau, T):
    if tau<0 or tau>2*tot_prop_time_frm_T(T):
        return -1
    elif tau<= (1/2)*tot_prop_time_frm_T(T):
        return T_gc_yr*sinh(tau/T_gc_yr)
    else:
        return 2*T-T_gc_yr*sinh(2*arcsinh(T/T_gc_yr)-tau/T_gc_yr)

# computes coord dist covered x(t) (in ly) in terms of coord time (in yr) and half-time T (in yr)
def dist_cov_ly_frm_t_yr(t, T):
    if t<0 or t>2*T:
        return -1
    elif t<=T:
        return L_gc_ly*(sqrt((t/T_gc_yr)**2+1)-1)
    else:
        return L_gc_ly*(2*sqrt((T/T_gc_yr)**2+1)-sqrt(((2*T-t)/T_gc_yr)**2+1)-1)

# computes current speed
def speed_now(t,T):
    if t<=T:
        return (t/T_gc_yr)/sqrt((t/T_gc_yr)**2+1)
    else:
        return ((2*T-t)/T_gc_yr)/sqrt(((2*T-t)/T_gc_yr)**2+1)

# computes current Lorentz factor
def gamma_now(t,T):
    if t<=T:
        return sqrt((t/T_gc_yr)**2+1)
    else:
        return sqrt(((2*T-t)/T_gc_yr)**2+1)

# computes coord normalized (/c) velocity beta(t) in terms of t and T (both in yr)
def beta_frm_t_yr(t, T):
    if t<0 or t>2*T:
        return -1
    elif t<=T:
        return (t/T_gc_yr)/sqrt((t/T_gc_yr)**2+1)
    else:
        return ((2*T-t)/T_gc_yr)/sqrt(((2*T-t)/T_gc_yr)**2+1)

def trip_info(plx_des, trip_ratio):
    # distance to destination in ly
    L_ly = Plx_mas_to_Dist_ly(plx_des)

    T = travel_half_time_yr(L_ly)
    tau = trip_ratio * tot_prop_time_frm_T(T)
    t = coord_time_frm_prop_time(tau,T)
    x_t, beta_t = dist_cov_ly_frm_t_yr(t,T), beta_frm_t_yr(t,T)
    return L_ly, T, t, tau, x_t, beta_t
    
def trip_into_string_normal(des_name, L_ly, T, t, tau, x_t, beta_t, trip_ratio):

    speed_str = '{:.3f}'.format(speed_now(t,T))
    if speed_str == '1.000': speed_str = '>0.999'
    
    return 'Dest.: ' + des_name + '\n'\
               +'Dist.: '+'{:.3f}'.format(L_ly)+' ly\n'\
               +'Tot. trip dur. (coord. time): '+'{:.3f}'.format(2*T)+' yr\n'\
               +'Dist. covered: '+'{:.3f}'.format(x_t)+' ly = ' + '{:.3f}'.format((x_t/L_ly)*100) +  '% tot.\n'\
               +'Speed $v=$'+speed_str+' c\n'\
               +'Lorentz factor $\\gamma=$'+'{:.3f}'.format(gamma_now(t,T))+'\n'\
               +'Coord. time $t=$'+'{:.3f}'.format(t)+' yr = '+'{:.2f}'.format((t/(2*T))*100)+'% tot.\n'\
               +'Proper time $\\tau=$'+'{:.3f}'.format(tau)+' yr = '+'{:.2f}'.format(trip_ratio*100)+'% tot.\n'

def trip_info_string_latex(des_name, L_ly, T, t, tau, x_t, beta_t, trip_ratio):
    
    speed_str = '{:.3f}'.format(speed_now(t,T))
    if speed_str == '1.000': speed_str = '>0.999'
    
    return 'Dest.: ' + des_name + '\n'\
               +'Dist.: '+'{:.3f}'.format(L_ly)+' ly\n'\
               +'Tot. trip dur. (coord. time): '+'{:.3f}'.format(2*T)+' yr\n'\
               +'Dist. covered: '+'{:.3f}'.format(x_t)+' ly = ' + '{:.3f}'.format((x_t/L_ly)*100) +  '% tot.\n'\
               +'Speed v='+speed_str+' c\n'\
               +'Lorentz factor: '+'{:.3f}'.format(gamma_now(t,T))+'\n'\
               +'Coord. time t='+'{:.3f}'.format(t)+' yr = '+'{:.2f}'.format((t/(2*T))*100)+'% tot.\n'\
               +'Proper time tau='+'{:.3f}'.format(tau)+' yr = '+'{:.2f}'.format(trip_ratio*100)+'% tot.\n'