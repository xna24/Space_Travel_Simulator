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
## File name: main.py                                                        ##
## File content: producing sky maps                                          ##
###############################################################################

# --------------
# load libraries
# --------------

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import os
import time

from numpy import sin, cos, tan, sqrt, pi, log, dot, sinh, tanh, arccos, log10, exp
from numpy.linalg import norm
csc = lambda x: 1/sin(x)
sec = lambda x: 1/cos(x)
arctan2 = lambda x,y: np.arctan2(y,x)
arcsinh = lambda x: log(x+sqrt(x**2+1))

from travel import trip_info, trip_into_string_normal, trip_info_string_latex, gamma_now
from geometry import long_lat_to_UVW, rotator, UVW_to_theta_phi, theta_phi_to_unit_vec, \
                     UVW_to_long_lat
from physics import Lorentz_aber, bb_temp_adjust, Planck_law_small_asymp, \
                    Planck_law_large_asymp, T_large, l_cent, gamma_frm_beta, \
                    log10_alpha_small, log10_alpha_large
from astronomy import Plx_mas_to_Dist_ly, BV_to_temp, mag_to_size, \
                      residual_mV_term
from settings import settings_read

work_folder = os.getcwd()

# ----------------------
# constants
# ----------------------
px = 1/plt.rcParams['figure.dpi']

RA_hr_in_deg, RA_min_in_deg, RA_sec_in_deg = 360/24, 360/1440, 360/86400
DE_min_in_deg, DE_sec_in_deg = 1/60, 1/3600

deg_in_rad, hr_in_rad, minute_in_rad = pi/180, pi/12, pi/(12*60)

# ------------------
# Functions: general
# ------------------

# if a string represent a float number, return True, otherwise False
def is_num(s: str) -> bool: 
    try:
        float(s)
        return True
    except ValueError:
        return False

strip_it = lambda s: s.strip()

def str_num_add(lst):
    r = ''
    for j in range(len(lst)):
        x = lst[j]
        if x < 0 : 
            r += '-'+'{:.3f}'.format(-x)
        else:
            r += '+'+'{:.3f}'.format(x)
    return r[1:]

# ----------------------------
# Load data: B-V color to sRGB
# ----------------------------

temp_to_RGB_file = os.path.join(work_folder, 'temp_to_RGB.csv')

if not os.path.isfile(temp_to_RGB_file):
    sys.exit('\'temp_to_RGB.csv\' missing')
else:
    temp_to_RGB_DF = pd.read_csv(temp_to_RGB_file, sep = '\t')
    temp_std_list = temp_to_RGB_DF['temp'].to_numpy()
    R_std_list = temp_to_RGB_DF['R'].to_numpy()
    G_std_list = temp_to_RGB_DF['G'].to_numpy()
    B_std_list = temp_to_RGB_DF['B'].to_numpy()

# ----------------------------------------------------#
# Load data: Johnson V filter                         #
# ----------------------------------------------------#
# alpha ~ \int_0^\infty g(l)*B_T(l)dl                 #
# where B_T(l) is blackbody radiation at temperature T#
# normalized so alpha(5800K)=1                        #
# g(l) is the profile of Johnson V filter             #
# normalized so int_0^\infty g(l)dl=1                 #
# for T very large/small, g can be replaced by a Dirac#
# delta at roughly l_cent, so we use asymptotics of   #
# Planck's law to estimate alpha                      #
# ----------------------------------------------------#

# load data file for the interpolation values of alpha function
log10_alpha_file = os.path.join(work_folder, 'log10_alpha_data.csv')
if not os.path.isfile(log10_alpha_file):
    sys.exit('\'log10_alpha_data.csv\' missing')
else:
    DF_log10_alpha = pd.read_csv(log10_alpha_file, sep='|', index_col=None)

for_alpha_Tlst = DF_log10_alpha['T_K'].tolist()
for_alpha_log10alphalst = DF_log10_alpha['log10_alpha'].tolist()

def log10_alpha(T_K):
    if T_K < 100: # you're going to get divide by zero error for T_K too low, not addressing for the moment
        return log10_alpha_small(T_K)
    elif T_K > T_large: # note there is a slight drop in valid at T_large, not going to cause anything super weird
        return log10_alpha_large(T_K)
    else:
        return np.interp(T_K, for_alpha_Tlst, for_alpha_log10alphalst)

# Trip information and string to print on sky map

def trip_info_string(des_name, L_ly, T, t, tau, x_t, beta_t, trip_ratio, use_latex):
    if use_latex:
        return trip_into_string_normal(des_name, L_ly, T, t, tau, x_t, beta_t, trip_ratio)
    else:
        return trip_info_string_latex(des_name, L_ly, T, t, tau, x_t, beta_t, trip_ratio)

# ------------------------------------------------------------------------------#
# Master function                                                               #
# ------------------------------------------------------------------------------#
# Requires global variable: tot_star_computed, tot_star_to_compute
# X_data, Y_data, rad_lst, color_lst should be np.array of size DF_star.shape[0], 
# to be written in and passed to matplotlib.axes.Axes.plot. They are respectively: 
# horizontal and vertical coordinates, size of marker and color of marker.
# 
# DF_stars is a DataFrame with key columns in "col_names" list, containing :
# HIP, TYC, longitude, latitude, parallax, parallax validity (1=valid or 0=invalid),
# apparent magnitude (at Earth), B-V color
# All values in this DataFrame is assumed to be strings, HIP number contain no white spaces
# columns of long, lat, plx, plx_validity, mag, B-V must all have value (i.e. float(..) can be applied)
# Note we only travel to stars in Hipparcos catalog with more accurate parallaxes
# compute_list (for MasterF_debug only): list of row indices to compute in DF_stars
# Return a tuple of trip current status plus a count of visible stars
# debug_flag: if True, will print info for each star computed in "compute_list"
#             do NOT set to True if you're using it on all the stars, or else it prints too much
# -----------------------------------------------------------------------

# Update:
# MasterF rewritten: it now accesses more global variables
#                    in particular, l_star_LIST, b_star_LIST, plx_LIST
#                                   vmag_LIST, bv_LIST, plx_val_LIST, all float
#                    another change: no need for col_names!
#                    all these lists need to have same length LIST_length, also part of inputs
# I've temporarily taken away debug things

def MasterF(ldes, bdes, plx_des, trip_ratio, LIST_length, \
            azm_range_deg, alt_range_deg, vs_cutoff):
    
    # only do so for global variables we alter
    global tot_star_computed, computation_time_spent  
    global l_star_LIST, b_star_LIST, plx_LIST, vmag_LIST, plx_val_LIST
    global D_ly_LIST, T_emit_LIST
    global X_data, Y_data, rad_lst, color_lst

    
    azm_min_deg, azm_max_deg = azm_range_deg
    alt_min_deg, alt_max_deg = alt_range_deg

    L_ly, T, t, tau, x_t, beta_t = trip_info(plx_des, trip_ratio)
    
    count = 0

    for j in range(LIST_length):
        
        vmag = vmag_LIST[j]
        
        D_ly = D_ly_LIST[j]
        
        [Urot, Vrot, Wrot] = UVWrot_LIST[j]
        
        if star_catalog == 'HIP_TYC' or star_catalog == 'custom':
            if plx_val_LIST[j] == 1:
                x1, x2, x3 = Urot-x_t, Vrot, Wrot
                D_new_sq = x1**2+x2**2+x3**2 # sqrt costs time
            else:
                x1, x2, x3 = Urot, Vrot, Wrot # in this case D_new_sq won't even be used
        elif star_catalog == 'HIP':
            if plx_LIST[j] > 0:
                x1, x2, x3 = Urot-x_t, Vrot, Wrot
                D_new_sq = x1**2+x2**2+x3**2 # sqrt costs time
            else:
                x1, x2, x3 = Urot, Vrot, Wrot # in this case D_new_sq won't even be used
        
        theta, phi = UVW_to_theta_phi(x1, x2, x3) # spherical coordinate with velocity due north
        theta_new = Lorentz_aber(beta_t, theta) # Lorentz transform corrected, no longer instantaneous false "visual" spatial transform
        e1, e2, e3 = theta_phi_to_unit_vec(theta_new, phi)
        ls, bs = UVW_to_long_lat(e1, e2, e3)
        
        tot_star_computed += 1
        if tot_star_computed % 2500 == 0 and \
           not (tot_star_computed % (2500*50) == 0):
            print('\u2588', end='', flush=True) # \u2588 is a block in Unicode
        elif tot_star_computed % (2500*50) == 0:
            computation_time_spent = time.time() - computation_time_start
            print('\u2588| ' + \
            '{:.2f}'.format(100*tot_star_computed/tot_star_to_compute)+'%, '\
            +'{:.2f}'.format(computation_time_spent)+'s')

        if ls > azm_max_deg*deg_in_rad or ls < azm_min_deg*deg_in_rad: 
            if debug_flag: print('Azimuth out of view. (|Azimuth/deg|>'+'{:.3f}'.format(azm_max_deg)+')')
            continue
        if bs > alt_max_deg*deg_in_rad or bs < alt_min_deg*deg_in_rad: 
            if debug_flag: print('Altitude out of view. (|Altitude/deg|>'+'{:.3f}'.format(alt_max_deg)+')')
            continue
        
        if star_catalog == 'HIP_TYC':
            if plx_val_LIST[j] == 1:
                mag_diff_dist = 5*(0.5*log10(D_new_sq)-log10(D_ly))
            else:
                mag_diff_dist = 0
        elif star_catalog == 'HIP':
            if plx > 0:
                mag_diff_dist = 5*(0.5*log10(D_new_sq)-log10(D_ly))
            else:
                mag_diff_dist = 0

        # at this point the star is for sure in the field of vision
        X_data[count] = ls/deg_in_rad
        Y_data[count] = bs/deg_in_rad

        T_emit = T_emit_LIST[j]
        Ts = bb_temp_adjust(T_emit, beta_t, theta)
        alpha_diff = 2.5*log10_alpha(T_emit) - 2.5*log10_alpha(Ts)

        vs = vmag + alpha_diff + mag_diff_dist + residual_mV_term(beta_t, theta)

        if vs >= vs_cutoff: continue
        
        size = mag_to_size(vs,star_sz_mult,star_max_size,star_sz_pw)
        rad_lst[count] = size

        R_val, G_val, B_val = np.interp(Ts, temp_std_list, R_std_list),\
                              np.interp(Ts, temp_std_list, G_std_list),\
                              np.interp(Ts, temp_std_list, B_std_list)
        color_lst[count] = (R_val, G_val, B_val)

        count += 1
    
    return L_ly, T, t, tau, x_t, beta_t, count

# Rewritten from the update function fed to matplotlib.animation.FuncAnimation
# Need these global variables defined before calling:
# ldes, bdes, plx_des: info about destination
# trip_ratio: amount of trip completed in proper time
# fig, axs are matplotlib.figure.Figure, matplotlib.axes._axes.Axes objects
# DF, a pandas DataFrame containing star information
# col_names, list of strings, useful columns in DF (see MasterF)
# azm_range, alt_range, vs_cutoff (see MasterF)
# X_data, Y_data, rad_lst, color_lst: four arrays of zeros of length DF.shape[0]

def sky_map():
    global fig, axs
    fig.clear()

    axs = fig.add_subplot(111)
    axs.set_facecolor((0,0,0))
    axs.set_xlim(azm_range)
    axs.set_ylim(alt_range)

    if use_latex:
        axs.set_xlabel('azimuth (${}^\\circ$)')
        axs.set_ylabel('altitude (${}^\\circ$)')
    else:
        axs.set_xlabel('azimuth (deg)')
        axs.set_ylabel('altitude (deg)')

    axs.xaxis.set_label_coords(0.5,0.075)
    axs.yaxis.set_label_coords(0.96,0.5)

    axs.xaxis.label.set_color('w')
    axs.yaxis.label.set_color('w')
    axs.yaxis.tick_right()

    axs.xaxis.label.set_color('w')

    axs.tick_params(axis='x', direction='in', colors='w', pad=-20)
    axs.tick_params(axis='y', direction='in', colors='w', pad=-35)

    L_ly, T, t, tau, x_t, beta_t, count = MasterF(ldes, bdes, plx_des, trip_ratio, LIST_length,\
                                                  azm_range, alt_range, vs_cutoff)
    
    info_str = trip_info_string(des_name, L_ly, T, t, tau, x_t, beta_t, trip_ratio, use_latex)
    
    scat = axs.scatter(X_data[:count], Y_data[:count], marker = 'o', \
                       c = color_lst[:count], s = rad_lst[:count])
    txt = axs.text(0,1, info_str, color='w',\
                horizontalalignment='left', verticalalignment = 'top', transform=axs.transAxes)
     
    arr_img = plt.imread(os.path.join(work_folder, 'color_legend.png'))
    im = OffsetImage(arr_img)
    ab = AnnotationBbox(im, (0.975,1), xycoords = 'axes fraction', box_alignment=(1,1), pad = 0)
    axs.add_artist(ab)
    
    star_mag_legend_axes_fractions = [[0.04,0.15],[0.08,0.15],[0.12,0.15],[0.16,0.15],[0.2,0.15],\
                                      [0.04,0.1],[0.08,0.1],[0.12,0.1],[0.16,0.1],[0.2,0.1]]
    star_mag_legend_X_coord = [x[0]*(azimuth_max-azimuth_min)+azimuth_min for x in star_mag_legend_axes_fractions]
    star_mag_legend_Y_coord = [x[1]*(altitude_max-altitude_min)+altitude_min for x in star_mag_legend_axes_fractions]
    star_mag_legend_Y_minus_coord = [(x[1]-0.025)*(altitude_max-altitude_min)+altitude_min for x in star_mag_legend_axes_fractions]
    mag_lst = [x-4 for x in range(10)]
    star_mag_legend_txt_str = [str(x-4) for x in range(10)]
    star_mag_legend_size = [mag_to_size(x,star_sz_mult,star_max_size,star_sz_pw) for x in mag_lst]
    
    star_mag_legend_pt = axs.scatter(star_mag_legend_X_coord, star_mag_legend_Y_coord, c = 'w', s = star_mag_legend_size)
    
    star_mag_legend_txt = [axs.text(star_mag_legend_X_coord[j],\
                                    star_mag_legend_Y_minus_coord[j],\
                                    star_mag_legend_txt_str[j], c = 'w', fontsize=12,\
                                    horizontalalignment='center',\
                                    verticalalignment='center') for j in range(10)]
    star_mag_legend_txt += [axs.text(0.12*(azimuth_max-azimuth_min)+azimuth_min,\
                                     0.175*(altitude_max-altitude_min)+altitude_min, \
                                     'Visual Magnitude', c='w', fontsize = 15,\
                                    horizontalalignment='center',\
                                    verticalalignment='center')]
    pass
    

# --------------------------------------------#
# Runtime                                     #
# --------------------------------------------#    

# --------------------------------------------#
# main trip_ratio output_file                 #
# --------------------------------------------#

overriden = False
debug_flag = False
cmd_args = sys.argv[1:]
if len(cmd_args) >= 1:
    if not is_num(cmd_args[0]):
        sys.exit('Trip ratio not a number.')
    else:
        trip_ratio = float(cmd_args[0])
        if trip_ratio < 0.001 or trip_ratio > 0.999:
            sys.exit('Trip ratio too close to 0 or 1 or out of range.')
    if len(cmd_args) >= 2:
        sky_map_file = os.path.join(work_folder, cmd_args[1])
    else:
        sky_map_file = 'sky_map.png'
    overriden = True

if len(cmd_args) >= 2:
    if cmd_args[1] == '-debug': 
        debug_flag = True
        overriden = True
        print('Debug mode')
        if len(cmd_args) < 3: sys.exit('HIP of star to debug not given.')
        HIP_star = cmd_args[2]
        if not is_num(HIP_star): sys.exit('HIP of star to debug not a number.')
        print('Debugging star HIP '+HIP_star)
        HIP_star = int(HIP_star)
        

settings_file = os.path.join(work_folder, 'settings.ini')

# load parameters from "settings.ini"
use_latex, star_catalog, pic_width, pic_height, HIP_des, des_name, \
col_names, star_sz_mult, star_max_size, star_sz_pw, azimuth_span, altitude_span,\
vs_cutoff, animated, trip_ratio_lst, file_head, file_tail, frame_extn,\
ldes, bdes, plx_des = settings_read(settings_file)

if debug_flag:
    print('------------- Begin: Settings -------------------')
    print('Using LaTeX: '+str(use_latex))
    print('Star catalog used: '+star_catalog)
    print('Sky map dimension: '+str(pic_width)+'x'+str(pic_height))
    print('Destination name: '+des_name)
    if HIP_des != 0: print('HIP of destination: ' + str(HIP_des))
    print('DataFrame column names: '+str(col_names))
    print('Star size multiplier: '+str(star_sz_mult))
    print('Star size max.: '+str(star_max_size))
    print('Star size power-law power.:'+str(star_sz_pw))
    print('Azimuth span (deg): '+str(azimuth_span))
    print('Altitude span (deg): '+str(altitude_span))
    print('Vmag cutoff: '+str(vs_cutoff))
    print('Animated? '+str(animated))
    
    if animated:
        print('Trip ratio list: '+str(trip_ratio_lst))
        print('File name head/tail: '+file_head+'/'+file_tail)
        print('Frame file extension: '+frame_extn)
    
    if star_catalog == 'custom':
        print('RA, Dec(deg) of dest.: '+str(ldes/deg_in_rad)+', '+str(bdes/deg_in_rad))
        print('Destination parallax(mas):'+str(plx_des))
        
    print('-------------   End: Settings -------------------')
    
    
# --------------------------------------------#
# Load star catalog                           #
# --------------------------------------------#
# use_HIP_TYC = True:                         #
# Load Hipparcos and Tycho catalog (combined) #
# use_HIP_TYC = False:                        #
# Load Hipparcos catalog only                 #
# --------------------------------------------#

if star_catalog == 'HIP_TYC':
    database_path = os.path.join(work_folder, 'HIP_TYC.csv')
    if not os.path.isfile(database_path):
        sys.exit('HIP_TYC.csv not found')
    else:
        DF = pd.read_csv(database_path, index_col = False, sep='|', low_memory=False, dtype=str)
elif star_catalog == 'HIP':
    database_path = os.path.join(work_folder, 'HIP.csv')
    if not os.path.isfile(database_path):
        sys.exit('HIP.csv not found')
    else:
        DF = pd.read_csv(database_path, sep='|', low_memory=False, dtype=str)
elif star_catalog == 'custom':
    if 'StarCatalogFileName' in setting_keys:
        database_path = os.path.join(work_folder, settings['StarCatalogFileName'])
        if not os.path.isfile(database_path):
            sys.exit('Custom star catalog file missing')
        else:
            DF = pd.read_csv(database_path, sep='|', low_memory=False, dtype=str)
    else:
        sys.exit('StarCatalog option invalid in \'settings.ini\'')

LIST_length = DF.shape[0]

if star_catalog == 'HIP_TYC':
    [l_col, b_col, plx_col, plx_val_col, vmag_col, bv_col] = col_names
elif star_catalog == 'HIP' or star_catalog == 'custom':
    [l_col, b_col, plx_col, vmag_col, bv_col] = col_names

l_star_LIST = [float(x)*deg_in_rad for x in DF[l_col].tolist()]
b_star_LIST = [float(x)*deg_in_rad for x in DF[b_col].tolist()]
plx_LIST = [float(x) for x in DF[plx_col].tolist()]
vmag_LIST = [float(x) for x in DF[vmag_col].tolist()]
bv_LIST = [float(x) for x in DF[bv_col].tolist()]
T_emit_LIST = [BV_to_temp(bv) for bv in bv_LIST]

if star_catalog == 'HIP_TYC':
    plx_val_LIST = [int(x) for x in DF[plx_val_col].tolist()]
else:
    plx_val_LIST = [1]*LIST_length
    
# zero and negative parallax will be treated as being infinitely far
if star_catalog == 'HIP_TYC' or star_catalog == 'custom':
    D_ly_LIST = [Plx_mas_to_Dist_ly(plx_LIST[j]) if (plx_val_LIST[j] == 1 and plx_LIST[j]>0) else 1 for j in range(LIST_length)]
elif star_catalog == 'HIP':
    D_ly_LIST = [Plx_mas_to_Dist_ly(plx_LIST[j]) if (plx_val_LIST[j] == 1 and plx_LIST[j]>0) else 1 for j in range(LIST_length)]

print('Star catalog loaded.\nTot. star count: '+str(DF.shape[0]))

# destination info
if star_catalog == 'HIP_TYC':
    DF_des = DF[DF['HIP']==str(HIP_des)]
    if DF_des.shape[0] == 0:
        sys.exit('Destination HIP number not found in star catalog.')
    ldes, bdes, plx_des = float(DF_des.iloc[0]['RAdeg'])*deg_in_rad, \
                          float(DF_des.iloc[0]['DEdeg'])*deg_in_rad,\
                          float(DF_des.iloc[0]['Plx'])
elif star_catalog == 'HIP':
    DF_des = DF[DF['hip_number']==str(HIP_des)]
    if DF_des.shape[0] == 0:
        sys.exit('Destination HIP number not found in star catalog.')
    ldes, bdes, plx_des = float(DF_des.iloc[0]['lii'])*deg_in_rad, \
                          float(DF_des.iloc[0]['bii'])*deg_in_rad, \
                          float(DF_des.iloc[0]['parallax'])

if Plx_mas_to_Dist_ly(plx_des) > 1000:
    print('Distance to destination > 1000 ly, mid-travel scene not going to look pretty!')

if ldes > pi: ldes -= 2*pi

UVWstar_LIST = [long_lat_to_UVW(l_star_LIST[j], b_star_LIST[j], D_ly_LIST[j]) for j in range(LIST_length)]
UVWrot_LIST = [rotator(ldes, bdes, UVWstar) for UVWstar in UVWstar_LIST]

if not HIP_des == 0:
    des_name_show = des_name + ' HIP '+str(int(HIP_des))
else:
    des_name_show = des_name

print('Destination: '+des_name_show)
print('Distance: '+'{:.3f}'.format(Plx_mas_to_Dist_ly(plx_des))+' ly')
print('Right Ascension(deg): '+'{:.3f}'.format(ldes/deg_in_rad)) # will be corrected in current HIP catalog file
print('Declination(deg): '+'{:.3f}'.format(bdes/deg_in_rad))

# prepare lists to be written into
X_data = [0]*DF.shape[0]
Y_data = [0]*DF.shape[0]
rad_lst = [0]*DF.shape[0]
color_lst = [0]*DF.shape[0]

azimuth_max, azimuth_min = azimuth_span, -azimuth_span
altitude_max, altitude_min = altitude_span, -altitude_span

azm_range = (azimuth_min,azimuth_max)
alt_range = (altitude_min,altitude_max)

# if we receive valid command-line parameters, we will draw a single sky map
if overriden: animated = False


# -------------------------------------#
# Plot sky map frame(s)                #
# -------------------------------------#

tot_star_computed = 0
computation_time_start = time.time()
computation_time_spent = 0

if animated and not debug_flag:

    tot_frame = len(trip_ratio_lst)
    tot_star_to_compute = tot_frame * DF.shape[0]
    print('Tot. frame to draw: '+str(tot_frame))
    
    plt.rcParams.update({'font.size': 20})
    
    for frame in range(len(trip_ratio_lst)):
        fig = plt.figure(figsize = (pic_width*px,pic_height*px), facecolor = (0,0,0))
        trip_ratio = trip_ratio_lst[frame]
        plt.rcParams.update({'font.size': 20})
        fig = plt.figure(figsize = (pic_width*px,pic_height*px), facecolor = (0,0,0))
        sky_map()
        digit_count = int(np.ceil(log10(len(trip_ratio_lst))))+1
        plt.savefig(os.path.join(work_folder, file_head+'_'+str(frame).zfill(digit_count)+'_'+file_tail+'.'+frame_extn), bbox_inches='tight')
        plt.close('all')
        computation_time_spent = time.time() - computation_time_start
        # print('\nTime spent: '+'{:.3f}'.format(computation_time_spent)+'sec')
elif not debug_flag:
    fig = plt.figure(figsize = (pic_width*px,pic_height*px), facecolor = (0,0,0))
    tot_star_to_compute = DF.shape[0]
    if not overriden: # if no command-line input detected
        trip_ratio = 0.2
        sky_map_file = 'sky_map_example.png'
    plt.rcParams.update({'font.size': 20})
    fig = plt.figure(figsize = (pic_width*px,pic_height*px), facecolor = (0,0,0))
    sky_map()
    plt.savefig(os.path.join(work_folder, sky_map_file), bbox_inches='tight')
    computation_time_spent = time.time() - computation_time_start
    # print('\nTime spent: '+'{:.3f}'.format(computation_time_spent)+'sec')
else: # we are debugging
    if star_catalog == 'HIP_TYC' or star_catalog == 'custom':
        index_debug_lst = DF.index[DF['HIP']==str(HIP_star)].tolist()
    else:
        index_debug_lst = DF.index[DF['hip_number']==str(HIP_star)].tolist()
    if len(index_debug_lst) == 0: sys.exit('Star HIP '+str(HIP_star)+' not found')
    
    L_ly, T, t, tau, x_t, beta_t, count = MasterF_debug(ldes, bdes, plx_des, trip_ratio, \
                                                  DF, index_debug_lst, col_names, \
                                                  azm_range, alt_range, vs_cutoff)
                                                  
