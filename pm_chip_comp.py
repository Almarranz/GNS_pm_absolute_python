#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 12:59:21 2025

@author: amartinez
"""

import numpy as np
from astropy import units as u
from astropy.time import Time
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import sys
from matplotlib import rcParams
from astroquery.gaia import Gaia
import astroalign as aa
import time
from astropy.coordinates import match_coordinates_sky
from compare_lists import compare_lists
import IPython
import copy
from skimage import data
from skimage import transform
import math
import glob
import skimage as ski
from astropy.table import Table, hstack
from scipy.stats import sigmaclip
from astropy.stats import sigma_clip
import os
import cluster_finder
from filters import filter_gaia_data
from filters import filter_hosek_data
from filters import filter_gns_data
from filters import filter_vvv_data
import Polywarp as pw
from alignator import alignator
from plot_machine import diff_hist
# %%plotting parametres
from matplotlib import rc
from matplotlib import rcParams
rcParams.update({'xtick.major.pad': '7.0'})
rcParams.update({'xtick.major.size': '7.5'})
rcParams.update({'xtick.major.width': '1.5'})
rcParams.update({'xtick.minor.pad': '7.0'})
rcParams.update({'xtick.minor.size': '3.5'})
rcParams.update({'xtick.minor.width': '1.0'})
rcParams.update({'ytick.major.pad': '7.0'})
rcParams.update({'ytick.major.size': '7.5'})
rcParams.update({'ytick.major.width': '1.5'})
rcParams.update({'ytick.minor.pad': '7.0'})
rcParams.update({'ytick.minor.size': '3.5'})
rcParams.update({'ytick.minor.width': '1.0'})
rcParams.update({'font.size': 20})
rcParams.update({'figure.figsize':(10,5)})
rcParams.update({
    "text.usetex": False,
    "font.family": "sans",
    "font.sans-serif": ["Palatino"]})
plt.rcParams["mathtext.fontset"] = 'dejavuserif'
rc('font',**{'family':'serif','serif':['Palatino']})
plt.rcParams.update({'figure.max_open_warning': 0})# a warniing for matplot lib pop up because so many plots, this turining it of
# Enable automatic plotting mode
# IPython.get_ipython().run_line_magic('matplotlib', 'auto')
IPython.get_ipython().run_line_magic('matplotlib', 'inline')



field_one, chip_one, field_two, chip_two,t1,t2,max_sig = np.loadtxt('/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/lists/fields_and_chips.txt', 
                                                       unpack=True)


pruebas2 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2absolute_python/pruebas/'
pruebas1 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/pruebas/'
gaia_list1 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/Gaia_lists/'
gaia_list2 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2absolute_python/Gaia_lists/'
tmp1 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/tmp/'
tmp2 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2absolute_python/tmp/'
Hos_cat = '/Users/amartinez/Desktop/PhD/Arches_and_Quintuplet_Hosek/'
gaia_hos = fits.open(Hos_cat + 'quintuplet_HST_gaia_matched_pm_EOM.fits')
pm_folder_off ='/Users/amartinez/Desktop/PhD/HAWK/pm_gns1_gns2_absolute_python/'
# field_one = field_one.astype(int)
# chip_one = chip_one.astype(int)
# field_two = field_two.astype(int)
# chip_two = chip_two.astype(int)
field_one = 100
chip_one = 1
field_two = 20
chip_two = 1

max_sig = 0.5
max_sep = 0.01*u.arcsec#!!!
max_deg = 3
d_pm = 150  


chip_a = 1
chip_b = 4
sig_cl = 3

ca = Table.read(pm_folder_off + 'pm_GaiaRF_ep1_f%sc%s_ep2_f%sc%sdeg%s_dmax%s_sxy%s.txt'%(field_one, chip_a, field_two, chip_two,max_deg-1,d_pm,max_sig), format = 'ascii')
cb = Table.read(pm_folder_off + 'pm_GaiaRF_ep1_f%sc%s_ep2_f%sc%sdeg%s_dmax%s_sxy%s.txt'%(field_one, chip_b, field_two, chip_two,max_deg-1,d_pm,max_sig), format = 'ascii')






ca_co = SkyCoord(ra =ca['ra1'],dec = ca['Dec1'], unit = 'degree')
cb_co = SkyCoord(ra =cb['ra1'],dec = cb['Dec1'], unit = 'degree')

idx,d2d,d3d = ca_co.match_to_catalog_sky(cb_co,nthneighbor=1)# ,nthneighbor=1 is for 1-to-1 match
sep_constraint = d2d < max_sep
ca_match = ca[sep_constraint]
cb_match = cb[idx[sep_constraint]]

fig, ax = plt.subplots(1,1)
ax.scatter(ca['ra1'],ca['Dec1'])
ax.scatter(cb['ra1'],cb['Dec1'])
ax.scatter(ca_match['ra1'],ca_match['Dec1'])

d_pmx = ca_match['pm_RA'] - cb_match['pm_RA']
d_pmy = ca_match['pm_Dec'] - cb_match['pm_Dec']

mask_x, lx_lim,hx_lim = sigma_clip(d_pmx, sigma=sig_cl, masked = True, return_bounds= True)
mask_y, ly_lim,hy_lim = sigma_clip(d_pmy, sigma=sig_cl, masked = True, return_bounds= True)
mask_xy = np.logical_and(np.logical_not(mask_x.mask), np.logical_not(mask_y.mask))

d_pmx_cl = d_pmx[mask_xy]
d_pmy_cl = d_pmy[mask_xy]


fig, (ax,ax1)= plt.subplots(1,2)
leg_x = '$\overline{\Delta \mu_{RA}}$'
leg_y = '$\overline{\Delta \mu_{Dec}}$'
lab_x = '$\Delta \mu_{RA} [mas/yr]$'
lab_y = '$\Delta \mu_{Dec} [mas/yr]$'
ax.hist(d_pmx_cl, label='%s = %.2f\n$\sigma$ = %.2f'%(leg_x,np.mean(d_pmx_cl),np.std(d_pmx_cl)))
ax.legend()
ax1.hist(d_pmy_cl, label='%s = %.2f\n$\sigma$ = %.2f'%(leg_y,np.mean(d_pmy_cl),np.std(d_pmy_cl)))
ax.legend()
ax1.legend()
ax.hist(d_pmx, histtype = 'step', ls = 'dashed', color = 'k')
ax.legend()
ax1.hist(d_pmy, histtype = 'step', ls = 'dashed', color = 'k')
ax.legend()
ax1.legend()




