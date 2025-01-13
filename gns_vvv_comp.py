#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 10:15:54 2025

@author: amartinez
"""

import numpy as np
import matplotlib.pyplot as plt
from compare_lists import compare_lists 
from astropy.table import Table
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
import matplotlib.colors as colors_plt
from skimage import data
from skimage import transform
import astroalign as aa
from astroquery.gaia import Gaia
import skimage as ski
import sys
from astropy.stats import sigma_clip
from astropy.io import fits
from astropy.coordinates import SkyCoord
import Polywarp as pw
from astroquery.gaia import Gaia
from astropy import units as u
import cluster_finder
import pandas as pd
import copy
import cluster_finder
from filters import filter_gaia_data
from filters import filter_hosek_data
from filters import filter_gns_data
from filters import filter_vvv_data
# %% 
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
# %%
pm_folder = '/Users/amartinez/Desktop/PhD/HAWK/pm_gns1_gns2_absolute_python/'

field_one, chip_one, field_two, chip_two,t1,t2,max_sig = np.loadtxt('/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/lists/fields_and_chips.txt', 
                                                       unpack=True)
field_one = field_one.astype(int)
chip_one = chip_one.astype(int)
field_two = field_two.astype(int)
chip_two = chip_two.astype(int)

# field_one = 7
# chip_one = 2
# field_two = 7
# chip_two = 1
# max_sig = 0.5#TODO

max_deg = 3
d_pm = 150
gns1_pm = Table.read(pm_folder + 'pm_GaiaRF_ep1_f%sc%s_ep2_f%sc%sdeg%s_dmax%s_sxy%s.txt'%(field_one, chip_one, field_two, chip_two,max_deg-1,d_pm,max_sig), format = 'ascii')


#Gaia comparison

search_r = 50*u.arcsec

ra_c = np.mean(gns1_pm['ra1'])
dec_c = np.mean(gns1_pm['Dec1'])


# Comparison with VVV 
VVV = Table.read('/Users/amartinez/Desktop/PhD/Catalogs/VVV/b333/PMS/b333.dat', format = 'ascii')


gns1_coor = SkyCoord(ra = gns1_pm['ra1'],dec = gns1_pm['Dec1'],unit = 'degree',
                     frame = 'fk5', obstime = 'J2015.4301')
rad = 200/3600

delta_ra = (VVV['ra'] - ra_c) * np.cos(np.radians(dec_c))
delta_dec = VVV['dec'] - dec_c

# Calculate the angular distance
angular_distance = np.sqrt(delta_ra**2 + delta_dec**2)

# Select rows where the distance is within the radius
within_radius = angular_distance <= rad

vvv_c = VVV[within_radius]

Ks_max = None
Ks_min = None
vvv_c = filter_vvv_data(vvv_c,
                    pmRA = 'good',
                    pmDE = None,
                    epm = 0.95,
                    ok = 'yes',
                    max_Ks = Ks_max,
                    min_Ks = Ks_min,
                    center = 'yes'
                    )



fig, ax = plt.subplots(1,1)
ax.scatter(vvv_c['ra'],vvv_c['dec'])
ax.scatter(gns1_pm['ra1'], gns1_pm['Dec1'], alpha = 0.1)

# Moves coordinates to GNS1 obstime???
tvvv = 2012.29578304
# vvv_c['ra'] = vvv_c['ra'] + (t1-tvvv) * (vvv_c['pmRA']/1000.0)/3600.0 * np.cos(vvv_c['dec']*np.pi/180.)
# vvv_c['ra'] = vvv_c['ra'] + (t1-tvvv) * (vvv_c['pmRA']/1000.0)/3600.0 
# vvv_c['dec'] = vvv_c['dec'] + (t1-tvvv) * (vvv_c['pmDEC']/1000.0)/3600.0
vvv_coord = SkyCoord(ra = vvv_c['ra'], dec = vvv_c['dec'], unit = 'degree',
                      frame = 'icrs', obstime = 'J2012.29578304')
# vvv_coord = SkyCoord(ra = vvv_c['ra'], dec = vvv_c['dec'], unit = 'degree',
#                      frame = 'icrs', obstime = 'J2015.4301')


max_sep = 0.035*u.arcsec#!!!

idx,d2d,d3d = vvv_coord.match_to_catalog_sky(gns1_coor,nthneighbor=1)# ,nthneighbor=1 is for 1-to-1 matchsep_constraint = d2d < max_sep
sep_constraint = d2d < max_sep
gns_vvv = gns1_pm[idx[sep_constraint]]
vvv_match = vvv_c[sep_constraint]



sig_cl = 3
dKs = gns_vvv['Ks1']-vvv_match['Ks']
mask_m, l_lim,h_lim = sigma_clip(dKs, sigma=sig_cl, masked = True, return_bounds= True)


fig, (ax, ax1) = plt.subplots(1,2)
ax.scatter(gns_vvv['ra1'],gns_vvv['Dec1'])
ax.scatter(vvv_match['ra'],vvv_match['dec'],s =1)
ax.set_title('Macthes = %s'%(len(vvv_match)))
ax1.hist(dKs,bins = 'auto', label = '$\Delta$Ks = %.2f\n$\sigma$ = %.2f'%(np.mean(dKs),np.std(dKs)))
ax1.axvline(l_lim, ls = 'dashed', color = 'r', label = '$\pm$ %s$\sigma$'%(sig_cl))
ax1.axvline(h_lim, ls = 'dashed', color = 'r')
ax1.legend(loc = 4, fontsize = 12)

vvv_match = vvv_match[np.logical_not(mask_m.mask)]
gns_vvv = gns_vvv[np.logical_not(mask_m.mask)]
ax1.set_title(f'{sig_cl}$\sigma$ Macthes= %s'%(len(vvv_match)))




diff_pmx = gns_vvv['pm_RA'] - vvv_match['pmRA']
diff_pmy = gns_vvv['pm_Dec'] - vvv_match['pmDEC']




mask_pmx, l_lim,h_lim = sigma_clip(diff_pmx, sigma=3, masked = True, return_bounds= True)
mask_pmy, l_lim,h_lim = sigma_clip(diff_pmy, sigma=3, masked = True, return_bounds= True)

mask_pm = np.logical_and(np.logical_not(mask_pmx.mask),np.logical_not(mask_pmy.mask))
diff_pmx_clip = diff_pmx[mask_pm]
diff_pmy_clip = diff_pmy[mask_pm]

# %
fig, (ax,ax1) = plt.subplots(1,2)
ax.set_title('Matching = %s'%(len(gns_vvv['pm_RA'])))
ax1.set_title('VVV comparison')
ax.hist(diff_pmx_clip, histtype = 'step', label = '$\overline{\Delta x}$ =%.2f\n$\sigma$ = %.2f '%(np.nanmean(diff_pmx_clip),np.nanstd(diff_pmx_clip)))
ax.hist(diff_pmx, histtype = 'step', color ='k', alpha = 0.3)

ax1.hist(diff_pmy_clip, histtype = 'step',color = 'orange', label = '$\overline{\Delta x}$ =%.2f\n$\sigma$ = %.2f '%(np.nanmean(diff_pmy_clip),np.nanstd(diff_pmy_clip)))
ax1.hist(diff_pmy, histtype = 'step',color = 'k', alpha = 0.3, ls = 'dashed')
ax.legend()
ax1.legend()
ax.set_xlabel('$\Delta \mu_{\parallel}$')
ax1.set_xlabel('$\Delta \mu_{\perp}$')




















