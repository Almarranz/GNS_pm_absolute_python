#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 10:30:27 2022

@author: amartinez
"""

# Generates offsets for Gaia stars and astroaligned with xy coordinates in GNS1

# Here we are going to align GNS (1 and 2) to Gaia reference frame for the 
# each of tha gns epochs
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
import astroquery
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
from scipy import stats
from matplotlib import colors 
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



field_one = 60
chip_one = 0
GNS_1='/Users/amartinez/Desktop/PhD/HAWK/GNS_1/lists/%s/chip%s/'%(field_one, chip_one)

gns1=  Table.read(GNS_1 + 'stars_calibrated_H_chip%s.ecsv'%(chip_one),  format = 'ascii.ecsv')
factor = -1
pix_scale =0.5*0.1064
p_mas = pix_scale*1000#Trnasform pixeles into mas
gns1['x'] = gns1['x']*factor*p_mas
gns1['y'] = gns1['y']*p_mas

min_exp = gns1['nexp'] > np.max(gns1['nexp'])*0.25
gns1 = gns1[min_exp]

num_bins = 200
s_xy = np.sqrt(gns1['sl']**2 + gns1['sb']**2)

# Create a meshgrid for plotting

# %
# Plot the result
fig, (ax,ax2) = plt.subplots(2,1, figsize =(6,6))
statistic, x_edges, y_edges, binnumber = stats.binned_statistic_2d(gns1['l'], gns1['b'], s_xy, statistic='median', bins=(num_bins,int(num_bins/2)))
X, Y = np.meshgrid(x_edges, y_edges)

# c = ax.pcolormesh(X, Y, statistic.T, cmap='Spectral',norm = colors.LogNorm())
c = ax.pcolormesh(X, Y, statistic.T, cmap='Spectral_r', vmin = 0.01,vmax = 0.05)
cb1 = fig.colorbar(c, ax=ax, label=r'$\sqrt{sl^2 +sb^2}$ [arcsec]', shrink = 0.8,)
# ax.colorbar(c, ax=ax, label=r'$\sqrt{sl^2 +sb^2}$ [arcsec]', shrink = 0.5)
ax.set_xlabel('l(º)')
ax.set_ylabel('b (º)')
ax.axis('scaled')

statistic, x_edges, y_edges, binnumber = stats.binned_statistic_2d(gns1['x'], gns1['y'], s_xy, statistic='median', bins=(num_bins,int(num_bins/2)))
X, Y = np.meshgrid(x_edges, y_edges)

c2 = ax2.pcolormesh(X, Y, statistic.T, cmap='Spectral_r',vmin = 0.01,vmax = 0.05)
cb2  = fig.colorbar(c2, ax=ax2, label=r'$\sqrt{sl^2 +sb^2}$ [arcsec]', shrink = 0.8, )
# ax2.colorbar(c2, ax=ax, label=r'$\sqrt{sl^2 +sb^2}$ [arcsec]', shrink = 0.5)
ax2.set_xlabel('x [mas]')
ax2.set_ylabel('y[mas]')
ax2.axis('scaled')



fig,ax = plt.subplots(1,1)
N=1
# ax.scatter(gns1['x'][::N],gns1['y'][::N])
ax.scatter(gns1['l'][::N],gns1['b'][::N])

mx = (np.min(gns1['x']) + np.max(gns1['x']))/2
my = (np.min(gns1['y']) + np.max(gns1['y']))/2

chip1 = (gns1['x'] > np.min(gns1['x'])) & (gns1['x'] < mx*1.1) & (gns1['y'] >0) & (gns1['y'] <my*0.9) 
chip2 = (gns1['x'] > mx*0.95) & (gns1['x'] < mx*-2) & (gns1['y'] >0) & (gns1['y'] <my*0.9) 
chip3 = (gns1['x'] > mx*0.95) & (gns1['x'] < mx*-2) & (gns1['y'] >my*1.1) & (gns1['y'] <my*2) 
chip4 = (gns1['x'] > np.min(gns1['x']))  & (gns1['x'] < mx*1.1) & (gns1['y'] >my*1.1) & (gns1['y'] <my*2) 

all_chips= chip1 | chip2 | chip3 | chip4
# all_chips= chip2 | chip3

# gns1 = gns1[chip1]
gns1 = gns1[all_chips]
# gns1 = gns1[np.logical_not(all_chips)]


# ax.scatter(gns1['x'][::N],gns1['y'][::N],s =3)
ax.scatter(gns1['l'][::N],gns1['b'][::N],s =3)
ax.axis('equal')

radec_ = SkyCoord(l = gns1['l'], b = gns1['b'], frame = 'galactic', ).fk5
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source" # Select early Data Release 3
# Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source" # Select early Data Release 3
Gaia.ROW_LIMIT = -1  # Ensure the default row limit.
coord = SkyCoord(ra=np.mean(radec_.ra), dec=np.mean(radec_.dec), unit='degree', frame='icrs',equinox ='J2000',obstime='J2016.0')


transf = 'affine'#!!!
# transf = 'similarity'#!!!1
# transf = 'polynomial'#!!!

order_trans = 1
e_pm = 0.3#!!!
max_sep = 0.1*u.arcsec#!!!Frist Macthig distance with Gaia 
d_m = 200 #!!! distance in mas using for compare_lists
max_deg = 3#!!! Maximun degree using for the interative polynomial alignment
clip_in_alig = 'yes'#!!! Perform clipping inside the interative polynomial alignment
bad_sig  = 3
align_1 = 'Polywarp'
# align_1 = '2DPoly'
f_mode = 'NWnC'

size = 0.35*u.deg
gaia= Gaia.query_object_async(coord, width= 0.5*size, height=size)
# gaia_ = Gaia.cone_search((coord.ra, coord.dec),0.5*size)
# gaia_ = j.get_results()
id_gaia = np.arange(len(gaia))

gaia['id'] = id_gaia
# If you get gaia_ from the Gaia website use "duplicated_source= False"
gaia_good = filter_gaia_data(
    gaia_table=gaia,
    astrometric_params_solved=31,
    duplicated_source= False,
    parallax_over_error_min=-10,
    astrometric_excess_noise_sig_max=2,
    phot_g_mean_mag_min= None,
    phot_g_mean_mag_max= 13,
    pm_min=0,
    pmra_error_max=None,
    pmdec_error_max=None
    )

gaia = gaia_good

gaia_coord = SkyCoord(l = gaia['l'], b = gaia['b'], frame = 'galactic')
coord = coord.transform_to('galactic')
ragai1_off,decgai1_off = coord.spherical_offsets_to(gaia_coord.frame)

gaia['x'] = ragai1_off.to(u.mas)
gaia['y'] = decgai1_off.to(u.mas)


lw_g = gaia_coord.l.wrap_at('361d').degree
# lw_g = gaia_coord.l.wrap_at(360)

# %%
fig, ax = plt.subplots(1,1)
ax.scatter(gns1['l'],gns1['b'], label = 'GNS1')
# ax.scatter(lw_g, gaia['b'])


cut = (gaia['l']<max(gns1['l']))&(gaia['l']>min(gns1['l']))&(gaia['b']<max(gns1['b']))&(gaia['b']>min(gns1['b']))
gaia = gaia[cut]

ax.scatter(gaia['l'],gaia['b'], label =f'Gaia = {len(gaia)}') 
ax.axis('scaled')



ga_c = SkyCoord(l = gaia['l'], b = gaia['b'], unit = 'degree',frame = 'galactic' )
gn_c= SkyCoord(l= gns1['l'], b=gns1['b'],unit = 'degree', frame = 'galactic')
# gns1_coor_match = gns1_coor_fg.transform_to('icrs')
idx,d2d,d3d = ga_c.match_to_catalog_sky(gn_c)# ,nthneighbor=1 is for 1-to-1 match
sep_constraint = d2d < max_sep
ga_m = gaia[sep_constraint]
gn_m = gns1[idx[sep_constraint]]

ax.scatter(ga_m['l'],ga_m['b'],s = 3, label = f'Gaia matches = {len(ga_m)}')

ax.legend()


diff_l = (gn_m['l'] - ga_m['l'])*u.deg.to(u.arcsec)
diff_b = (gn_m['b'] - ga_m['b'])*u.deg.to(u.arcsec)

mask_l, l_liml,h_liml = sigma_clip(diff_l, sigma=3, masked = True, return_bounds= True)
mask_b, l_limb,h_limb = sigma_clip(diff_b, sigma=3, masked = True, return_bounds= True)

mask_lb = np.logical_and(np.logical_not(mask_l.mask), np.logical_not(mask_b.mask))

diff_lc= diff_l[mask_lb]
diff_bc= diff_b[mask_lb]
# %%
fig, (ax, ax2) = plt.subplots(1,2)
ax.set_title(f'Matches = {len(diff_lc)}')
ax.hist(diff_lc, histtype = 'step',linewidth =1, label = fr'$\Delta$ l = {np.mean(diff_lc): .2f}$\pm$ {np.std(diff_lc): .2f} arcsec' )
ax2.hist(diff_bc, histtype = 'step',linewidth =1,label = fr'$\Delta$ b = {np.mean(diff_bc): .2f}$\pm$ {np.std(diff_bc): .2f} arcsec')
ax.hist(diff_l,color = 'k', alpha = 0.05)
ax2.hist(diff_b,color = 'k', alpha = 0.05)
ax.axvline(l_liml, color = 'r', ls = 'dashed')
ax.axvline(h_liml, color = 'r', ls = 'dashed')
ax2.axvline(l_limb, color = 'r', ls = 'dashed')
ax2.axvline(h_limb, color = 'r', ls = 'dashed')
ax.legend()
ax2.legend()
ax.set_xlabel(r'$\Delta l$[arcsec]')
ax2.set_xlabel(r'$\Delta b$[arcsec]')

# %%


# gn_m = gn_m[mask_lb]
# ga_m = ga_m[mask_lb]


# =============================================================================
# xy_gn = np.array([gn_m['x'],gn_m['y']]).T
# xy_ga = np.array([ga_m['x'],ga_m['y']]).T
# =============================================================================
xy_gn = np.array([gn_m['l'],gn_m['b']]).T
xy_ga = np.array([ga_m['l'],ga_m['b']]).T

# %%
fig, ax = plt.subplots(1,1)
ax.scatter(xy_gn[:,0],xy_gn[:,1])
ax.scatter(xy_ga[:,0],xy_ga[:,1],s = 1)
ax.axis('equal')



if transf == 'polynomial':
    p = ski.transform.estimate_transform(transf,
                                        xy_gn, 
                                        xy_ga, order = order_trans)
else:    
    p = ski.transform.estimate_transform(transf,
                                    xy_gn, 
                                    xy_ga)



print(p)

# %%
fig, ax = plt.subplots(1,1)
# ax.scatter(gn_m['x'],gn_m['y'], marker = '*', label = 'GNS1')
# ax.scatter(ga_m['x'],ga_m['y'], marker = '*', label = 'Gaia')
ax.scatter(gn_m['l'],gn_m['b'], marker = '*', label = 'GNS1')
ax.scatter(ga_m['l'],ga_m['b'], marker = '*', label = 'Gaia')
xy_gn_t = p(xy_gn)
ax.scatter(xy_gn_t[:,0],xy_gn_t[:,1], label = 'GNS1_t',s=2)
ax.legend()
gn_xy = np.array([gns1['x'],gns1['y']]).T
gn_xyt = p(gn_xy)
gns1['x'] = gn_xyt[:,0] 
gns1['y'] = gn_xyt[:,1] 
gns1['x'] = gn_xyt[:,0] 
gns1['y'] = gn_xyt[:,1] 
ax.axis('equal')


sys.exit(278)
s_ls = compare_lists(np.array([gns1['x'],gns1['y']]).T, np.array([gaia['x'],gaia['y']]).T, d_m)
# %%
gns1 = alignator(1, gns1, gaia, s_ls, d_m, max_deg =max_deg, clipping = clip_in_alig, align_by = align_1, f_mode = f_mode, plot = 'yes')




































