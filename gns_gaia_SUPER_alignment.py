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



field_one = 100
chip_one = 0
GNS_1='/Users/amartinez/Desktop/PhD/HAWK/GNS_1/lists/%s/chip%s/'%(field_one, chip_one)

gns1=  Table.read(GNS_1 + 'stars_calibrated_H_chip%s.ecsv'%(chip_one),  format = 'ascii.ecsv')

radec_ = SkyCoord(l = gns1['l'], b = gns1['b'], frame = 'galactic', ).fk5
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source" # Select early Data Release 3
# Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source" # Select early Data Release 3
Gaia.ROW_LIMIT = -1  # Ensure the default row limit.
coord = SkyCoord(ra=np.mean(radec_.ra), dec=np.mean(radec_.dec), unit='degree', frame='icrs',equinox ='J2000',obstime='J2016.0')




e_pm = 0.3

# j = Gaia.cone_search_async(coord, )
# coord = SkyCoord(ra=10.625, dec=41.2, unit=(u.degree, u.degree), frame='icrs')
# j = Gaia.cone_search_async(coord,100*u.arcsec)
# j = Gaia.cone_search_async(coord,100*u.arcsec)
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

gaia_coord = SkyCoord(l = gaia['l'].value*u.deg, b = gaia['b'].value*u.deg, frame = 'galactic')

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


max_sep = 0.5*u.arcsec
ga_c = SkyCoord(l = gaia['l'], b = gaia['b'], unit = 'degree',frame = 'galactic' )
gn_c= SkyCoord(l= gns1['l'], b=gns1['b'],unit = 'degree', frame = 'galactic')
# gns1_coor_match = gns1_coor_fg.transform_to('icrs')
idx,d2d,d3d = ga_c.match_to_catalog_sky(gn_c,nthneighbor=1)# ,nthneighbor=1 is for 1-to-1 match
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


gn_m = gn_m[mask_lb]
ga_m = ga_m[mask_lb]


radec_ = radec_.transform_to('icrs')
ragai1_off,decgai1_off = radec_.spherical_offsets_to(GaiaCoord.frame)

xy_gn = np.array([gn_m['x'],gn_m['y']]).T
xy_ga = np.array([ga_m['x'],ga_m['y']]).T



















































