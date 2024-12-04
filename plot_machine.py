#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:49:42 2024

@author: amartinez
"""

import numpy as np
from astropy import units as u
from astropy.time import Time
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from scipy.stats import sigmaclip
from astropy.stats import sigma_clip
from matplotlib import rc
from matplotlib import rcParams
import IPython
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

def diff_hist(l1_x,l2_x,l1_y,l2_y,sig_cl, clipping = None, variable = None):
    diff_x = l1_x - l2_x
    diff_y = l1_y - l2_y
    
    leg_x = '$\overline{\Delta x}$'
    leg_y = '$\overline{\Delta y}$'
    lab_x = '$\Delta x$'
    lab_y = '$\Delta y$'

    if variable is not None:
        if variable == 'coordinates':
            leg_x = '$\overline{\Delta RA}$'
            leg_y = '$\overline{\Delta Dec}$'
            lab_x = '$\Delta RA$ [mas/yr]'
            lab_y = '$\Delta Dec$ [mas/yr]'
        if variable == 'pm':
            leg_x = '$\overline{\Delta \mu_{RA}}$'
            leg_y = '$\overline{\Delta \mu_{Dec}}$'
            lab_x = '$\mu_{RA} [mas/yr]$'
            lab_y = '$\mu_{Dec} [mas/yr]$'

    # sig_cl = 3#!!!
    mask_x, lx_lim,hx_lim = sigma_clip(diff_x, sigma=sig_cl, masked = True, return_bounds= True)
    mask_y, ly_lim,hy_lim = sigma_clip(diff_y, sigma=sig_cl, masked = True, return_bounds= True)
 
   # mask_xy = mask_x & mask_y # mask_xy = np.logical(mx, my)
    mask_xy = np.logical_and(np.logical_not(mask_x.mask), np.logical_not(mask_y.mask))


    fig, (ax,ax1) = plt.subplots(1,2)
    ax.hist(diff_x, histtype = 'step', color ='#1f77b4', 
            label = '%s = %.2f\n$\sigma = %.2f$'%(leg_x, np.mean(diff_x),np.std(diff_x)) )
    ax.axvline(lx_lim, ls = 'dashed', color = 'r', lw = 3)
    ax.axvline(hx_lim, ls = 'dashed', color = 'r', lw = 3)
    ax1.hist(diff_y, histtype = 'step', color ='#ff7f0e', 
            label = '%s = %.2f\n$\sigma = %.2f$'%(leg_y, np.mean(diff_y),np.std(diff_y)) )
    ax1.axvline(ly_lim, ls = 'dashed', color = 'r', lw = 3)
    ax1.axvline(hy_lim, ls = 'dashed', color = 'r', lw = 3)
    
    
    if np.all(mask_xy) == False:
        diff_xm = diff_x[mask_xy]    
        diff_ym = diff_y[mask_xy]    
        ax.hist(diff_xm, color ='k', alpha = 0.5,
                label = '%s = %.2f\n$\sigma = %.2f$'%(leg_x,np.mean(diff_xm),np.std(diff_xm)) )
        ax1.hist(diff_ym, color ='k', alpha = 0.5,
                label = '%s = %.2f\n$\sigma = %.2f$'%(leg_y, np.mean(diff_ym),np.std(diff_ym)) )

    ax.legend()
    ax1.legend()
    ax.set_xlabel(lab_x)
    ax1.set_xlabel(lab_y)
    
    # return mask_xy