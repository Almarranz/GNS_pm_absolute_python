#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 17:45:11 2024

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


# Pixel scale = 0.0135"/pixel



field = 2

# folder_ls = '/Users/amartinez/Desktop/for_people/for_Laly/fields/field_%s/'%(field)
folder_ls = '/Users/amartinez/Desktop/for_people/for_Laly/fields/field%s/'%(field)
pruebas = '/Users/amartinez/Desktop/for_people/for_Laly/pruebas/'

lw = 'IB224'
# l_a = 'IB_calib_field2_2011.txt'
# l_a = 'IB_calib_field2_2011.txt'
l_a = '2011/%s_calib_field%s_2011.txt'%(lw,field)
l_b = '2013/%s_calib_field%s_2013.txt'%(lw,field)

l1 = Table.read(folder_ls +l_a, format = 'ascii') 
l2 = Table.read(folder_ls +l_b, format = 'ascii') 
sig_cl_ls = np.arange(1,3,0.1)

# with open(pruebas + 'sig_and_com.txt', 'w') as file:
#     file.write('# sig_cl num of mathches\n')


    

sig_cl =3
d_m = 2
deg = 1
print(f'Degree = {deg} sig = {sig_cl}')
comom_ls = []
dic_xy = {}
loop = -1
max_deg = 2
# for loop in range(20):
while deg < max_deg:
    loop += 1 
    l1_xy = np.array([l1['x'],l1['y']]).T
    l2_xy = np.array([l2['x'],l2['y']]).T
    comp = compare_lists(l1_xy,l2_xy,d_m)
    if len(comom_ls) >1:
        if comom_ls[-1] <= comom_ls[-2]:
            l1['x'] = dic_xy[f'trans_{loop-2}'][:,0]
            l1['y'] = dic_xy[f'trans_{loop-2}'][:,1]
            comom_ls =[]
            dic_xy = {}
            deg += 1
            loop = -1
            continue
            
    comom_ls.append(len(comp))
    print(f'Common in loop {loop}, degree {deg} = %s'%(len(comp['ind_1'])))
    # if loop == 1:
    #     with open(pruebas + 'sig_and_com.txt', 'a') as file:
    #         file.write('%.1f %.0f\n'%(sig_cl, len(comp['ind_1'])))

    l1_com = l1[comp['ind_1']]
    l2_com = l2[comp['ind_2']]
    
    diff_mag = l1_com['%s_diff'%(lw)] - l2_com['%s_diff'%(lw)] 
    # diff_mag1 = l1_com['IB230_diff'] - l2_com['IB230_diff'] 
    diff_x =  l2_com['x'] - l1_com['x'] 
    diff_y =  l2_com['y'] - l1_com['y'] 
    diff_xy = (diff_x**2 + diff_y**2)**0.5
    mask_m, l_lim,h_lim = sigma_clip(diff_mag, sigma=sig_cl, masked = True, return_bounds= True)
    
    
    
    # l2_clip = l2_com
    # l1_clip = l1_com
    
    l1_clip = l1_com[~mask_m.mask]
    l2_clip = l2_com[~mask_m.mask]
    
    fig, (ax,ax1) = plt.subplots(1,2)
    fig.suptitle(f'Degree = {deg}. Loop = {loop}')
    ax.set_xlabel('$\Delta$ IB224')
    ax.hist(diff_mag, label = 'matches = %s\ndist = %.2f arcsec'%(len(comp['ind_1']), d_m*0.0135))
    ax.axvline(l_lim, color = 'red', ls = 'dashed', label = '$\pm$%s$\sigma$'%(sig_cl))
    ax.axvline(h_lim, color = 'red', ls = 'dashed')
    ax.legend()
    
    ax1.hist(diff_x, label = '$\overline{\Delta x} = %.2f\pm%.2f$'%(np.mean(diff_x),np.std(diff_x)), histtype = 'step')
    ax1.hist(diff_y, label = '$\overline{\Delta y} = %.2f\pm%.2f$'%(np.mean(diff_y),np.std(diff_y)), histtype = 'step')

    # ax1.hist(diff_xy, label = 'Distance',color = 'grey', alpha = 0.3)
    
    # ax1.hist(comp['dist'], color = 'pink', alpha = 0.5)
    # ax1.hist(comp['dist'], color = 'red',  histtype = 'step', lw=2)
    ax1.set_xlabel('$\Delta$ pixel')
    # ax1.set_xlabel('$\Delta$ IB230')
    # ax1.hist(diff_mag1)
    ax1.legend()
    
    
    
    
    xy_1c = np.array([l1_clip['x'],l1_clip['y']]).T
    xy_2c = np.array([l2_clip['x'],l2_clip['y']]).T
    
    # Kx,Ky=pw.polywarp(xy_1c[:,0],xy_1c[:,1],xy_2c[:,0],xy_2c[:,1],degree=1)
    Kx,Ky=pw.polywarp(xy_2c[:,0],xy_2c[:,1],xy_1c[:,0],xy_1c[:,1],degree=deg)
    
    xi=np.zeros(len(l1))
    yi=np.zeros(len(l1))
    
    for k in range(deg+1):
                for m in range(deg+1):
                    xi=xi+Kx[k,m]*l1['x']**k*l1['y']**m
                    yi=yi+Ky[k,m]*l1['x']**k*l1['y']**m
    # dic_xy[f'trans_{loop+1}'] = np.array([xi,yi]).T
    dic_xy[f'trans_{loop+1}'] = np.array([xi,yi]).T
    l1['x'] = xi
    l1['y'] = yi
# sys.exit()
# =============================================================================
# fig, (ax,ax1) = plt.subplots(1,2)
# fig.suptitle(f'Degree = {deg}. Loop = {loop}')
# ax.set_xlabel(f'$\Delta$ {lw}')
# ax.hist(diff_mag, label = 'matches = %s\ndist = %.2f arcsec'%(len(comp['ind_1']), d_m*0.0135))
# ax.axvline(l_lim, color = 'red', ls = 'dashed', label = '$\pm$%s$\sigma$'%(sig_cl))
# ax.axvline(h_lim, color = 'red', ls = 'dashed')
# ax.legend()
# 
# ax1.hist(diff_x, label = '$\overline{\Delta x} = %.2f\pm%.2f$'%(np.mean(diff_x),np.std(diff_x)), histtype = 'step')
# ax1.hist(diff_y, label = '$\overline{\Delta y} = %.2f\pm%.2f$'%(np.mean(diff_y),np.std(diff_y)), histtype = 'step')
# 
# # ax1.hist(diff_xy, label = 'Distance',color = 'grey', alpha = 0.3)
# 
# # ax1.hist(comp['dist'], color = 'pink', alpha = 0.5)
# # ax1.hist(comp['dist'], color = 'red',  histtype = 'step', lw=2)
# ax1.set_xlabel('$\Delta$ pixel')
# # ax1.set_xlabel('$\Delta$ IB230')
# # ax1.hist(diff_mag1)
# ax1.legend()
# =============================================================================

# %%
l1_xy = np.array([l1['x'],l1['y']]).T
l2_xy = np.array([l2['x'],l2['y']]).T
comp = compare_lists(l1_xy,l2_xy,d_m)
print(len(comp))    

comp_pm = compare_lists(l1_xy,l2_xy,3)

pm_x = (l2['x'][comp_pm['ind_2']] - l1['x'][comp_pm['ind_1']])*0.0135/2*1000
pm_y = (l2['y'][comp_pm['ind_2']] - l1['y'][comp_pm['ind_1']])*0.0135/2*1000
l1_c = l1[comp_pm['ind_1']]
l2_c = l2[comp_pm['ind_2']]

fig, (ax, ax1) = plt.subplots(1,2)

ax.hist(pm_x,bins = 'auto',histtype = 'step', label = '$\mu_{x} =%.2f \pm %.2f $'%(np.mean(pm_x),np.std(pm_x)))
ax1.hist(pm_y,bins = 'auto',histtype = 'step', label = '$\mu_{y} =%.2f \pm %.2f $'%(np.mean(pm_y),np.std(pm_y)))
ax.set_xlabel('$\mu_{x}$ [mas]')
ax1.set_xlabel('$\mu_{y}$ [mas]')
ax.legend()
ax1.legend()
# %%
fig, (ax,ax1) = plt.subplots(1,2)
cb = ax.hist2d(pm_x, pm_y, bins = 20,norm=colors_plt.LogNorm())
cb1 = ax1.hist2d(l2['x'], l2['y'], bins =20    ,norm=colors_plt.LogNorm())
# ax.scatter(pm_x,pm_y,s=1,color = 'r')
fig.colorbar(cb[3], ax = ax)
fig.colorbar(cb1[3], ax = ax1)
ax.set_xlabel('$\mu_{x} [mas]$')
ax.set_ylabel('$\mu_{y} [mas]$')
ax1.set_ylabel('$y [pix]$')
ax1.set_xlabel('$x [pix]$')
# %%
knn = 10
sim_by= 'kernnel'
lim = 'minimun'
modes = ['pm_xy_color','pm_xy','pm_color','pm']
# mode = 'pm_xy'
# clst = cluster_finder.finder(pm_x, pm_y, l2_c['x'], l2_c['y'], l2_c['RA'], l2_c['Dec'], modes[1], 
#                              l2_c['H_diff'],l2_c['IB236_diff'],
#                              knn,sim_by,lim)
clst = cluster_finder.finder(pm_x, pm_y, l2_c['x'], l2_c['y'], l2_c['RA'], l2_c['Dec'], modes[1], 
                              l2_c['%s_diff'%(lw)],l2_c['%s_diff'%(lw)],
                              knn,sim_by,lim)

# cluster_finder.finder(pmra, pmdec, x, y, Ra, Dec, clustered_by, color_A, color_B, samples_dist, gen_sim, sim_lim)


# %%
# =============================================================================
# # %%
# # GAIA comparison
# 
# ra_m = np.mean(l1['RA'])
# dec_m = np.mean(l1['Dec'])
# 
# cent = SkyCoord(ra_m, dec_m, unit = 'degree')
# 
# 
# 
# ima = folder_ls + 'IB224_astr_field2_2011.fits'
# 
# head = fits.getheader(ima)
# width = head['NAXIS1']*0.0135*u.arcsec
# height = head['NAXIS2']*0.0135*u.arcsec
# 
# Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source" # Select early Data Release 3
# Gaia.ROW_LIMIT = -1  # Ensure the default row limit.
# coord = SkyCoord(ra=cent.ra, dec=cent.dec, unit='degree', frame='icrs',equinox ='J2000',obstime='J2011.0')
# 
# 
# j = Gaia.cone_search_async(coord, np.sqrt((width)**2+height**2))
# 
# gaia = j.get_results()
# # %%
# with open(pruebas + f'gaia_field{field}.txt','w' ) as f:
#     f.write('')
# 
# with open(pruebas + f'gaia_field{field}.txt','w' ) as f:
#     for i in range(len(gaia)):
#         f.write('%s %s\n'%(gaia['ra'][i],gaia['dec'][i]))
#     
# =============================================================================






































