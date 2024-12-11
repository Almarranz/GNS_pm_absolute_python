#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 20:39:16 2024

@author: amartinez
"""

import numpy as np
from compare_lists import compare_lists
import Polywarp as pw
import matplotlib.pyplot as plt
import sys
from astropy.stats import sigma_clip


def gns_alignator(mode,gns_var,gns_cons, max_deg_gns, d_m,  plot = None, clipping = None, sig_clip = None):
    loop = 0
    comom_ls = []
    dic_xy = {}
    dic_Kx ={}
    dic_xy_final = {}
    
   
    mode = 'two_to_one'
    
    if mode == 'two_to_one':
       x_1 = f'x{1}'
       y_1 = f'y{1}'
       x_2 = f'x{2}'
       y_2 = f'y{2}'
       H1 = 'H1'
       H2 = 'H2'
    else:
        x_2 = f'x{1}'
        y_2 = f'y{1}'
        x_1 = f'x{2}'
        y_1 = f'y{2}'
        H1 = 'H2'
        H2 = 'H1'
        
    
    deg_gns = 1
    d_m = 15
    d_m_pm = 300
    while deg_gns < max_deg_gns:
        loop += 1 
        l1_xy = np.array([gns_cons[x_1],gns_cons[y_1]]).T
        l2_xy = np.array([gns_var[x_2],gns_var[y_2]]).T
        comp = compare_lists(l1_xy,l2_xy,d_m)
        if len(comom_ls) >1:
            if comom_ls[-1] <= comom_ls[-2]:
                try:
                    gns_var[x_2] = dic_xy[f'trans_{loop-2}'][:,0]
                    gns_var[y_2] = dic_xy[f'trans_{loop-2}'][:,1]
                    dic_xy_final['xy_deg%s'%(deg_gns)] = np.array([dic_xy[f'trans_{loop-2}'][:,0],dic_xy[f'trans_{loop-2}'][:,1]]).T            
                    comom_ls =[]
                    dic_xy = {}
                    dic_Kx = {}
                    deg_gns += 1
                    print(f'Number of common star in loop {loop-1} lower tha in loop {loop-2}.\nJupping to degree {deg_gns} ')
                    loop = -1
                    continue
                except:
                    gns_var[x_2] = dic_xy_final[f'xy_deg{deg_gns -1}'][:,0]
                    gns_var[y_2] = dic_xy_final[f'xy_deg{deg_gns -1}'][:,1]
                    print(f'Number of common star with polynomial degere {deg_gns} decreases after a single iteration.\nUsing the last iteration of degree {deg_gns -1} ')
                    break
                
        comom_ls.append(len(comp))
        print(f'Common in loop {loop}, degree {deg_gns} = %s'%(len(comp['ind_1'])))
        # if loop == 1:
        #     with open(pruebas + 'sig_and_com.txt', 'a') as file:
        #         file.write('%.1f %.0f\n'%(sig_cl, len(comp['ind_1'])))

        l1_com = gns_cons[comp['ind_1']]
        l2_com = gns_var[comp['ind_2']]
       
        
        diff_mag = l1_com[H1] - l2_com[H2] 
        
       
        # diff_mag1 = l1_com['IB230_diff'] - l2_com['IB230_diff'] 
        diff_x =  l2_com[x_2] - l1_com[x_1] 
        diff_y =  l2_com[y_2] - l1_com[y_1] 

        
        # diff_xy = (diff_x**2 + diff_y**2)**0.5
        mask_m, l_lim,h_lim = sigma_clip(diff_mag, sigma=sig_clip, masked = True, return_bounds= True)
 
        
        
        # l2_clip = l2_com
        # l1_clip = l1_com
        
        l1_clip = l1_com[np.logical_not(mask_m.mask)]
        l2_clip = l2_com[np.logical_not(mask_m.mask)]
        
        if plot == 'yes':
            fig, (ax,ax1) = plt.subplots(1,2)
            fig.suptitle(f'Degree = {deg_gns}. Loop = {loop}')
            ax.set_xlabel('$\Delta$ H')
            ax.hist(diff_mag, label = 'matches = %s\ndist = %.0f mas'%(len(comp['ind_1']), d_m))
            ax.axvline(l_lim, color = 'red', ls = 'dashed', label = '$\pm$%s$\sigma$'%(sig_clip))
            ax.axvline(h_lim, color = 'red', ls = 'dashed')
            ax.legend()
            
            ax1.hist(diff_x, label = '$\overline{\Delta x} = %.2f\pm%.2f$'%(np.mean(diff_x),np.std(diff_x)), histtype = 'step')
            ax1.hist(diff_y, label = '$\overline{\Delta y} = %.2f\pm%.2f$'%(np.mean(diff_y),np.std(diff_y)), histtype = 'step')
        
            ax1.set_xlabel('$\Delta$ pos [mas]')
            ax1.legend()
        
        
        
        
        xy_1c = np.array([l1_clip[x_1],l1_clip[y_1]]).T
        xy_2c = np.array([l2_clip[x_2],l2_clip[y_2]]).T
        

        Kx,Ky=pw.polywarp(xy_1c[:,0],xy_1c[:,1],xy_2c[:,0],xy_2c[:,1],degree=deg_gns)
        
        xi=np.zeros(len(gns_var))
        yi=np.zeros(len(gns_var))
        
        for k in range(deg_gns+1):
                    for m in range(deg_gns+1):
                        xi=xi+Kx[k,m]*gns_var[x_2]**k*gns_var[y_2]**m
                        yi=yi+Ky[k,m]*gns_var[x_2]**k*gns_var[y_2]**m
        dic_xy[f'trans_{loop+1}'] = np.array([xi,yi]).T
        
        # print(Kx[0][0])
        gns_var[x_2] = xi
        gns_var[y_2] = yi

        
