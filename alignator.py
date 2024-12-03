#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 11:44:21 2024

@author: amartinez
"""
import numpy as np
from compare_lists import compare_lists
import Polywarp as pw
import matplotlib.pyplot as plt
import sys
from astropy.stats import sigma_clip

def alignator(survey,gns,gaia,s_ls, d_m,max_deg, deg):
    # Proceed with the iterative alignment process from an initial common lists
    
    id1 = f'x{survey}'
    id2 = f'y{survey}'
    loop =1
    deg = 1

    comom_ls =[]
    comom_ls.append(len(s_ls))
    dic_xy = {} 
    dic_xy['trans_0'] = np.array([gns[id1],gns[id2]]).T
    # for loop in range(1,10):
    max_deg = 4#!!!
    sig_cl = 3#!!!
    while deg < max_deg:
        
        Kx,Ky=pw.polywarp(s_ls['l2_x'],s_ls['l2_y'],s_ls['l1_x'],s_ls['l1_y'],degree=deg)
        
        xi=np.zeros(len(gns))
        yi=np.zeros(len(gns))
        
        for k in range(deg+1):
                    for m in range(deg+1):
                        xi=xi+Kx[k,m]*gns[id1]**k*gns[id2]**m
                        yi=yi+Ky[k,m]*gns[id1]**k*gns[id2]**m
        dic_xy[f'trans_{loop}'] = np.array([xi,yi]).T
        
        # print(Kx[0][0])
        gns[id1] = xi
        gns[id2] = yi
        
        # l1_xy = np.array([gns[id1],gns[id2]]).T
        # g1_xy = np.array([gaia['x'], gaia['y']]).T
        
        l_xy = np.array([gns[id1],gns[id2]]).T
        ga_xy = np.array([gaia['x'], gaia['y']]).T
        
        
        
        s_ls = compare_lists(l_xy,ga_xy,d_m)
        print(f'\nCommon GNS1 and Gaia after loop{loop} = {len(s_ls)}')
        
        diff_x =s_ls['l2_x'] - s_ls['l1_x']  
        diff_y =s_ls['l2_y'] - s_ls['l1_y']  
        
        mask_x, lx_lim,hx_lim = sigma_clip(diff_x, sigma=sig_cl, masked = True, return_bounds= True)
        mask_y, ly_lim,hy_lim = sigma_clip(diff_y, sigma=sig_cl, masked = True, return_bounds= True)

        
        comom_ls.append(len(s_ls))
        if comom_ls[-1] <= comom_ls[-2] :
            if len(comom_ls)>2:
                
                gns[id1] = dic_xy[f'trans_{loop-1}'][:,0]
                gns[id2] = dic_xy[f'trans_{loop-1}'][:,1]
                deg +=1
                print(30*'-'+f'\nBreaking after loop = {loop-1}\nStarting alignment degree = {deg}\n' + 30*'-')
                continue
            else:
                print(f'Polynomial of degree {deg} does not work')
                gns[id1] = dic_xy['trans_0'][:,0]
                gns[id2] = dic_xy['trans_0'][:,1]
                break
           
        
        
        fig, (ax,ax1)  = plt.subplots(1,2)
        ax1.set_title('Degree = %s'%(deg))
        ax.set_title(f'Loop = {loop}. Matching stars = {len(s_ls)}')
        ax.hist(diff_x, histtype = 'step', label = '$\overline{x} = %.2f$\n$\sigma x$ =%.2f'%(np.mean(diff_x),np.std(diff_x)))
        ax.axvline(lx_lim, color = 'red', ls = 'dashed', label = '$\pm$%s$\sigma$'%(sig_cl))
        ax.axvline(hx_lim, color = 'red', ls = 'dashed')
        ax1.axvline(ly_lim, color = 'red', ls = 'dashed')
        ax1.axvline(hy_lim, color = 'red', ls = 'dashed')
        ax.legend()
        ax1.hist(diff_y, histtype = 'step', label = '$\overline{y} = %.2f$\n$\sigma x$ =%.2f'%(np.mean(diff_y),np.std(diff_y)))
        ax1.legend()
        
        loop +=1     
    return gns