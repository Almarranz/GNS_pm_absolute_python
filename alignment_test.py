#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 10:58:23 2024

@author: amartinez
"""

import numpy as np
from compare_lists import compare_lists
import Polywarp as pw
import matplotlib.pyplot as plt
import sys
from astropy.stats import sigma_clip
from astropy.modeling.models import Polynomial2D
from astropy.modeling.fitting import LinearLSQFitter
from astropy.modeling import models, fitting
from astropy.table import Table, hstack

import numpy as np

field_one, chip_one, field_two, chip_two,t1,t2,max_sig = np.loadtxt('/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/lists/fields_and_chips.txt', 
                                                       unpack=True)
field_one = field_one.astype(int)
chip_one = chip_one.astype(int)
field_two = field_two.astype(int)
chip_two = chip_two.astype(int)


# # Number of points
# num_points = 50000

# # Original distribution in the x, y plane
# x = np.random.uniform(-100000, 100000, num_points)
# y = np.random.uniform(-50000, 50000, num_points)

GNS_1off='/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/lists/%s/chip%s/'%(field_one, chip_one)
GNS_2off='/Users/amartinez/Desktop/PhD/HAWK/GNS_2absolute_python/lists/%s/chip%s/'%(field_two, chip_two)
gns2 = Table.read(GNS_2off + 'stars_calibrated_H_chip%s_on_gns1_f%sc%s_sxy%s.txt'%(chip_two,field_one,chip_one,max_sig), format = 'ascii')
gns1 = Table.read(GNS_2off + 'stars_calibrated_H_chip%s_on_gns1_f%sc%s_sxy%s.txt'%(chip_two,field_one,chip_one,max_sig), format = 'ascii')



num_points = len(gns2)

# Displacements: Gaussian distribution
mean_shift = 0  # Mean of the Gaussian
std_shift = 15  # Standard deviation of the Gaussian
shifts = np.random.normal(mean_shift, std_shift, num_points)

# Angles for displacement direction
angles = np.random.uniform(0, 2 * np.pi, num_points)

# Calculate the displacement components
dx = shifts * np.cos(angles)
dy = shifts * np.sin(angles)

# New distribution with Gaussian shifts
x_sh = gns2['x2'] + dx
y_sh = gns2['y2'] + dy


gns1['x2'] = x_sh
gns1['y2'] = y_sh

fig, ax = plt.subplots(1,1)
ax.scatter(gns2['x2'][::100],gns2['y2'][::100] , label = 'Original') 
ax.scatter(gns1['x2'][::100],gns1['y2'][::100] ,marker = 'x', label = 'Shifted') 

bs = 20
fig, ax = plt.subplots(1,1)
ax.hist(gns2['x2']-gns1['x2'],bins = bs)
sys.exit(74)
# %%
# Alignment of the simulated date using different methods

dis_min = 10
xy = np.array([x,y]).T 
xy_sh = np.array([x_sh,y_sh]).T 
ls_com = compare_lists(xy, xy_sh, dis_min)


loop = 0
comom_ls = []
dic_xy = {}
dic_Kx ={}
dic_xy_final = {}

deg_gns = 1
max_deg = 3
while deg_gns < max_deg:
   
    l1_xy = xy
    l2_xy = xy_sh
    comp = compare_lists(l1_xy,l2_xy,dis_min)
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























