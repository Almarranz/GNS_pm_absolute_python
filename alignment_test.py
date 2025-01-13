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
from gns_alignator import gns_alignator
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
std_shift = 5  # Standard deviation of the Gaussian
shifts = np.random.normal(mean_shift, std_shift, num_points)

# Angles for displacement direction
angles = np.random.uniform(0, 2 * np.pi, num_points)

# Calculate the displacement components
dx = shifts * np.cos(angles)
dy = shifts * np.sin(angles)

# New distribution with Gaussian shifts
x_sh = gns2['x2'] + dx
y_sh = gns2['y2'] + dy


gns1['x1'] = x_sh
gns1['y1'] = y_sh

gns1['H1'] = gns1['H2']
gns1['dx1'] = gns1['dx2']
gns1['dy1'] = gns1['dy2']


fig, ax = plt.subplots(1,1)
ax.scatter(gns2['x2'][::100],gns2['y2'][::100] , label = 'Original') 
ax.scatter(gns1['x1'][::100],gns1['y1'][::100] ,marker = 'x', label = 'Shifted') 

l1 = np.array([gns2['x2'],gns2['y2']]).T
l2 = np.array([gns1['x1'],gns1['y1']]).T
ls_com = compare_lists(l1, l2, 20)

diff_mag = gns2['H2'][ls_com['ind_1']]-gns1['H1'][ls_com['ind_2']]

bs = 'auto'
fig, ax = plt.subplots(1,1)
ax.hist(diff_mag, bins = 100, label = '$\overline{\Delta H}$ = %.3f\n$\sigma$ = %.2f'%(np.mean(diff_mag),np.std(diff_mag))) 
ax.set_xlabel('$\Delta H$')
ax.legend()

fig, (ax,ax1) = plt.subplots(1,2)
ax.hist(gns2['x2']-gns1['x1'],bins = bs)
ax1.hist(gns2['y2']-gns1['y1'],bins = bs)
ax.set_xlabel('$\Delta x$')
ax1.set_xlabel('$\Delta y$')
# sys.exit(77)
# %%
# Alignment of the simulated date using different methods
min_dis = 20
max_deg = 3
align_by = 'Polywarp'
# align_by = '2DPoly'
gns2_c = gns2[ls_com['ind_1']]
gns1_c = gns1[ls_com['ind_2']]
gns1 = gns_alignator(mode = 'one_to_two',  gns_cons = gns2, gns_var = gns1,
                      align_by = align_by,f_mode = 'W',
                      d_m = min_dis, max_deg_gns = max_deg, sig_clip = 0.1,
                      plot = 'yes')




