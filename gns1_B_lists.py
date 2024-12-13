#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 16:19:45 2022

@author: amartinez
"""

# Generates the GNS1 second reduction with the Ks and H magnitudes

import numpy as np
import matplotlib.pyplot as plt
import astroalign as aa
from astropy.io.fits import getheader
from astropy.io import fits
from scipy.spatial import distance
import pandas as pd
import sys
import time
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time

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
#%%
# field_one = 7
# chip_one = 4
# field_two = 7
# chip_two = 1

# field_one = 10
# chip_one = 2
# field_two = 4
# chip_two = 3


# field_one = 16
# chip_one = 2
# field_two = 7
# chip_two = 1


field_one = 100#This field B6
chip_one = 4
field_two = 20
chip_two = 1

max_sig = 0.5
# max_sig = 2

red_pers = 'RS'
# red_pers = 'AL'
if field_one == 7 or field_one == 12 or field_one == 10 or field_one == 16:
    t1 = Time(['2015-06-07T00:00:00'],scale='utc')
elif field_one == 60:
    t1 = Time(['2016-06-13T00:00:00'],scale='utc')
elif field_one ==  100:
    t1 = Time(['2016-05-20T00:00:00'],scale='utc')
else:
    print(f'NO time detected for this field_one = {field_one}')
    sys.exit()
if field_two == 7 or field_two == 5:
    t2 = Time(['2022-05-27T00:00:00'],scale='utc')
elif field_two == 4:
    t2 = Time(['2022-04-05T00:00:00'],scale='utc')
elif field_two == 20:
    t2 = Time(['2022-07-25T00:00:00'],scale='utc')
else:
    print(f'NO time detected for this field_two = {field_two}')
    sys.exit()




# Arches and Quintuplet coordinates for plotting and check if it will be covered.
# Choose Arches or Quituplet central coordinates #!!!
# arch = SkyCoord(ra = '17h45m50.65020s', dec = '-28d49m19.51468s', equinox = 'J2000').galactic
# arch_ecu = SkyCoord(ra = '17h45m50.65020s', dec = '-28d49m19.51468s', equinox = 'J2000')
arch =  SkyCoord('17h46m15.13s', '-28d49m34.7s', frame='icrs',obstime ='J2016.0').galactic#Quintuplet
arch_ecu =  SkyCoord('17h46m15.13s', '-28d49m34.7s', frame='icrs',obstime ='J2016.0')#Quintuplet


GNS_1='/Users/amartinez/Desktop/PhD/HAWK/GNS_1/lists/%s/chip%s/'%(field_one, chip_one)
GNS_2='/Users/amartinez/Desktop/PhD/HAWK/GNS_2/lists/%s/chip%s/'%(field_two, chip_two)


with open(GNS_2 + 'Reduction_by.txt', 'w') as f:
    f.close()
if field_two == '5':
    if red_pers == 'RS':
        GNS_2='/Users/amartinez/Desktop/PhD/HAWK/GNS_2/lists/%s_RS/chip%s/'%(field_two, chip_two)
        with open(GNS_2 + 'Reduction_by.txt', 'a') as f:
            f.write('RS')
            f.close()
            
    elif red_pers == 'AL':
        GNS_2='/Users/amartinez/Desktop/PhD/HAWK/GNS_2/lists/%s_AL/chip%s/'%(field_two, chip_two)


GNS_2absolute = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2absolute_python/lists/%s/chip%s/'%(field_two, chip_two)
GNS_1absolute = '/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/lists/%s/chip%s/'%(field_one, chip_one)

np.savetxt('/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/lists/fields_and_chips.txt',
           np.array([field_one, chip_one, field_two, chip_two,t1.decimalyear[0],t2.decimalyear[0],max_sig]).reshape(1, -1), fmt = 4*'%.0f ' +3*'%.4f ')



pruebas1 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/pruebas/'
pruebas2 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2abolute_python/pruebas/'

# 0     1    2   3   4   5   6    7    8    9   
# ra1, dec1, x1, y1, f1, H1, dx1, dy1, df1, dH1 = np.loadtxt(GNS + 'stars_calibrated_H_chip1.txt', unpack = True)
# gns1_H = np.loadtxt(GNS_1 + 'stars_calibrated_H_chip%s.txt'%(chip_one))
# gns1_K = np.loadtxt(GNS_1 + 'stars_calibrated_Ks_chip%s.txt'%(chip_one))

gns1_H = Table.read(GNS_1 + 'stars_calibrated_H_chip%s.txt'%(chip_one), names = ('ra1',	'Dec1',	'x1',	'y1',	'f1',	'H1',	'dx1',	'dy1',	'df1',	'dH1'), format = 'ascii')
gns1_K = Table.read(GNS_1 + 'stars_calibrated_Ks_chip%s.txt'%(chip_one), names = ('ra1',	'Dec1',	'x1',	'y1',	'f1',	'Ks1',	'dx1',	'dy1',	'df1',	'dKs1'),format = 'ascii')

gns1H_coor = SkyCoord(ra = gns1_H['ra1'],dec = gns1_H['Dec1'], unit = 'degree', frame = 'fk5',equinox ='J2000',obstime='J2015.4')
gns1K_coor = SkyCoord(ra = gns1_K['ra1'],dec = gns1_K['Dec1'], unit = 'degree', frame = 'fk5',equinox ='J2000',obstime='J2018.4')

max_sep = 0.050 * u.arcsec
idx,d2d,d3d = gns1H_coor.match_to_catalog_sky(gns1K_coor)# ,nthneighbor=1 is for 1-to-1 match
sep_constraint = d2d < max_sep
gns1H_match = gns1_H[sep_constraint]
gns1K_match = gns1_K[idx[sep_constraint]]





gns1H_match['Ks1'] = gns1K_match['Ks1']
gns1H_match['dKs1'] = gns1K_match['dKs1']

gns1_all = gns1H_match

# unc_cut = np.where(np.sqrt(gns1_all[:,1]**2 + gns1_all[:,3]**2)<max_sig)
unc_cut = np.where((gns1_all['dx1']<max_sig) & (gns1_all['dy1']<max_sig))
gns1 = gns1_all[unc_cut]
# np.savetxt(GNS_1off +'stars_calibrated_HK_chip%s_sxy%s.txt'%(chip_one,max_sig), gns1, fmt ='%.8f',
#                  header = 'ra1 0, dec1 1, x1 2, y1 3, f1 4, H1 5, dx1 6, dy1 7, df1 8, dH1 9 ,Ks 10, dKs 11')
#Cut GNS1 on GNS2 and GNS2 on GNS1. This way we will aligning the exact same
# areas in both surveys and using the same Gaia stars, so both alignment will
# suffer the same kind of distorsions

gns1_gal = SkyCoord(ra = gns1['ra1'], dec = gns1['Dec1'], 
                    unit = 'degree', frame = 'fk5', equinox = 'J2000',
                    obstime = 'J2015.43').galactic

# %%
gns2_all = Table.read(GNS_2 + 'stars_calibrated_H_chip%s.txt'%(chip_two), names = ('ra2',	'Dec2',	'x2',	'y2',	'f2',	'H2',	'dx2',	'dy2',	'df2',	'dH2'), format = 'ascii')
unc_cut2 = np.where((gns2_all['dx2']<max_sig) & (gns2_all['dy2']<max_sig))
# gns2 = gns2_all[unc_cut2]
gns2 = gns2_all
gns2_gal = SkyCoord(ra = gns2['ra2'], dec = gns2['Dec2'], 
                    unit = 'degree', frame = 'fk5', equinox = 'J2000',
                    obstime = 'J2022.4').galactic


# %%

l2 = gns2_gal.l.wrap_at('360d')
l1 = gns1_gal.l.wrap_at('360d')
fig, ax = plt.subplots(1,1,figsize =(10,10))
ax.scatter(l1, gns1_gal.b,label = 'GNS_1 Fied %s, chip %s'%(field_one,chip_one),zorder=3)
ax.scatter(l2, gns2_gal.b,s = 1,color = 'r',label = 'GNS_2 Fied %s, chip %s'%(field_two,chip_two))

ax.invert_xaxis()
ax.legend()
ax.set_xlabel('l[deg]', fontsize = 40)
ax.set_ylabel('b [deg]', fontsize = 40)

# f_work = '/Users/amartinez/Desktop/PhD/Thesis/document/mi_tesis/tesis/Future_work/'
# plt.savefig(f_work+ 'gsn1_gns2_fields.png', bbox_inches='tight')

# %%

if arch.b.value < 0.1:
    clus_name = 'Quintuplet'
else: 
    clus_name = 'Arches'


# ax.scatter(arch.l, arch.b,s = 200, label = clus_name)


buenos1 = np.where((gns1_gal.l>min(gns2_gal.l)) & (gns1_gal.l<max(gns2_gal.l)) &
                   (gns1_gal.b>min(gns2_gal.b)) & (gns1_gal.b<max(gns2_gal.b)))

gns1 = gns1[buenos1]
gns1['ID'] = np.arange(len(gns1))


# np.savetxt(GNS_1relative +'stars_calibrated_HK_chip%s_on_gns2_f%sc%s_sxy%s.txt'%(chip_one,field_two,chip_two,max_sig), 
#            gns1, fmt ='%.8f',
                  # header = 'ra1 0, dec1 1, x1 2, y1 3, f1 4, H1 5, dx1 6, dy1 7, df1 8, dH1 9 ,Ks 10, dKs 11, ID 12')
gns1.write(GNS_1absolute + 'stars_calibrated_HK_chip%s_on_gns2_f%sc%s_sxy%s.txt'%(chip_one,field_two,chip_two,max_sig), format = 'ascii', overwrite = True)

buenos2 = np.where((gns2_gal.l>min(gns1_gal.l)) & (gns2_gal.l<max(gns1_gal.l)) &
                   (gns2_gal.b>min(gns1_gal.b)) & (gns2_gal.b<max(gns1_gal.b)))

gns2 = gns2[buenos2]
gns2['ID'] = np.arange(len(gns2))
# %%
fig, ax = plt.subplots(1,1,figsize =(10,10))

ax.scatter(gns2['ra2'], gns2['Dec2'],label = 'GNS_2 Fied %s, chip %s'%(field_two,chip_two))
ax.scatter(gns1['ra1'], gns1['Dec1'],label = 'GNS_1 Fied %s, chip %s'%(field_one,chip_one))
# ax.scatter(arch_ecu.ra, arch_ecu.dec,s = 200, label = clus_name)
ax.legend()
ax.set_xlabel('Ra(deg)')
ax.set_ylabel('Dec(deg)')



# np.savetxt(GNS_2relative +'stars_calibrated_H_chip%s_on_gns1_f%sc%s_sxy%s.txt'%(chip_two,field_one,chip_one,max_sig)
#            ,gns2, fmt ='%.8f',
#             header = 'ra1 0, dec1 1, x1 2, y1 3, f1 4, H1 5, dx1 6, dy1 7, df1 8, dH1 9 , ID 10')
gns2.write(GNS_2absolute +'stars_calibrated_H_chip%s_on_gns1_f%sc%s_sxy%s.txt'%(chip_two,field_one,chip_one,max_sig), format = 'ascii', overwrite = True)







