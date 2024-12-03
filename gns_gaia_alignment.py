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
# %%
# =============================================================================
# Gaia data:https://gea.esac.esa.int/archive/
# =============================================================================
# Load coordinates for the gns1 pointing
# field_one = 10
# chip_one = 2
# field_two = 4
# chip_two = 3
# field_one = 7
# chip_one = 4
# field_two = 7
# chip_two = 1
pruebas2 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2absolute_python/pruebas/'
pruebas1 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/pruebas/'
gaia_list1 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/Gaia_lists/'
gaia_list2 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2absolute_python/Gaia_lists/'
tmp1 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/tmp/'
tmp2 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2absolute_python/tmp/'
Hos_cat = '/Users/amartinez/Desktop/PhD/Arches_and_Quintuplet_Hosek/'
gaia_hos = fits.open(Hos_cat + 'quintuplet_HST_gaia_matched_pm_EOM.fits')
pm_folder_off ='/Users/amartinez/Desktop/PhD/HAWK/pm_gns1_gns2_absolute_python/'


field_one, chip_one, field_two, chip_two,t1,t2,max_sig = np.loadtxt('/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/lists/fields_and_chips.txt', 
                                                       unpack=True)
field_one = field_one.astype(int)
chip_one = chip_one.astype(int)
field_two = field_two.astype(int)
chip_two = chip_two.astype(int)


GNS_1='/Users/amartinez/Desktop/PhD/HAWK/GNS_1/lists/%s/chip%s/'%(field_one, chip_one)

GNS_1off='/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_python/lists/%s/chip%s/'%(field_one, chip_one)
GNS_2off='/Users/amartinez/Desktop/PhD/HAWK/GNS_2absolute_python/lists/%s/chip%s/'%(field_two, chip_two)
# gns1_im = GNS_1 +  'field%s_H_gns1.fits'%(field_one)
gns1_im = GNS_1 +   'calibrated_stars_%s.fits'%(chip_one)

# delete previous catalogs
pm_delete_id = pm_folder_off + 'ID_pm_GaiaRF_ep1_f%sc%s_ep2_f%sc%s**.txt'%(field_one, chip_one, field_two, chip_two)
pm_delete = pm_folder_off + 'pm_GaiaRF_ep1_f%sc%s_ep2_f%sc%s**.txt'%(field_one, chip_one, field_two, chip_two)

max_sig = 0.5


for file in glob.glob(pm_delete):
    try:
        os.remove(file)
        print(f"Deleted: {file}")
    except:
        print('yomamma')
for file in glob.glob(pm_delete_id):
    try:
        os.remove(file)
        print(f"Deleted: {file}")
    except:
        print('yomamma')



max_sig = 0.5#TODO
# max_sig = 1#TODO



bad1 = []  
bad2 = []  
#


# bad_all = bad1 + bad2
# bad_all = np.unique(bad1 + bad2)
# bad1 = list(bad_all)
# bad2 = list(bad_all)


# Save this lists!
# They make the alignment for f10_c2, f4_c3 (Quintuplet) works (carefull with star 42)
# bad1 = [14, 28, 54, 19, 37, 47,58,56,55,23,9,11,7,40,0,33,12,49,42]  
# bad2 = [14, 28, 54, 19, 37, 47,58,56,55,23,9,11,7,40,0,33,12,49,]  

np.savetxt(tmp1 + 'bad1_f%sc%s.txt'%(field_one,chip_one),np.array(bad1).T, fmt='%.i')
np.savetxt(tmp2 + 'bad2_f%sc%s.txt'%(field_two,chip_two),np.array(bad2).T, fmt='%.i')

if len(bad1) >0:
    clip_bad1 = 'yes'#TODO
else:
    clip_bad1 = 'no'#TODO
# ra1 0, dec1 1, x1 2, y1 3, f1 4, H1 5, dx1 6, dy1 7, df1 8, dH1 9 ,Ks 10, dKs 11
# gns1 = np.loadtxt(GNS_1off +'stars_calibrated_HK_chip%s_on_gns2_f%sc%s_sxy%s.txt'%(chip_one,field_two,chip_two,max_sig))
gns1 = Table.read(GNS_1off +'stars_calibrated_HK_chip%s_on_gns2_f%sc%s_sxy%s.txt'%(chip_one,field_two,chip_two,max_sig), format = 'ascii')


p_mas = 0.5*0.106*1000#Trnasform pixeles into mas


# gns1['x1'] = gns1['x1']
gns1['x1'] = gns1['x1']*-1
gns1['x1'] = gns1['x1']*p_mas
gns1['y1'] = gns1['y1']*p_mas
# gns1[:,2] = gns1[:,2]
# np.savetxt(GNS_1 +'stars_calibrated_HK_chip%s_sxy%s.txt'%(chip_one,max_sig),gns1, fmt ='%.8f',header='x, dx, y, dy, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK')
gns1_coor =  SkyCoord(ra=gns1['ra1'], dec=gns1['Dec1'],unit = 'degree' ,frame = 'fk5',equinox ='J2000',obstime='J2015.5')
gns1_coor =gns1_coor.transform_to('icrs')


# Here we are going get the informations in the header of GNS1
f =fits.open(gns1_im)
w = WCS(f[0].header)    
header = fits.getheader(gns1_im)
# x_cent, y_cent =header['NAXIS1']/2,header['NAXIS2']/2
x_cent, y_cent = np.mean(gns1['ra1']), np.mean(gns1['Dec1'])
radec_ = SkyCoord(ra = x_cent, dec = y_cent, unit = 'degree', frame = 'fk5', equinox = 'J2000', obstime ='J2015.5')
with open(pruebas1 + 'ZP_centroid_gns1_f%sc%s.reg'%(field_one, chip_one),'w') as file:
    file.write('# Region file format: DS9 version 4.1'+"\n"+'global color=blue dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+"\n"+'fk5'+'\n')
with open(pruebas1 + 'ZP_centroid_gns1_f%sc%s.reg'%(field_one, chip_one),'a') as file:
    file.write('circle(%s,%s,0.5")'%(x_cent, y_cent))




if field_one == 7 or field_one == 12 or field_one == 10 or field_one == 16:
    t_gns1 = Time(['2015-06-07T00:00:00','2016-01-01T00:00:00'],scale='utc')

if field_two == 7 or field_two == 5:
    t_gns2 = Time(['2022-05-27T00:00:00','2016-01-01T00:00:00'],scale='utc')
if field_two == 4:
    t_gns2 = Time(['2022-04-05T00:00:00','2016-01-01T00:00:00'],scale='utc')


del_t = t_gns2[0]-t_gns1[0]
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source" # Select early Data Release 3
# Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source" # Select early Data Release 3
Gaia.ROW_LIMIT = -1  # Ensure the default row limit.
coord = SkyCoord(ra=radec_.ra, dec=radec_.dec, unit='degree', frame='icrs',equinox ='J2000',obstime='J2016.0')


ancho2 = header['NAXIS1']
alto2 = header['NAXIS2']
width2 = ancho2*0.053*u.arcsec
height2 = alto2*0.053*u.arcsec


# =============================================================================
# j = Gaia.cone_search_async(coord, np.sqrt((width2)**2+height2**2)/2)
# # j = Gaia.cone_search_async(coord, 1000*u.arcsec)
# gaia_ = j.get_results()
# gaia_.write(pruebas1 + 'gaia1_f%sc%s.txt'%(field_one,chip_one), format = 'ascii', overwrite = True)
# =============================================================================

gaia_ = Table.read(pruebas1 + 'gaia1_f%sc%s.txt'%(field_one,chip_one), format = 'ascii')

# gaia_com_ind = np.isin(gaia_['source_id'],gaia_hos[1].data['source_id_2'])

e_pm = 0.3
# gaia_com = gaia_[gaia_com_ind]
# If you get gaia_ from the Gaia website use "duplicated_source= False"
gaia_com = gaia_
# gaia_good = filter_gaia_data(
#     gaia_table=gaia_,
#     astrometric_params_solved=31,
#     duplicated_source= False,
#     parallax_over_error_min=-10,
#     astrometric_excess_noise_sig_max=2,
#     phot_g_mean_mag_min= None,
#     phot_g_mean_mag_max= 13,
#     pm_min=0,
#     pmra_error_max=e_pm,
#     pmdec_error_max=e_pm
#     )
# If you get gaia_ from a stored gaia Table, use "duplicated_source= 'False'"
gaia_good = filter_gaia_data(
    gaia_table=gaia_,
    astrometric_params_solved=31,
    duplicated_source= 'False',
    parallax_over_error_min=-10,
    astrometric_excess_noise_sig_max=2,
    phot_g_mean_mag_min= None,
    phot_g_mean_mag_max= 13,
    pm_min=0,
    pmra_error_max=e_pm,
    pmdec_error_max=e_pm
    )

gaia_coord =  SkyCoord(ra=gaia_good['ra'], dec=gaia_good['dec'],unit = 'degree',frame = 'icrs',obstime='J2016.0')


np.savetxt(pruebas1 + 'gaia.txt', np.array([gaia_good['ra'],gaia_good['dec']]).T)
# %
delta_t1 = t_gns1[0] - t_gns1[1]
delta_t2 = t_gns2[0] - t_gns2[1]
# delta_t1 =  t_gns1[0] - t_gns1[0]
# delta_t2 =  t_gns2[0] - t_gns2[0]

GaiaCoord = SkyCoord(ra=gaia_good['ra'],
                   dec=gaia_good['dec'],
                   unit = 'degree',
                   frame='icrs',
                     equinox = 'J2000',
                    obstime='J2016.0')
if field_one == 7 and chip_one ==4:
    # center_arc = SkyCoord(ra = '17h45m50.65020s', dec = '-28d49m19.51468s', equinox = 'J2000') if choosen_cluster =='Arches' else SkyCoord('17h46m14.68579s', '-28d49m38.99169s', equinox = 'J2000')#Quintuplet
    radec_ = SkyCoord(ra = '17h45m50.65020s', dec = '-28d49m19.51468s', equinox = 'J2000')
    print('Choosing Arches center as reference point')
if field_one == 10 and chip_one ==2:
    # center_arc = SkyCoord(ra = '17h45m50.65020s', dec = '-28d49m19.51468s', equinox = 'J2000') if choosen_cluster =='Arches' else SkyCoord('17h46m14.68579s', '-28d49m38.99169s', equinox = 'J2000')#Quintuplet
    radec_ = SkyCoord('17h46m14.68579s', '-28d49m38.99169s', equinox = 'J2000')#Quintuplet
    print('Choosing Arches center as reference point')

# Asigns Offsets coordinates to Gaia stars, moves then  and return the
# corresponding ra dec coordenates using spherical_offsets_by method
radec_ = radec_.transform_to('icrs')
ragai1_off,decgai1_off = radec_.spherical_offsets_to(GaiaCoord.frame)
ragai1_off = (ragai1_off.to(u.mas)).value + (np.array(gaia_good['pmra'])*delta_t1.to(u.yr)).value
decgai1_off = (decgai1_off.to(u.mas)).value + (np.array(gaia_good['pmdec'])*delta_t1.to(u.yr)).value
GaiaGNSCoord1 = radec_.spherical_offsets_by(ragai1_off*u.mas, decgai1_off*u.mas)
# Could be this transformation of Gaia stars (from the sky to the plane, and back) a source of error??

# %

# Propagation of error in gaia
dx_des = np.sqrt(gaia_good['ra_error']**2 + (delta_t1.to(u.yr).value*gaia_good['pmra_error'])**2)
dy_des = np.sqrt(gaia_good['dec_error']**2 + (delta_t1.to(u.yr).value*gaia_good['pmra_error'])**2)
dxy_des = np.sqrt(dx_des**2 + dy_des**2) 
# %
# Here we are going to cut the Gaia stars over the area of GNS1
gaia_coord_gal = gaia_coord.galactic
gns1_coord_gal = gns1_coor.galactic

lg = gaia_coord_gal.l.wrap_at('360d')
l1 = gns1_coord_gal.l.wrap_at('360d')


buenos1 = np.where((gaia_coord_gal.l>min(gns1_coord_gal.l)) & (gaia_coord_gal.l<max(gns1_coord_gal.l)) &
                   (gaia_coord_gal.b>min(gns1_coord_gal.b)) & (gaia_coord_gal.b<max(gns1_coord_gal.b)))

fig, ax = plt.subplots(1,1,figsize =(6,6))
# ax.scatter(lg, gaia_coord_gal.b, label ='Gaia stars')
ax.scatter(l1, gns1_coord_gal.b, label = 'GNS1 (f%s c%s)'%(field_one, chip_one) )  
ax.scatter(lg[buenos1], gaia_coord_gal[buenos1].b, label ='Gaia over GNS1')
ax.invert_xaxis()
ax.set_xlabel('l[deg]', fontsize = 20)
ax.set_ylabel('b [deg]',fontsize = 20)
ax.axis('equal')
ax.legend()



# Each gaia stars get its own ID
ga1_id = np.arange(len(gaia_good[buenos1]))
gaia_np1 =np.array([GaiaGNSCoord1[buenos1].ra.value, GaiaGNSCoord1[buenos1].dec.value,
                     dx_des[buenos1], dy_des[buenos1],
                     ragai1_off[buenos1], decgai1_off[buenos1],
                     np.array(gaia_good['pmra'][buenos1]), np.array(gaia_good['pmdec'][buenos1]),
                     gaia_good['pmra_error'][buenos1].value, gaia_good['pmdec_error'][buenos1].value,
                     gaia_good['phot_g_mean_mag'][buenos1].value,
                     ga1_id]).T
# Saves two differnte lists: one with All gaia stars and the other after the 
# 3sigma clipping done when comparing with Gaia.
gaia1 = Table(gaia_np1, names = ('ra',	'dec',	'dra',	'ddec',	'x',	'y',	'pmra',	'pmdec',	'dpmra',	'dpmdec','phot_g_mean_mag','gaia1_id'))

gaia1.write(GNS_1off + 'ALL_gaia_refstars_on_gns1_f%sc%s.txt'%(field_one,chip_one), format = 'ascii', overwrite = True) 

with open(gaia_list1+ 'gaia_gns1_f%sc%s_on_gns2_f%sc%s_sxy%s.reg'%(field_one,chip_one,field_two,chip_two,max_sig), 'w') as f:
    f.write('# Region file format: DS9 version 4.1'+"\n"+'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+"\n"+'fk5'+'\n')
    f.close
for gs in range(len(ga1_id)):
    with open(gaia_list1+ 'gaia_gns1_f%sc%s_on_gns2_f%sc%s_sxy%s.reg'%(field_one,chip_one,field_two,chip_two,max_sig), 'a') as f:
        f.write('circle(%s,%s,0.5") \n point(%s,%s) # point=cross \n # text(%s,%s) font="helvetica 24 normal roman" text={%.0f} \n'%(gaia1['ra'][gs], gaia1['dec'][gs],
                                                                        gaia1['ra'][gs], gaia1['dec'][gs],
                                                                        gaia1['ra'][gs]+0.00013, gaia1['dec'][gs]+0.00013,
                                                                        gs))

if clip_bad1 == 'yes':#TODO
    del_1 = np.isin(gaia1['gaia1_id'], bad1)
    gaia1 = gaia1[~del_1]


    # gaia_np1 = np.delete(gaia_np1,bad1, axis =0)
# np.savetxt(GNS_1off + 'gaia_refstars_on_gns1_f%sc%s.txt'%(field_one,chip_one), gaia_np1, fmt ='%.8f', 
#            header = 'ra, dec, dra(mas), ddec(mas), x, y, pmra(mas/yr), pmdec, dpmra(mas/yr), dpmdec,dradec')
gaia1.write(GNS_1off + 'gaia_refstars_on_gns1_f%sc%s.txt'%(field_one,chip_one),format = 'ascii',overwrite = True)

ID_gns1 = np.arange(len(gns1[buenos1]))

# with open(pruebas1+ 'gns1_around_gaiaf%sc%s_sxy%s.reg'%(field_one,chip_one,max_sig), 'w') as f:
#     f.write('# Region file format: DS9 version 4.1'+"\n"+'global color=blue dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+"\n"+'fk5'+'\n')
#     f.close

# for gst in range(len(gaia_np1)):
#     m_point = SkyCoord(ra =[gaia_np1[gst][0]], dec  =[gaia_np1[gst][1]],unit ='degree', frame ='icrs', equinox = 'J2000', obstime = 'J2015.43')
#     idx, d2d, d3d = match_coordinates_sky(m_point, gns1_coor)
#     # idxc, group_md, d2d,d3d =  ap_coord.search_around_sky( m_point ,gns1_coor, 0.6*u.arcsec)
#     with open(pruebas1+ 'gns1_around_gaiaf%sc%s_sxy%s.reg'%(field_one,chip_one,max_sig), 'a') as f:
#         f.write('point(%s,%s) # point = cross \n # text(%s,%s) text={%.0f} \n'%(float(gns1_coor.ra.value[idx]) ,float(gns1_coor.dec.value[idx]),
#                                                                                 float(gns1_coor.ra.value[idx]) ,float(gns1_coor.dec.value[idx]),
#                                                                                int(ID_gns1[idx])))
#         f.close 

# Asigns Offsets coordinates to Gaia stars, moves then (epoch 2) and return the
# corresponding ra dec coordenates using spherical_offsets_by method
radec_ = radec_.transform_to('icrs')
ragai2_off,decgai2_off = radec_.spherical_offsets_to(GaiaCoord.frame)
ragai2_off = (ragai2_off.to(u.mas)).value + (np.array(gaia_good['pmra'])*delta_t2.to(u.yr)).value
decgai2_off = (decgai2_off.to(u.mas)).value + (np.array(gaia_good['pmdec'])*delta_t2.to(u.yr)).value
GaiaGNSCoord2 = radec_.spherical_offsets_by(ragai2_off*u.mas, decgai2_off*u.mas)



dx_des = np.sqrt(gaia_good['ra_error']**2 + (delta_t2.to(u.yr).value*gaia_good['pmra_error'])**2)
dy_des = np.sqrt(gaia_good['dec_error']**2 + (delta_t2.to(u.yr).value*gaia_good['pmra_error'])**2)
dxy_des = np.sqrt(dx_des**2 + dy_des**2) 
ga2_id = np.arange(len(gaia_good))
# %
# Here we are going to cut the Gaia stars over the are of GNS2

gns2 = Table.read(GNS_2off + 'stars_calibrated_H_chip%s_on_gns1_f%sc%s_sxy%s.txt'%(chip_two,field_one,chip_one,max_sig), format = 'ascii')
# gns2['x2'] = gns2['x2']#TODO
gns2['x2'] = gns2['x2']*-1#TODO
gns2['x2'] = gns2['x2']*p_mas
gns2['y2'] = gns2['y2']*p_mas#TODO
# gns2[:,3] = gns2[:,3]*-1#TODO
gns2_coor = SkyCoord(ra= gns2['ra2']*u.degree, dec=gns2['Dec2']*u.degree, frame = 'fk5', equinox = 'J2000',obstime='J2022.4')
gns2_coor = gns2_coor.transform_to('icrs')

gaia_coord_gal = gaia_coord.galactic
gns2_coord_gal = gns2_coor.galactic

lg = gaia_coord_gal.l.wrap_at('360d')
l2 = gns2_coord_gal.l.wrap_at('360d')


buenos2 = np.where((gaia_coord_gal.l>min(gns2_coord_gal.l)) & (gaia_coord_gal.l<max(gns2_coord_gal.l)) &
                   (gaia_coord_gal.b>min(gns2_coord_gal.b)) & (gaia_coord_gal.b<max(gns2_coord_gal.b)))

fig, ax = plt.subplots(1,1,figsize =(10,10))
ax.scatter(lg, gaia_coord_gal.b, label = 'Gaia stars')
ax.scatter(l2, gns2_coord_gal.b, label = 'GNS2 (f%s, c%s'%(field_two, chip_two))
ax.scatter(lg[buenos2], gaia_coord_gal[buenos2].b, label ='Gaia over GNS2')
ax.set_xlabel('l(deg)')
ax.set_ylabel('b(deg)')
ax.invert_xaxis()
ax.axis('equal')
ax.legend()


ga2_id = np.arange(len(gaia_good[buenos2]))

gaia_np2 = np.array([GaiaGNSCoord2[buenos2].ra.value,GaiaGNSCoord2[buenos2].dec.value,
                   dx_des[buenos2],dy_des[buenos2],
                   ragai2_off[buenos2], decgai2_off[buenos2],
                   np.array(gaia_good['pmra'][buenos2]), np.array(gaia_good['pmdec'][buenos2]),
                   gaia_good['pmra_error'][buenos2].value,gaia_good['pmdec_error'][buenos2].value,
                   gaia_good['phot_g_mean_mag'][buenos2].value,
                   ga2_id]).T
gaia2 = Table(gaia_np2, names = ('ra',	'dec',	'dra',	'ddec',	'x',	'y',	'pmra',	'pmdec',	'dpmra',	'dpmdec','phot_g_mean_mag','gaia2_id'))

gaia2.write(GNS_2off + 'ALL_gaia_refstars_on_gns2_f%sc%s_gns1_f%sc%s.txt'%(field_two,chip_two,field_one,chip_one), format = 'ascii', overwrite = True) 


with open(gaia_list2+ 'gaia_gns2_f%sc%s_on_gns_f%sc%s_sxy%s.reg'%(field_two,chip_two,field_one,chip_one,max_sig), 'w') as f:
    f.write('# Region file format: DS9 version 4.1'+"\n"+'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+"\n"+'fk5'+'\n')
    f.close
for gs in range(len(gaia_np2)):
    with open(gaia_list2+ 'gaia_gns2_f%sc%s_on_gns_f%sc%s_sxy%s.reg'%(field_two,chip_two,field_one,chip_one,max_sig), 'a') as f:
        f.write('circle(%s,%s,0.5") \n point(%s,%s) # point=cross \n # text(%s,%s) font="helvetica 24 normal roman" text={%.0f} \n'%(gaia_np2[gs][0], gaia_np2[gs][1],
                                                                        gaia_np2[gs][0], gaia_np2[gs][1],
                                                                        gaia_np2[gs][0]+0.00013, gaia_np2[gs][1]+0.00013,
                                                                        gs))

if len(bad2) >0:
    clip_bad2 = 'yes'#TODO
else:
    clip_bad2 = 'no'#TODO                 
if clip_bad2 == 'yes':#TODO
    del_2 = np.isin(gaia2['gaia2_id'], bad2)
    gaia2 = gaia2[~del_2]
 
gaia2.write(GNS_2off + 'gaia_refstars_on_gns2_f%sc%s_gns1_f%sc%s.txt'%(field_two,chip_two,field_one,chip_one), format = 'ascii', overwrite = True)


# %%
# We select GNS1 foregroud stars and astroaling their pixels coordenates with 
# the Gaia stars offsets

fg = (gns1['H1']-gns1['Ks1'])<1.3    
gns1_fg = gns1[fg]
# gns1_fg = gns1


max_sep = 0.08*u.arcsec#!!!
gaia1_coord = SkyCoord(ra = gaia1['ra'], dec = gaia1['dec'], unit = 'degree',frame = 'icrs',obstime='J2016.0' )
gns1_coor_fg= SkyCoord(ra= gns1_fg['ra1']*u.degree, dec=gns1_fg['Dec1']*u.degree, frame = 'fk5',obstime='J2015.43')
gns1_coor_match = gns1_coor_fg.transform_to('icrs')
idx,d2d,d3d = gaia1_coord.match_to_catalog_sky(gns1_coor_fg,nthneighbor=1)# ,nthneighbor=1 is for 1-to-1 match
sep_constraint = d2d < max_sep
gaia1_match = gaia1[sep_constraint]
gns1_match = gns1_fg[idx[sep_constraint]]


xy_gns1 = np.array([gns1_match['x1'],gns1_match['y1']]).T
xy_gaia1 = np.array([gaia1_match['x'],gaia1_match['y']]).T

# p = ski.transform.estimate_transform('polynomial',
#                                               xy_gns1, 
#                                               xy_gaia1, order = 1)
# p = ski.transform.estimate_transform('affine',
#                                               xy_gns1, 
#                                               xy_gaia1)
p = ski.transform.estimate_transform('similarity',
                                              xy_gns1, 
                                              xy_gaia1)

fig, ax = plt.subplots(1,1)
ax.scatter(gns1_match['x1'],gns1_match['y1'], marker = '*', label = 'GNS1')
ax.scatter(gaia1_match['x'],gaia1_match['y'], marker = '*', label = 'Gaia')
ax.legend()
gns1_xy = np.array([gns1['x1'],gns1['y1']]).T
gns1_xyt = p(gns1_xy)
gns1['x1'] = gns1_xyt[:,0] 
gns1['y1'] = gns1_xyt[:,1] 

# sys.exit(491)

 # %%
fig, ax = plt.subplots(1,1)
ax.set_title(f'GNS1 f{field_one} c{chip_one}')
ax.scatter(gns1['x1'],gns1['y1'])
ax.scatter(gaia1['x'],gaia1['y'])

# pix_scale = 0.1064*0.53
# d_m = max_sep.value*pix_scale
d_m = 30 #!!! pixeles are in mas
s_ls = compare_lists(np.array([gns1['x1'],gns1['y1']]).T, np.array([gaia1['x'],gaia1['y']]).T, d_m)
ax.scatter(s_ls['l1_x'],s_ls['l1_y'],s=20, marker = 'x', label = f'Matching = %s\nbefore = {len(xy_gns1)}'%(len(s_ls['l1_x'])))
ax.legend()
print(30*'-'+f'\nCommon GNS1 and Gaia after initial transformation = f{len(s_ls)}')
# %%
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

gns1 = alignator(1, gns1, gaia1, s_ls, 20, 4, 1)
# gns2 = alignator(2, gns2, gaia2, s_ls, 20, 4, 1)


# # %%
# loop =1
# deg = 1

# comom_ls =[]
# comom_ls.append(len(s_ls))
# dic_xy = {} 
# dic_xy['trans_0'] = np.array([gns1['x1'],gns1['y1']]).T
# # for loop in range(1,10):
# max_deg = 4#!!!
# sig_cl = 3#!!!
# # %%
# # Alignemnt of GNS1
# while deg < max_deg:
    
#     Kx,Ky=pw.polywarp(s_ls['l2_x'],s_ls['l2_y'],s_ls['l1_x'],s_ls['l1_y'],degree=deg)
    
#     xi=np.zeros(len(gns1))
#     yi=np.zeros(len(gns1))
    
#     for k in range(deg+1):
#                 for m in range(deg+1):
#                     xi=xi+Kx[k,m]*gns1['x1']**k*gns1['y1']**m
#                     yi=yi+Ky[k,m]*gns1['x1']**k*gns1['y1']**m
#     dic_xy[f'trans_{loop}'] = np.array([xi,yi]).T
    
#     # print(Kx[0][0])
#     gns1['x1'] = xi
#     gns1['y1'] = yi
    
#     l1_xy = np.array([gns1['x1'],gns1['y1']]).T
#     g1_xy = np.array([gaia1['x'], gaia1['y']]).T
#     # l2_xy = np.array([gns2['x2'],gns2['y2']]).T
    
    
#     s_ls = compare_lists(l1_xy,g1_xy,d_m)
#     print(f'\nCommon GNS1 and Gaia after loop{loop} = {len(s_ls)}')
    
#     diff_x =s_ls['l2_x'] - s_ls['l1_x']  
#     diff_y =s_ls['l2_y'] - s_ls['l1_y']  
    
#     mask_x, lx_lim,hx_lim = sigma_clip(diff_x, sigma=sig_cl, masked = True, return_bounds= True)
#     mask_y, ly_lim,hy_lim = sigma_clip(diff_y, sigma=sig_cl, masked = True, return_bounds= True)

    
#     comom_ls.append(len(s_ls))
#     if comom_ls[-1] <= comom_ls[-2] :
#         if len(comom_ls)>2:
            
#             gns1['x1'] = dic_xy[f'trans_{loop-1}'][:,0]
#             gns1['y1'] = dic_xy[f'trans_{loop-1}'][:,1]
#             deg +=1
#             print(30*'-'+f'\nBreaking after loop = {loop-1}\nStarting alignment degree = {deg}\n' + 30*'-')
#             continue
#         else:
#             print(f'Polynomial of degree {deg} does not work')
#             gns1['x1'] = dic_xy['trans_0'][:,0]
#             gns1['y1'] = dic_xy['trans_0'][:,1]
#             break
       
    
    
#     fig, (ax,ax1)  = plt.subplots(1,2)
#     ax1.set_title('Degree = %s'%(deg))
#     ax.set_title(f'Loop = {loop}. Matching stars = {len(s_ls)}')
#     ax.hist(diff_x, histtype = 'step', label = '$\overline{x} = %.2f$\n$\sigma x$ =%.2f'%(np.mean(diff_x),np.std(diff_x)))
#     ax.axvline(lx_lim, color = 'red', ls = 'dashed', label = '$\pm$%s$\sigma$'%(sig_cl))
#     ax.axvline(hx_lim, color = 'red', ls = 'dashed')
#     ax1.axvline(ly_lim, color = 'red', ls = 'dashed')
#     ax1.axvline(hy_lim, color = 'red', ls = 'dashed')
#     ax.legend()
#     ax1.hist(diff_y, histtype = 'step', label = '$\overline{y} = %.2f$\n$\sigma x$ =%.2f'%(np.mean(diff_y),np.std(diff_y)))
#     ax1.legend()
    
#     loop +=1     

     
# %%
# Alignment of GNS
    


