#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 17:08:14 2024

@author: federicagucci
"""
import matplotlib.pyplot as plt
import xarray as xr
import seaborn as sns

plt.rcParams.update({'font.size': 14})

#snohats 
ds_1c_sh = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/Submission/Code&Data/SnoHATS_stable_1c.nc')
ds_2c_sh = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/Submission/Code&Data/SnoHATS_stable_2c.nc')
ds_3c_sh = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/Submission/Code&Data/SnoHATS_stable_3c.nc')

#metcrax stable
ds_1c = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/DataSamu/Metcrax_anisotropy_samu/MetCrax_stable_1c.nc')
ds_2c = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/DataSamu/Metcrax_anisotropy_samu/MetCrax_stable_2c.nc')
ds_3c = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/DataSamu/Metcrax_anisotropy_samu/MetCrax_stable_3c.nc')

#metcrax unstable
ds_1c_us = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/DataSamu/Metcrax_anisotropy_samu/MetCrax_unstable_1c.nc')
ds_2c_us = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/DataSamu/Metcrax_anisotropy_samu/MetCrax_unstable_2c.nc')
ds_3c_us = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/DataSamu/Metcrax_anisotropy_samu/MetCrax_unstable_3c.nc')


#################################################################################
## SnoHATS: Plot Pdfs of uu,vv,ww and uw,vw,uv in a panel plot
#           (each row corresponds to 1c,2c,3c)
################################################################################

fig, ax = plt.subplots(3,2, figsize =(16,6), constrained_layout = True, sharex = 'col')

for a, ds, label in zip(ax, [ds_1c_sh,ds_2c_sh,ds_3c_sh], ['1C', '2C', '3C']):
    ds = ds.stack(index=('time','sonic'))
    
    #densities
    for var in ['uu','vv', 'ww']:
        sns.kdeplot(x=ds[var]/ds.tke, ax=a[0], label = var, fill = True, linewidth = 2.5, alpha = 0.3)
    
    for var in ['uw','vw', 'uv']:
        sns.kdeplot(x=ds[var]/ds.tke, ax=a[1], label = var, fill = True, linewidth = 2.5, alpha = 0.3)

 #labels and visuals
    a[0].legend(fontsize = 15)
    a[0].set_xlim(-0.2,2.2)
    a[0].set_ylabel('Pdf', fontsize = 15, rotation = 90, labelpad = 0)
    a[0].yaxis.set_label_coords(-0.05, 0.5)
    a[0].annotate(label, xy=(-0.12, 1), xycoords='axes fraction', fontsize=25,
                    horizontalalignment='left', verticalalignment='top')
    a[0].grid(True)
 
    a[1].legend(fontsize = 15)
    a[1].set_xlim(-1.2,1.2)
    a[1].axvline(x = 0, color = 'black',lw=2,linestyle=":")
    a[1].set_ylabel('', fontsize = 15, rotation = 90, labelpad = 0)
    a[1].grid(True)

   
ax[0,0].set_title('Variances',pad=20,fontsize=25)
ax[0,1].set_title('Covariances',pad=20,fontsize=25)
ax[0,0].annotate('a)', xy=(0, 1.5), xycoords='axes fraction', fontsize=30,
                horizontalalignment='left', verticalalignment='top')
ax[-1,0].set_xlabel(r'$u_i u_i/TKE$', fontsize=20, labelpad = 10)
ax[-1,1].set_xlabel(r'$u_i u_j/TKE$', fontsize=20, labelpad = 10)

#################################################################################
## METCRAX STABLE: Plot Pdfs of uu,vv,ww and uw,vw,uv in a panel plot
#                  (each row corresponds to 1c,2c,3c)
################################################################################
fig, ax = plt.subplots(3,2, figsize =(16,6), constrained_layout = True, sharex = 'col')

for a, ds, label in zip(ax, [ds_1c,ds_2c,ds_3c], ['1C', '2C', '3C']):
    ds = ds.stack(index=('time','heights'))
    
    #densities
    for var in ['uu','vv', 'ww']:
        sns.kdeplot(x=ds[var]/ds.tke, ax=a[0], label = var, fill = True, linewidth = 2.5, alpha = 0.3)
    
    for var in ['uw','vw', 'uv']:
        sns.kdeplot(x=ds[var]/ds.tke, ax=a[1], label = var, fill = True, linewidth = 2.5, alpha = 0.3)

 #labels and visuals
    a[0].legend(fontsize = 15)
    a[0].set_xlim(-0.2,2.2)
    a[0].set_ylim(0,6)
    a[0].set_ylabel('Pdf', fontsize = 15, rotation = 90, labelpad = 0)
    a[0].yaxis.set_label_coords(-0.05, 0.5)
    a[0].annotate(label, xy=(-0.12, 1), xycoords='axes fraction', fontsize=25,
                    horizontalalignment='left', verticalalignment='top')
    a[0].grid(True)
 
    a[1].legend(fontsize = 15)
    a[1].set_xlim(-1.2,1.2)
    a[1].axvline(x = 0, color = 'black',lw=2,linestyle=":")
    a[1].set_ylabel('', fontsize = 15, rotation = 90, labelpad = 0)
    a[1].grid(True)
    
   
ax[0,0].set_title('Variances',pad=20,fontsize=25)
ax[0,1].set_title('Covariances',pad=20,fontsize=25)
ax[0,0].annotate('b)', xy=(0, 1.5), xycoords='axes fraction', fontsize=30,
                horizontalalignment='left', verticalalignment='top')
ax[-1,0].set_xlabel(r'$u_i u_i/TKE$', fontsize=20, labelpad = 10)
ax[-1,1].set_xlabel(r'$u_i u_j/TKE$', fontsize=20, labelpad = 10)

#################################################################################
## METCRAX UNSTABLE: Plot Pdfs of uu,vv,ww and uw,vw,uv in a panel plot
#                    (each row corresponds to 1c,2c,3c)
################################################################################
ffig, ax = plt.subplots(3,2, figsize =(16,6), constrained_layout = True, sharex = 'col')

for a, ds, label in zip(ax, [ds_1c_us,ds_2c_us,ds_3c_us], ['1C', '2C', '3C']):
    ds = ds.stack(index=('time','heights'))
    
    #densities
    for var in ['uu','vv', 'ww']:
        sns.kdeplot(x=ds[var]/ds.tke, ax=a[0], label = var, fill = True, linewidth = 2.5, alpha = 0.3)
    
    for var in ['uw','vw', 'uv']:
        sns.kdeplot(x=ds[var]/ds.tke, ax=a[1], label = var, fill = True, linewidth = 2.5, alpha = 0.3)

 #labels and visuals
    a[0].legend(fontsize = 15)
    a[0].set_xlim(-0.2,2.2)
    a[0].set_ylim(0,8)
    a[0].set_ylabel('Pdf', fontsize = 15, rotation = 90, labelpad = 0)
    a[0].yaxis.set_label_coords(-0.05, 0.5)
    a[0].annotate(label, xy=(-0.12, 1), xycoords='axes fraction', fontsize=25,
                    horizontalalignment='left', verticalalignment='top')
    a[0].grid(True)
 
    a[1].legend(fontsize = 15)
    a[1].set_xlim(-1.2,1.2)
    a[1].axvline(x = 0, color = 'black',lw=2,linestyle=":")
    a[1].set_ylabel('', fontsize = 15, rotation = 90, labelpad = 0)
    a[1].grid(True)
    
    
ax[0,0].set_title('Variances',pad=20,fontsize=25)
ax[0,1].set_title('Covariances',pad=20,fontsize=25)
ax[0,0].annotate('c)', xy=(0, 1.5), xycoords='axes fraction', fontsize=30,
                horizontalalignment='left', verticalalignment='top')
ax[-1,0].set_xlabel(r'$u_i u_i/TKE$', fontsize=20, labelpad = 10)
ax[-1,1].set_xlabel(r'$u_i u_j/TKE$', fontsize=20, labelpad = 10)


