#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 23:32:36 2024

@author: federicagucci
"""

import matplotlib.pyplot as plt
import xarray as xr
import seaborn as sns
import numpy as np

#snohats 
ds_1c_sh = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/Submission/Code&Data/SnoHATS_stable_1c.nc')
ds_2c_sh = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/Submission/Code&Data/SnoHATS_stable_2c.nc')
ds_3c_sh = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/Submission/Code&Data/SnoHATS_stable_3c.nc')

ds_sh = ds_1c_sh.stack(index=('time','sonic')).dropna(dim = 'index')

#metcrax stable
ds_1c = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/DataSamu/Metcrax_anisotropy_samu/MetCrax_stable_1c.nc')
ds_2c = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/DataSamu/Metcrax_anisotropy_samu/MetCrax_stable_2c.nc')
ds_3c = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/DataSamu/Metcrax_anisotropy_samu/MetCrax_stable_3c.nc')

ds = ds_1c.stack(index=('time','heights')).dropna(dim = 'index')

#################################################################################
## Build 2D histograms
#################################################################################

#snohats 

fig = plt.figure(figsize =(12,5))
gs = fig.add_gridspec(2, 4, height_ratios=(0.5,2),width_ratios=(0.5,2,2,2),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.12, hspace=0.05)
Hmin = 1
Hmax = 300
nbin = 30

y = ds_sh.phi.isel({'eigidx':0})

for col,var,xlim,colors in zip([1,2,3],['uu','vv', 'uv'],[[0,2],[0,2],[-1,1]],['C0','C1','C2']):
    ax = fig.add_subplot(gs[1, col])
    ax.tick_params('y', labelleft=False)
    x = ds_sh[var]/ds_sh.tke
    
    #2D histograms
    H, yedges, xedges = np.histogram2d(y, x, bins=nbin)
    H[H==0.0] = np.nan
    pc=ax.pcolormesh(xedges, yedges, H, cmap='viridis',vmin=Hmin,vmax=Hmax)
    ax.axhline(0, color='black',ls='--',lw=1,alpha=0.5)
    ax.set_xlabel(str(var)+'/TKE',fontsize=15)
    ax.set_xlim(xlim)
    ax.set_ylim((-90,90))
    ax.grid(alpha=0.5)
    
    
    ax2 = fig.add_subplot(gs[0, col],xticklabels=[])
    #ax2.tick_params('y', labelleft=False)
    ax2.set_ylim((0,2))
    ax2.set_xlim(xlim)
    ax2.grid(alpha=0.5)
    
    #density plots at margins
    p1=sns.kdeplot(x=x, ax=ax2,color=colors, shade = True)
    p1.set(ylabel = None)
    p1.tick_params(labelleft=False)

    
    if(col == 1):

        ax2.tick_params('y', labelleft=True)
        ax2.annotate('a)', xy=(-0.6, 1.7), xycoords='axes fraction', fontsize=30, 
                        horizontalalignment='left', verticalalignment='top')
        p1.set_ylabel('Pdf',fontsize=15)
        
# Adding the colorbar
cbar_ax = fig.add_axes([0.92, 0.1, 0.01, 0.62])
cbar = fig.colorbar(pc, cax=cbar_ax)
cbar.set_label(label='Occurrence',fontsize=15)
cbar.ax.yaxis.set_ticks_position('right')

ax_M = fig.add_subplot(gs[1, 0],sharey=ax)
ax_M.tick_params('y', labelleft=True)
ax_M.set_ylabel('$\phi$',labelpad=8,fontsize=20, rotation = 0)
ax_M.grid(alpha=0.5)  
ax_M.invert_xaxis()
ax_M.axhline(0, color='black',ls='--',lw=1,alpha=0.5)
ax_M.set_xlabel('Pdf',fontsize=15)  
ax_M.set_xlim((0.011,0))
sns.kdeplot(y=y, ax=ax_M,color='grey', shade = True)


#################################################################################

#metcrax stable
fig = plt.figure(figsize =(12,5))
gs = fig.add_gridspec(2, 4, height_ratios=(0.5,2),width_ratios=(0.5,2,2,2),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.12, hspace=0.05)
Hmin = 1
Hmax = 30 #30 metcrax stable # 300 snohats
nbin = 30

y = ds.phi.isel({'eigidx':0})

for col,var,xlim,colors in zip([1,2,3],['uu','vv', 'uv'],[[0,2],[0,2],[-1,1]],['C0','C1','C2']):
    ax = fig.add_subplot(gs[1, col])
    ax.tick_params('y', labelleft=False)
    x = ds[var]/ds.tke
    
    #2D histograms
    H, yedges, xedges = np.histogram2d(y, x, bins=nbin)
    H[H==0.0] = np.nan
    pc=ax.pcolormesh(xedges, yedges, H, cmap='viridis',vmin=Hmin,vmax=Hmax)
    ax.axhline(0, color='black',ls='--',lw=1,alpha=0.5)
    ax.set_xlabel(str(var)+'/TKE',fontsize=15)
    ax.set_xlim(xlim)
    ax.set_ylim((-90,90))
    ax.grid(alpha=0.5)
    
    
    ax2 = fig.add_subplot(gs[0, col],xticklabels=[])
    ax2.set_ylim((0,2))
    ax2.set_xlim(xlim)
    ax2.grid(alpha=0.5)
    
    #density plots at margins
    p1=sns.kdeplot(x=x, ax=ax2,color=colors, shade = True)
    p1.set(ylabel = None)
    p1.tick_params(labelleft=False)

    
    if(col == 1):

        ax2.tick_params('y', labelleft=True)
        ax2.annotate('b)', xy=(-0.6, 1.7), xycoords='axes fraction', fontsize=30,

                        horizontalalignment='left', verticalalignment='top')
        p1.set_ylabel('Pdf',fontsize=15)
        
# Adding the colorbar
cbar_ax = fig.add_axes([0.92, 0.1, 0.01, 0.62])
cbar = fig.colorbar(pc, cax=cbar_ax)
cbar.set_label(label='Occurrence',fontsize=15)
cbar.ax.yaxis.set_ticks_position('right')

ax_M = fig.add_subplot(gs[1, 0],sharey=ax)
ax_M.tick_params('y', labelleft=True)
ax_M.set_ylabel('$\phi$',labelpad=8,fontsize=20, rotation = 0)
ax_M.grid(alpha=0.5)  
ax_M.invert_xaxis()
ax_M.axhline(0, color='black',ls='--',lw=1,alpha=0.5)
ax_M.set_xlabel('Pdf',fontsize=15)  
ax_M.set_xlim((0.011,0))
sns.kdeplot(y=y, ax=ax_M,color='grey', shade = True)
