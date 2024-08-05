#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 19:17:50 2024

@author: federicagucci
"""

import seaborn as sns
from matplotlib import pyplot as plt 
import numpy as np
import pandas as pd
import xarray as xr

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
## SnoHATS: Plot Pdf of beta for 1c,2c,3c in a panel plot
################################################################################
# Prepare dataset
################## 

ds = ds_1c_sh.stack(index=('time','sonic'))

df1 = pd.DataFrame(
    {'beta':ds.beta,
     'C':np.repeat('1C',len(ds.index)),
     'beta_median': np.repeat(ds.beta.median(skipna = True).data,len(ds.index))
     },)

ds = ds_2c_sh.stack(index=('time','sonic'))

df2 = pd.DataFrame(
    {'beta':ds.beta,
     'C':np.repeat('2C',len(ds.index)),
     'beta_median': np.repeat(ds.beta.median().data,len(ds.index))
     },)

ds = ds_3c_sh.stack(index=('time','sonic'))

df3 = pd.DataFrame(
    {'beta':ds.beta,
     'C':np.repeat('3C',len(ds.index)),
     'beta_median': np.repeat(ds.beta.median().data,len(ds.index))
     },)

df_sh = pd.concat([df1,df2,df3])

##################################################################################
# Plotting dataset
#################
palette_tab10 = sns.color_palette('Set2', 10)
palette = sns.color_palette([palette_tab10[1], palette_tab10[2], palette_tab10[5]])


#theme
sns.set_theme(style="whitegrid", rc={"axes.facecolor": (0, 0, 0, 0)})


#grid
g = sns.FacetGrid(df_sh, col = 'C', hue='C',palette = palette, aspect=7, height=0.8)

#density plots
g.map(sns.kdeplot, 'beta',
      bw_adjust =0.5,
      fill=True, alpha=0.7, linewidth=1.5,
      common_norm = False)
g.map(sns.kdeplot, 'beta',
      bw_adjust =0.5,#to increase or decrease the amount of smoothing.
      fill=False, color='black', linewidth=1.5,
      common_norm = False)

#beta from klipp
g.map(plt.axvline, x=17,
      lw=2, color = 'black', linestyle="dashdot")

#median lines
def vline(x,**kwargs):
    return plt.vlines(x,0,0.1, **kwargs) #3rd argument:length of median line

g.map(vline, 'beta_median',lw = 2,color='black')

# labels
for i, ax in enumerate(g.axes[:,0]):
    ax.set_ylabel('', fontsize=20, rotation = 0, labelpad = 25, loc ='bottom')
    ax.set_xlim([-40,50])
 
#ticks sizes 
for ax in g.axes[-1,:]:
    plt.setp(ax.get_xticklabels(), fontsize=20)
    ax.set_xlabel(r'$\beta$', fontweight='bold', fontsize=35)
    
#layout
g.fig.subplots_adjust(hspace=-0.4, wspace = 0.08) #this is for height spacing between rows
g.set_titles("")
g.set(yticks=[])
g.despine(left=True)

#titles
axs = g.axes
axs[0,0].set_title('1C', fontsize = 30)
axs[0,1].set_title('2C', fontsize = 30)
axs[0,2].set_title('3C', fontsize = 30)
g.fig.text(0.07,1.6, 
            #verticalalignment='center', #make sure it's aligned at center vertically
            s='a)',
            color = 'black',
            fontsize =50, 
            rotation=0)

#################################################################################
## METCRAX stable: Mountain plot with Pdf of beta 
#                  (Panel plot of 1c,2c,3c for each height w)
################################################################################
# Prepare dataset
################## 
heights = ds_1c.heights.data

#component index

data = xr.concat(   [ds_1c.assign_coords(C = '1c').stack(index = ('time','heights')).reset_index('index').dropna(dim='index'),
                    ds_2c.assign_coords(C = '2c').stack(index = ('time','heights')).reset_index('index').dropna(dim='index'),
                    ds_3c.assign_coords(C = '3c').stack(index = ('time','heights')).reset_index('index').dropna(dim='index')]
                , dim = 'index')

#pandas dataframe
df = data[['heights','beta', 'C']].to_dataframe()
df['beta_median'] = df.groupby(['heights','C'])['beta'].transform('median')



##################################################################################
# Plotting dataset
#################

#theme
sns.set_theme(style="whitegrid", rc={"axes.facecolor": (0, 0, 0, 0)})
#grid
g = sns.FacetGrid(df, row = 'heights', col = 'C', hue='C', palette = palette, row_order=heights[::-1], aspect=5, height=0.75)


#density plots
g.map(sns.kdeplot, 'beta',
      bw_adjust =0.5,
      fill=True, alpha=0.7, linewidth=1.5,
      common_norm = False)
g.map(sns.kdeplot, 'beta',
      bw_adjust =0.5,#to increase or decrease the amount of smoothing.
      fill=False, color='black', linewidth=1.5,
      common_norm = False)

#beta from klipp
g.map(plt.axvline, x=17,
      lw=2, color = 'black', linestyle="dashdot")

#median lines
def vline(x,**kwargs):
    return plt.vlines(x,0,0.035, **kwargs) #3rd argument:length of median line

g.map(vline, 'beta_median',lw = 2,color='black')

# labels
for i, ax in enumerate(g.axes[:,0]):
    ax.set_ylabel(heights[::-1][i], fontsize=20, rotation = 0, labelpad = 25, loc ='bottom')
    ax.set_xlim([-40,50]) 
#ticks sizes 
for ax in g.axes[-1,:]:
    plt.setp(ax.get_xticklabels(), fontsize=20)
    ax.set_xlabel(r'$\beta$', fontweight='bold', fontsize=35)
    
#layout
g.fig.subplots_adjust(hspace=-0.4)
g.set_titles("")
g.set(yticks=[])
g.despine(left=True)

#titles
axs = g.axes
axs[0,0].set_title('1C', fontsize = 30)
axs[0,1].set_title('2C', fontsize = 30)
axs[0,2].set_title('3C', fontsize = 30)

# overall ylabel
g.fig.text(0,0.5, 
            #verticalalignment='center', #make sure it's aligned at center vertically
            s='Z',
            color = 'black',
            fontsize =35, 
            rotation=0)

g.fig.text(0,1.02, 
            #verticalalignment='center', #make sure it's aligned at center vertically
            s='b)',
            color = 'black',
            fontsize =50, 
            rotation=0)

#################################################################################
## METCRAX unstable: Mountain plot with Pdf of beta 
#                  (Panel plot of 1c,2c,3c for each height w)
################################################################################
# Prepare dataset
#################
#component index

data = xr.concat(   [ds_1c_us.assign_coords(C = '1c').stack(index = ('time','heights')).reset_index('index').dropna(dim='index'),
                    ds_2c_us.assign_coords(C = '2c').stack(index = ('time','heights')).reset_index('index').dropna(dim='index'),
                    ds_3c_us.assign_coords(C = '3c').stack(index = ('time','heights')).reset_index('index').dropna(dim='index')]
                , dim = 'index')

#pandas dataframe
df = data[['heights','beta', 'C']].to_dataframe()
df['beta_median'] = df.groupby(['heights','C'])['beta'].transform('median')

df.loc[(df.heights == 3) & (df.C == '1c'), 'beta_median'] = -9999
df.loc[(df.heights == 3) & (df.C == '3c'), 'beta_median'] = -9999

#################################################################################
# Plotting dataset
#################

#theme
sns.set_theme(style="whitegrid", rc={"axes.facecolor": (0, 0, 0, 0)})
#grid
g = sns.FacetGrid(df, row = 'heights', col = 'C', hue='C', palette = palette, row_order=heights[::-1], aspect=5, height=0.75)


#density plots
g.map(sns.kdeplot, 'beta',
      bw_adjust =0.5,
      fill=True, alpha=0.7, linewidth=1.5,
      common_norm = False)
g.map(sns.kdeplot, 'beta',
      bw_adjust =0.5,#to increase or decrease the amount of smoothing.
      fill=False, color='black', linewidth=1.5,
      common_norm = False)

#beta from klipp
g.map(plt.axvline, x=17,
      lw=2, color = 'black', linestyle="dashdot")

#median lines
def vline(x,**kwargs):
    return plt.vlines(x,0,0.035, **kwargs) #3rd argument:length of median line

g.map(vline, 'beta_median',lw = 2,color='black')

# labels
for i, ax in enumerate(g.axes[:,0]):
    ax.set_ylabel(heights[::-1][i], fontsize=20, rotation = 0, labelpad = 25, loc ='bottom')
    ax.set_xlim([-15,25])#unstable
 
#ticks sizes 
for ax in g.axes[-1,:]:
    plt.setp(ax.get_xticklabels(), fontsize=20)
    ax.set_xlabel(r'$\beta$', fontweight='bold', fontsize=35)
    
#layout
g.fig.subplots_adjust(hspace=-0.4)
g.set_titles("")
g.set(yticks=[])
g.despine(left=True)

#titles
axs = g.axes
axs[0,0].set_title('1C', fontsize = 30)
axs[0,1].set_title('2C', fontsize = 30)
axs[0,2].set_title('3C', fontsize = 30)

# overall ylabel
g.fig.text(0,0.5, 
            #verticalalignment='center', #make sure it's aligned at center vertically
            s='Z',
            color = 'black',
            fontsize =35, 
            rotation=0)

g.fig.text(0,1.02, 
            #verticalalignment='center', #make sure it's aligned at center vertically
            s='c)',
            color = 'black',
            fontsize =50, 
            rotation=0)


