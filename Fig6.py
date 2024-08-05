#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 20:52:37 2024

@author: federicagucci
"""

import matplotlib.pyplot as plt
import xarray as xr
import seaborn as sns
import numpy as np

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


datasets = [[ds_1c_sh, ds_2c_sh, ds_3c_sh],[ds_1c, ds_2c, ds_3c],[ds_1c_us, ds_2c_us, ds_3c_us]]

#################################################################################
## Plot beta(zeta) for 1c,2c,3c for each dataset and plot dataset in a panel plot 
#################################################################################

color = sns.color_palette('Set2', 10)

#binmedians
min, max, n = -2, 3, 10
exp_bin = np.linspace(min, max, n)
exp_labels = (exp_bin[:-1] + exp_bin[1:]) / 2
bins_zeta = 10 ** (exp_bin)
labels_zeta = 10 ** exp_labels

fig, axs = plt.subplots(1,3,figsize =(16,6), constrained_layout=True)

for ax, data in zip(axs, datasets):
    for ds, c, label in zip(data, [color[1],color[2],color[5]], ['1C', '2C', '3C']):

        ds_groups = ds.groupby_bins(ds.zeta, bins = bins_zeta, labels = labels_zeta)
        bin_count = ds_groups.count().zeta
        ds_median = ds_groups.median().where(bin_count>10)
        ds_UQ = ds_groups.quantile(0.25).where(bin_count>10)
        ds_LQ = ds_groups.quantile(0.75).where(bin_count>10)

        bars = np.array([ds_median.beta.data-ds_LQ.beta.data,ds_UQ.beta.data-ds_median.beta.data])
        ax.errorbar(ds_median.zeta_bins, ds_median.beta, yerr = abs(bars),
                    color = c, linewidth = 3, label = label, capsize = 10, capthick = 3,
                   marker = 'o', markeredgecolor = 'black', markersize = 10)
        ax.legend(fontsize = 20)
        ax.set_xlabel(r'$\zeta$', fontsize = 40)
        ax.grid(True)
        ax.set_xscale('log')
        ax.set_ylim(-25,35)
        
        #klipp angle
        ax.hlines(17, 1e-2, 1e3, color = 'black', linestyle = 'dashed')
        
        #zero line
        ax.hlines(0, 1e-2, 1e3, color = 'black')
        
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(20)  

        ax.set_xlabel(r'$\zeta$', fontsize = 25)

axs[0].set_title('a)', fontsize = 30, loc = 'left', pad = 20)
axs[1].set_title('b)', fontsize = 30, loc = 'left', pad = 20)
axs[2].set_title('c)', fontsize = 30, loc = 'left', pad = 20)  
axs[0].set_ylabel(r'$\beta$', fontsize = 40, rotation = 0, labelpad = 20)
axs[2].set_xlabel(r'$-\zeta$', fontsize = 25)
plt.subplots_adjust(wspace=0.1)