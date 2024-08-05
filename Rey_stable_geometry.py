#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 16:26:43 2024

@author: federicagucci
"""
import xarray as xr
import numpy as np
import copy

#Eigenvectors: 
#   compute eigenvectors with our convention explained in Sec.3.1 of the paper
#   compute eigenvalues and sort them in descending order
#   Λ1: largest λ[:,:,0]   Λ2: medium λ[:,:,1]     Λ3:smallest λ[:,:,2]

def Eigenvectors(UU,VV,WW,UV,UW,VW,nbin):
    print(nbin)
    small = []
    medium = []
    large = []
    lambdas = []
    for i in range(0,nbin,1):
        Rey = np.array([            
             [UU[i], UV[i] , UW[i]],
             [UV[i], VV[i] ,VW[i]],
             [UW[i], VW[i],  WW[i]]           

        ])
        
        if(np.isnan(Rey).any()):
            
            eigenvS = [np.nan,np.nan,np.nan]
            eigenvM = [np.nan,np.nan,np.nan]
            eigenvL = [np.nan,np.nan,np.nan]
            eigvals = [np.nan,np.nan,np.nan]
        
        else:
            eigvals, eigvecs = np.linalg.eigh(Rey) #eigenvalues in ascending order (opposite to R)
            
            #eigenvalues of bij = (eigenvalues of Rey - 2TKE/3)/2TKE
            eigvals = (eigvals/np.sum(eigvals) - 1./3. )
            
            #bij and Rey same set of eigenvectors
            eigenvS=eigvecs[:,0]
            eigenvL=eigvecs[:,2]
            
            if(eigenvS[2]<0.0):
                eigenvS = eigenvS*(-1.0)
            if(eigenvL[0]<0.0):
                eigenvL = eigenvL*(-1.0)

            eigenvM=np.cross(eigenvL,eigenvS)*(-1.0)
            
        small.append(eigenvS)
        medium.append(eigenvM)
        large.append(eigenvL)
        
        lambdas.append(eigvals[::-1]) #sort in descending order
        
    #eigenvalues and eigenvectors of the Reynolds stress tensor   
    return np.array(large),np.array(medium),np.array(small), np.array(lambdas)

#EllipsoidAngles: 
#   compute the inclination angle beta (see Sec.3.1): 
#               the angle between a plane (spanned by Λ1,Λ2) and a vector (U) is the complement of the
#               angle between the normal to the plane (Λ3) and the vector
#   compute the orientation angle phi (see Sec.3.1) for the three eigenvectors
#               (in the paper only phi of Λ1 is used)
#   compute the angle theta (between U and the vertical) for the three eigenvectors
#               (not used in the paper)
def EllipsoidAngles(v):
    
    beta = 90. - (np.arccos(v.isel({'xyz' : 0, 'eigidx' : 2})))*180./np.pi
    
    phi = (np.arctan(v.isel({'xyz' : 1})/v.isel({'xyz' : 0})))*180./np.pi
    
    theta = (np.arccos(v.isel({'xyz' : 2})))*180./np.pi
    
    return beta, phi, theta

#xyBarycentricMap:
#   compute the coordinate x,y of the barycentric map explained in Sec.3 (Eq. 2)
def xyBarycentricMap(evals):
    #eigenvalues are given in descending order to the function

    xB = []
    yB = []
   
    xB = -2.*evals.isel({'eigidx' : 1}) + 0.5*(evals.isel({'eigidx' : 2}) + 1.)
    yB = (np.sqrt(3.)/2.)*(3.*evals.isel({'eigidx' : 2}) + 1.)
 
    return xB, yB

# FindBarycentricPoints:
#   define the 3 kite-shaped regions in the barycentric map (defined as sum of two triangles) (Sec.3, Fig.1)
#               corresponding to 1c,2c,isotropic Reynolds stress tensor
#
#   check if the coordinates x,y fall in these tegions
#   return nan if not, otherwise return the index (of time) when this occurs
def FindBarycentricPoints(p_x, p_y, anis, perc):
   
    triangle = {}
    center = {}

    triangle['x'] = np.array([0, 0.5 * np.cos(60 * np.pi / 180), 0.5, 1 - 0.5 * np.cos(60 * np.pi / 180), 1, 0.5, 0])
    triangle['y'] = np.array([0, 0.5 * np.sin(60 * np.pi / 180), np.sqrt(3) / 2, 0.5 * np.sin(60 * np.pi / 180), 0, 0, 0])

    center['x'] = 0.5
    center['y'] = np.sqrt(3) / 6

    c = np.sqrt(center['x']**2 + center['y']**2)

    pc = {'x': np.full((p_x.shape[0], 6), np.nan), 'y': np.full((p_x.shape[0], 6), np.nan)}

    if anis == "2c":
        # ... points of the first triangle for blue
        pc['x'][:, 0] = triangle['x'][0]
        pc['y'][:, 0] = triangle['y'][0]

        pc['x'][:, 1] = triangle['x'][1] * (1 - perc / 100)
        pc['y'][:, 1] = triangle['y'][1] * (1 - perc / 100)

        pc['x'][:, 2] = c * (1 - perc / 100) * np.cos(30 * np.pi / 180)
        pc['y'][:, 2] = c * (1 - perc / 100) * np.sin(30 * np.pi / 180)

        # ... points of the second triangle for blue
        pc['x'][:, 3] = pc['x'][:, 2]
        pc['y'][:, 3] = pc['y'][:, 2]

        pc['x'][:, 4] = triangle['x'][5] * (1 - perc / 100)
        pc['y'][:, 4] = triangle['y'][5]

        pc['x'][:, 5] = pc['x'][:, 0]
        pc['y'][:, 5] = pc['y'][:, 0]

    elif anis == "iso":
        # ... points of the first triangle for green

        pc['x'][:, 0] = center['x']
        pc['y'][:, 0] = triangle['y'][2] - c * (1 - perc / 100)

        pc['x'][:, 1] = triangle['x'][1] * (1 + perc / 100)
        pc['y'][:, 1] = triangle['y'][1] * (1 + perc / 100)

        pc['x'][:, 2] = triangle['x'][2]
        pc['y'][:, 2] = triangle['y'][2]
        
        # ... points of the second triangle for green

        pc['x'][:, 3] = pc['x'][:, 2]
        pc['y'][:, 3] = pc['y'][:, 2]

        pc['x'][:, 4] = 1 - (1 + perc / 100) + triangle['x'][3] * (1 + perc / 100)
        pc['y'][:, 4] = triangle['y'][3] * (1 + perc / 100)

        pc['x'][:, 5] = pc['x'][:, 0]
        pc['y'][:, 5] = pc['y'][:, 0]

    elif anis == "1c":
        # ... points of the first triangle for red
        pc['x'][:, 0] = c * (1 + perc / 100) * np.cos(30 * np.pi / 180)
        pc['y'][:, 0] = c * (1 - perc / 100) * np.sin(30 * np.pi / 180)

        pc['x'][:, 1] = 1 - (1 - perc / 100) + triangle['x'][3] * (1 - perc / 100)
        pc['y'][:, 1] = triangle['y'][3] * (1 - perc / 100)

        pc['x'][:, 2] = triangle['x'][4]
        pc['y'][:, 2] = triangle['y'][4]
        
        # ... points of the second triangle for red

        pc['x'][:, 3] = pc['x'][:, 2]
        pc['y'][:, 3] = pc['y'][:, 2]

        pc['x'][:, 4] = triangle['x'][5] * (1 + perc / 100)
        pc['y'][:, 4] = triangle['y'][5]

        pc['x'][:, 5] = pc['x'][:, 0]
        pc['y'][:, 5] = pc['y'][:, 0]

    alpha1 = ((pc['y'][:, 1] - pc['y'][:, 2]) * (p_x - pc['x'][:, 2]) + (pc['x'][:, 2] - pc['x'][:, 1]) * (p_y - pc['y'][:, 2])) / \
             ((pc['y'][:, 1] - pc['y'][:, 2]) * (pc['x'][:, 0] - pc['x'][:, 2]) + (pc['x'][:, 2] - pc['x'][:, 1]) * (pc['y'][:, 0] - pc['y'][:, 2]))
    beta1 = ((pc['y'][:, 2] - pc['y'][:, 0]) * (p_x - pc['x'][:, 2]) + (pc['x'][:, 0] - pc['x'][:, 2]) * (p_y - pc['y'][:, 2])) / \
            ((pc['y'][:, 1] - pc['y'][:, 2]) * (pc['x'][:, 0] - pc['x'][:, 2]) + (pc['x'][:, 2] - pc['x'][:, 1]) * (pc['y'][:, 0] - pc['y'][:, 2]))
    gamma1 = 1.0 - alpha1 - beta1

    alpha2 = ((pc['y'][:, 4] - pc['y'][:, 5]) * (p_x - pc['x'][:, 5]) + (pc['x'][:, 5] - pc['x'][:, 4]) * (p_y - pc['y'][:, 5])) / \
             ((pc['y'][:, 4] - pc['y'][:, 5]) * (pc['x'][:, 3] - pc['x'][:, 5]) + (pc['x'][:, 5] - pc['x'][:, 4]) * (pc['y'][:, 3] - pc['y'][:, 5]))
    beta2 = ((pc['y'][:, 5] - pc['y'][:, 3]) * (p_x - pc['x'][:, 5]) + (pc['x'][:, 3] - pc['x'][:, 5]) * (p_y - pc['y'][:, 5])) / \
            ((pc['y'][:, 4] - pc['y'][:, 5]) * (pc['x'][:, 3] - pc['x'][:, 5]) + (pc['x'][:, 5] - pc['x'][:, 4]) * (pc['y'][:, 3] - pc['y'][:, 5]))
    gamma2 = 1.0 - alpha2 - beta2

    ind1 = np.where((alpha1 > 0.0) & (beta1 > 0.0) & (gamma1 > 0.0))[0]
    ind2 = np.where((alpha2 > 0.0) & (beta2 > 0.0) & (gamma2 > 0.0))[0]
    ind = np.sort(np.concatenate((ind1, ind2)))


    return ind

#1 min dataset
ds = xr.open_dataset('/Volumes/remember/NikkiVercauteren/2ndWork/Submission/Code&Data/SnoHATS.nc')

#select only data in ds corresponding to stable cases
ds_stable = ds.where(ds.wT < 0.0)

Λ3 = np.zeros((len(ds_stable.time), len(ds_stable.sonic),3))
Λ2 = np.zeros((len(ds_stable.time), len(ds_stable.sonic),3))
Λ1 = np.zeros((len(ds_stable.time), len(ds_stable.sonic),3))
λ = np.zeros((len(ds_stable.time), len(ds_stable.sonic),3))

for sn in range(0,len(ds_stable.sonic),1):
   
    uu = ds_stable.uu[:,sn].data
    vv = ds_stable.vv[:,sn].data
    ww = ds_stable.ww[:,sn].data
    
    uv = ds_stable.uv[:,sn].data
    uw = ds_stable.uw[:,sn].data
    vw = ds_stable.vw[:,sn].data
    
    Λ1[:,sn,:], Λ2[:,sn,:], Λ3[:,sn,:], λ[:,sn,:] = Eigenvectors(uu,vv,ww,uv,uw,vw,len(uu))

ds_stable = ds_stable.assign(
    eigenvectors = (['eigidx', 'time','sonic','xyz'], [Λ1,Λ2,Λ3]),
    eigenvalues = (['time','sonic','eigidx'], λ )
    )    

ds_stable['eigenvectors'] = ds_stable.eigenvectors.transpose('time','sonic','eigidx','xyz')

β, Φ, θ = EllipsoidAngles(ds_stable.eigenvectors)

ds_stable = ds_stable.assign(
    beta = ([ 'time','sonic'], β.data),
    phi = ([ 'time','sonic','eigidx'], Φ.data),
    theta = ([ 'time','sonic','eigidx'],  θ.data)
    
    )  
#ds_stable.phi = ds_stable.phi.transpose('time','sonic','eigidx')
#ds_stable.theta = ds_stable.theta.transpose('time','sonic','eigidx')

xb, yb = xyBarycentricMap(ds_stable.eigenvalues)

ds_stable = ds_stable.assign(
    
    xb = ([ 'time','sonic'], xb.data),
    yb = ([ 'time','sonic'], yb.data)
    
    )  

#percentage used to define the kite-shaped region in the barycentric map
#   (convention used as in Stiperski and Calaf, QJRMS, 2018, Figure 4)
percentage = 30 

ds_1c = copy.deepcopy(ds_stable)
ds_2c = copy.deepcopy(ds_1c)
ds_3c = copy.deepcopy(ds_2c)

one = []
two = []
iso = []
for sn in range(0,len(ds_stable.sonic),1):
    
    #array with indices (of time) where 1c,2c, or isotropic Reynolds stress tensor are found 
    one_sn = FindBarycentricPoints(ds_stable.xb.isel({'sonic' : sn}),ds_stable.yb.isel({'sonic' : sn}),"1c",percentage) 
    two_sn = FindBarycentricPoints(ds_stable.xb.isel({'sonic' : sn}),ds_stable.yb.isel({'sonic' : sn}),"2c",percentage) 
    iso_sn = FindBarycentricPoints(ds_stable.xb.isel({'sonic' : sn}),ds_stable.yb.isel({'sonic' : sn}),"iso",percentage) 
    
    #invalidate data as NaN if they do not fall in kite-shaped regions
    trash_1c = np.setdiff1d(np.arange(len(ds_stable.time)), one_sn)
    trash_2c = np.setdiff1d(np.arange(len(ds_stable.time)), two_sn)
    trash_3c = np.setdiff1d(np.arange(len(ds_stable.time)), iso_sn)
    
    for ii in ['uu','vv','ww','uw','uv','vw','wT','tke','zeta','beta','xb','yb']:
        
        ds_1c[ii][trash_1c,sn] = np.nan
        ds_2c[ii][trash_2c,sn] = np.nan
        ds_3c[ii][trash_3c,sn] = np.nan
        
    for ii in ['phi','theta','eigenvalues']:
        ds_1c[ii][trash_1c,sn,:] = np.nan
        ds_2c[ii][trash_2c,sn,:] = np.nan
        ds_3c[ii][trash_3c,sn,:] = np.nan
    
    ds_1c['eigenvectors'][trash_1c ,sn,:,:] = np.nan
    ds_2c['eigenvectors'][trash_2c,sn,:,:] = np.nan
    ds_3c['eigenvectors'][trash_3c,sn,:,:] = np.nan
    
    one.append(one_sn)
    two.append(two_sn)
    iso.append(iso_sn)

ds_stable.to_netcdf('/Volumes/remember/NikkiVercauteren/2ndWork/Submission/Code&Data/SnoHATS_stable.nc')
ds_1c.to_netcdf('/Volumes/remember/NikkiVercauteren/2ndWork/Submission/Code&Data/SnoHATS_stable_1c.nc')
ds_2c.to_netcdf('/Volumes/remember/NikkiVercauteren/2ndWork/Submission/Code&Data/SnoHATS_stable_2c.nc')
ds_3c.to_netcdf('/Volumes/remember/NikkiVercauteren/2ndWork/Submission/Code&Data/SnoHATS_stable_3c.nc')
