#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 20:32:12 2024

@author: sps
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# Load PSF and DCD files
u = mda.Universe('data-p.psf', 'eq-1.dcd', format="LAMMPS")
num_frames = len(u.trajectory)  # Number of frames in the trajectory
num_bins = 2400

# iterate across trajectory
for ts in u.trajectory[0:1]:
    # Select atoms named 'NIT'
    ag = u.select_atoms('resname NIT',updating=True)
    # bin-width
    dr = 1

    # Load the distancesezzzzzz
    distances = ag.positions

    # convert to radial distance
    distance = np.linalg.norm(distances,axis=1)
    
    max = np.max(distance)
    #create bins
    bins = np.arange(0,np.max(distance)+dr,dr)
    # find out bin indices
    indices = np.digitize(distance, bins)
    count = np.bincount(indices)

    #array = [ts,count]
    #average = np.mean(count,axis = 0)   # print(sum(counts))
    
    # calculate volume
    vol = 4/3 * np.pi * np.array([i**3 for i in bins])
    
    # calculate the differential volume
    dvol = np.diff(vol)
     
    # calculate number density
    ndens = count[1:] / dvol
    mdens = ndens * 10**4 * 14 / 6.023
plt.ylim(0,1200)
plt.xlim(0, 40)
plt.grid()
plt.plot(bins[1:],mdens)

threshold = 0.25*np.max(mdens)

for i in range(20, int(max)): 
    if mdens[i] < threshold: 
        print('The radius is:', i, "Angstrom")
        break

radius  = np.zeros(len(u.trajectory))


def radii_calc():    
    for q in range(len(u.trajectory)): 
        for ts in u.trajectory[q:q+1]:       
            updated = ag.select_atoms('sphzone 28 resname NIT', updating = True)
            com = updated.center_of_mass()
            
            distances = updated.positions - com
            # convert to radial distance
            distance = np.linalg.norm(distances,axis=1)
            
            #create bins
            bins = np.arange(0,np.max(distance)+dr,dr)
            # find out bin indices
            indices = np.digitize(distance, bins)
            count = np.bincount(indices)

            # calculate volume
            vol = 4/3 * np.pi * np.array([i**3 for i in bins])
            
            # calculate the differential volume
            dvol = np.diff(vol)
            
            # calculate number density
            ndens = count[1:] / dvol
            mdens = ndens * 10**4 * 14 / 6.023
            threshold = 0.25*np.max(mdens)
            for i in range(20, int(np.max(distance))): 
                if mdens[i] < threshold: 
                    r = np.cbrt((3/(4*np.pi))*vol[i])
                    radius[q] = r
                    break 
            

radii_calc()