#!/usr/bin/env python
#@file: energy_landscape_metrics.py
#@author: Rebecca F. Alford (ralford3@jhu.edu)
#@brief: Given an energy landscape for a transmembrane peptide, calculate the minimum energy orientaiton and the approximate dG of transfer from water to bilayer

import sys, os
import numpy as np

def compute_dG_transfer_energy( zcoords, angles, energies ):

    # Find the index of z = 0, angle = 0
    z_is_zero_indices = set(np.where( abs(zcoords) < 1 )[0])
    angles_is_zero_indices = set(np.where( angles == 0.0 )[0])
    first_index =list( z_is_zero_indices.intersection(angles_is_zero_indices) )

    bilayer_energy = energies[first_index[0]]

    # Find the index of zmax, ange = 270; RS changing it to 90 to have it within 0-90 degrees
    angle_is_90_indices = set(np.where( angles == 90.0 )[0])
    z_is_max = set(np.where( zcoords == np.max(zcoords) )[0] )
    second_index = list( angle_is_90_indices.intersection( z_is_max ))
    if len(second_index) > 0 and second_index[0] in energies:
        solution_state_energy = energies[second_index[0]]
    else:
        print("Unable to compute dG transfer energy -- data not found.")
        return 999999

    return round( bilayer_energy - solution_state_energy, 2 )

def compute_minimum_energy_orientation( zcoords, angles, energies ):
    loc = np.argmin( energies )
    return zcoords[loc], angles[loc]
