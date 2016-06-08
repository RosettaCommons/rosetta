## (c) Copyright Rosetta Commons Member Institutions.
## (c) This file is part of the Rosetta software suite and is made available under license.
## (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
## (c) For more information, see http://www.rosettacommons.org. Questions about this can be
## (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

# Author: Rahel Frick (frick.rahel@gmail.com)
# Author: Jeliazko Jeliazkov (jeliazkov@jhu.edu)

import os
import sys

# coordinate stdevs and means for plotting ideal distributions
# TODO calculate these from angles.sc file
lhoc_mus = { "dist": 14.6, # angstroms
                "Hopen": 97.2, # degrees
                "Lopen": 99.4, # degrees
                "pack": -52.3, # degrees
                }
                
lhoc_sigmas = { "dist": 0.34, # angstroms
                "Hopen": 2.63, # degrees
                "Lopen": 1.92, # degrees
                "pack": 3.86, # degrees
                }

# coordinate abbreviations
PA = 'VL_VH_packing_angle'
HOA = 'VL_VH_opening_angle'
LOA = 'VL_VH_opposite_opening_angle'
D = 'VL_VH_distance'
total = 'total_score'
name = 'description'
coordinates = [PA, LOA, HOA, D]

# paths and command lines
rosetta_path = os.environ.get('ROSETTA')
if rosetta_path:
    angles_file = rosetta_path + '/tools/antibody/angles.sc'
    rosetta_LHOC = rosetta_path + '/main/source/bin/packing_angle.macosclangrelease'
    rosetta_database = rosetta_path + '/main/database/'
else:
    sys.exit("""
    Could not find environment variable ROSETTA.
    Please `export ROSETTA=/path/to/Rosetta`.
    Thanks. Exiting.
    """)


# color codes for models (What template do they come from?)
color_dict = {}
color_dict['0'] = "black"
color_dict['1'] = "grey"
color_dict['2'] = "blue"
color_dict['3'] = "cyan"
color_dict['4'] = "green"
color_dict['5'] = "yellow"
color_dict['6'] = "orange"
color_dict['7'] = "red"
color_dict['8'] = "magenta"
color_dict['9'] = "purple"



# xlim for plots
x_lower = { LOA: 85, HOA: 87, PA: -70, D: 13}
x_upper = { LOA: 110, HOA: 115, PA: -35, D: 17}

alpha1 = ['101', '106', '107']
alpha2 = ['206', '217', '218', '220', '221', '223', '226', '228']
