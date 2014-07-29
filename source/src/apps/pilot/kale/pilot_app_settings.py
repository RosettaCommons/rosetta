# -*- mode:python;indent-tabs-mode:nil;show-trailing-whitespace:t; -*-
# vim: ft=python
#
# Project settings for rosetta sources
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

# Do Not Modify this file and check it in.  First copy this file to 
# pilot_apps.src.settings.my, then modify it to include the sources of the 
# pilot applications you want built. Remember to create a Remember also to 
# setup the devel.src.settings.my file to include any experimental files that 
# your pilot application needs

sources = {
    'pilot/kale': [
        'examples/current_example',
        #'cyclic_poses/load_cyclic_pose',
        #'cyclic_poses/kick_cyclic_pose',
        #'detailed_balance/backbone_tests',
        #'detailed_balance/sidechain_tests',
        #'native_ensemble/native_ensemble',
        'native_ensemble/analysis/trajectory_movie',
        'native_ensemble/analysis/query_trajectory',
        'kic_refactor/KicSandbox',
        ]
}

include_path = []
library_path = []
libraries    = []
subprojects  = [
        'devel', 
        'protocols.7',
        'protocols.6', 
        'protocols_c.5', 
        'protocols_b.5', 
        'protocols_a.5', 
        'protocols_f.4', 
        'protocols_e.4', 
        'protocols_d.4', 
        'protocols_c.4', 
        'protocols_b.4', 
        'protocols_a.4',
        'protocols.3',
        'protocols_b.2',
        'protocols_a.2',
        'protocols.1',
        'core.5',
        'core.4',
        'core.3',
        'core.2',
        'core.1',
        'basic',
        'numeric',
        'utility',
        'ObjexxFCL',
        'z',
        'cppdb',
        'sqlite3',
]
