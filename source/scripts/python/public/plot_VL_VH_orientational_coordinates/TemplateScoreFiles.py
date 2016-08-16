#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# Author: Rahel Frick (frick.rahel@gmail.com)
# Author: Jeliazko Jeliazkov (jeliazkov@jhu.edu)

import os
import re
from constants import *
import subprocess
import numpy as np


class TemplateScoreFiles(object):
    '''Determine Angles from pdbs and give coordinates to ScoreFile'''

    def __init__(self, counter, tempfiles):
        self.counter = counter
        self.template_list = []
        self.tempfiles = tempfiles

    def calculate_angles(self):
        models = os.listdir(self.tempfiles[self.counter])
        print models
        os.chdir(self.tempfiles[self.counter])
        for template in models:
            print 'TEMPLATE: ', template
            if re.match('model-[0-9].relaxed.pdb', template): # this breaks often, let's think of a better way? - JJ
                score_file_name = template.rstrip('.pdb') + '_LHOCs.score'
                LHOC_cmd = '%s -database %s -s %s -pack_missing_sidechains false -out:file:score_only %s' %(rosetta_LHOC, rosetta_database, template, score_file_name)
                print "Executing %s..." %(LHOC_cmd),
                rc = subprocess.call(LHOC_cmd + " > /dev/null", shell=True)
                assert rc == 0, "Command %s failed!" %(LHOC_cmd)
                print "...done!"
                #oldname = template.rstrip('.pdb') +'_0001.pdb'
                newname = 'template_' + template[6]
                #template_info = [newname]

                # read new files
                f = open(score_file_name, 'r')
                dat = f.readlines()[2].split()
                d, hoa, loa, pa = 0, 0, 0, 0
                d = float(dat[8])
                hoa = float(dat[9])
                loa = float(dat[10])
                pa = float(dat[11])
                template_info = (newname, d, hoa, loa, pa)

                self.template_list.append(template_info)
        os.chdir('..')

        if self.template_list == []:
            print 'WARNING: Did not find grafted templates!'

        self.template_array =np.array(self.template_list, dtype= [('name', np.str_, 10), ('VL_VH_distance', 'float'), ('VL_VH_opening_angle', 'float'), ('VL_VH_opposite_opening_angle', 'float'), ('VL_VH_packing_angle', 'float')])
