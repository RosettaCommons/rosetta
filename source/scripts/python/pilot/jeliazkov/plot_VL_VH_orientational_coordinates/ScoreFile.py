## (c) Copyright Rosetta Commons Member Institutions.
## (c) This file is part of the Rosetta software suite and is made available under license.
## (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
## (c) For more information, see http://www.rosettacommons.org. Questions about this can be
## (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

# Author: Rahel Frick (frick.rahel@gmail.com)
# Author: Jeliazko Jeliazkov (jeliazkov@jhu.edu)

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from constants import *
import os
from TemplateScoreFiles import TemplateScoreFiles




class LHOCDecoys(object):

    def __init__(self, decoy_array):
        self.decoy_array = decoy_array
        self.template_no = decoy_array[name][6]

    def get_coordinate(self, coordinate):
        return self.decoy_array[coordinate]


class ScoreFile(object):

    def __init__(self, counter, infiles, outpath, names):

        self.infiles = [item.rstrip() for item in infiles]
        self.counter = counter

        #set name if possible
        self.name = ''
        if len(names) == len(infiles):
            self.name = names[self.counter].rstrip('\n') + '_'

        # change the path to the pdb files if necessary
        scorefile_name = self.infiles[self.counter].split('/')[-1]

        self.outpath = outpath
        #self.outpath = self.infiles[self.counter].rstrip(scorefile_name) + 'LHOC'


        try:
            os.mkdir(self.outpath)
        except OSError as e:
            pass

        print self.infiles
        print self.counter
        self.scores = np.genfromtxt(self.infiles[self.counter], names=True, dtype=None, skip_header=1) #, usecols=scorefile_columns_to_read)
        self.top_x_angle_list= []

    def plot_for_all_coordinates(self, tempfiles, angles_file):
        if len(tempfiles) != 0:
            self.load_template_data(tempfiles)

        self.load_PDB_angles(angles_file)
        for coord in coordinates:
            self.plot_hist_and_top_x(coord, tempfiles, angles_file)

    def load_template_data(self, tempfiles):
        tempfiles = [item.rstrip() for item in tempfiles]
        templates = TemplateScoreFiles(self.counter, tempfiles)
        templates.calculate_angles()
        self.template_array = templates.template_array

    def extract_top_x(self, coordinate, column_name=total, top_number=10):
        '''
        sorts the scorefile and extracts the 10 best decoys based on their total score.
        '''
        self.top_x_angle_list = []
        self.sorted_scores = np.sort(self.scores, order=column_name)
        for decoy in self.sorted_scores[:top_number]:
            self.top_x_angle_list.append(LHOCDecoys(decoy))

    def load_PDB_angles(self, angles_file):
        '''
        reads a file that contains only the LHOC coordinates found in crystal structures in the PDB.
        '''
        self.PDB_angles = np.genfromtxt(angles_file, delimiter=',', names=True, dtype=None)

    def plot_hist_and_top_x(self, coordinate, tempfiles, angles_file, ranked_by = total, top_number = 10):
        '''
        Plots histogram of LHOC coordinates found in crystal structures and a histogram of
        decoys that were sampled here. It also plots the LHOC coordinates of the 10 templates
        from multi grafting and the top 10 models found after H3 modeling.
        '''
        self.extract_top_x(coordinate)

        #scores for histogram
        all_scores = self.scores[coordinate]
        max_score = np.amax(all_scores)
        min_score = np.amin(all_scores)

        # histogram
        plt.clf()
        fig = plt.gcf()
        ax = plt.gca()
        fig.set_size_inches(10, 6)
        histogram = plt.hist(all_scores, bins=100, normed=True, histtype='stepfilled', linewidth=0.0, label='all decoys', color='0.6') # JJ
        PDB_hist = plt.hist(self.PDB_angles[coordinate], bins=100, normed=True, histtype='step',linewidth=1.0, label='PDB data', color='black') # JJ
        #PDB_hist = plt.hist(self.PDB_angles[coordinate], bins=100, normed=True, histtype='stepfilled', linewidth=0.0, label='PDB data', color='0.60')
        #histogram = plt.hist(all_scores, bins=100, normed=True, histtype='step', linewidth=2.0, label='all decoys', color='red')
        y, x, _ = PDB_hist
        x = (x[1:] + x[:-1])/2
        y_half_max = y.max()/2

        yy, xx, _ = histogram

        # top scoring decoys
        for decoy in self.top_x_angle_list:
            color = color_dict[decoy.template_no]
            x_val = decoy.get_coordinate(coordinate)
            label= decoy.get_coordinate(name)
            dist = np.abs(xx - x_val)
            i = np.where(dist == dist.min())
            plt.plot(x_val, yy[i], color=color, marker = 'd', label='%s: %s' %(label, x_val))
            plt.annotate(label, xy=(x_val, yy[i]),xytext=(-4,70), textcoords='offset points', rotation='vertical', size = 'small', fontweight='medium')

        # templates
        if  len(tempfiles) != 0:
            for template in self.template_array:
                color = color_dict[template['name'][-1]]
                x_val = template[coordinate]
                label = template['name']
                plt.plot(x_val, y_half_max * 0.1, color=color, marker = 'o', label = '%s: %s' %(label, x_val))


        # plot
        ymax = y.max()
        if yy.max() > y.max():
            ymax = yy.max()

        plt.ylim(0,ymax * 1.3)
        plt.xlim(x_lower[coordinate], x_upper[coordinate])

        unit = ' [$^\circ$]'
        if coordinate == D:
            unit = r'[$\AA$]'

        # for display purposes
        coord_dict = {  'VL_VH_packing_angle': 'Packing Angle',
                        'VL_VH_opposite_opening_angle':'Light Opening Angle',
                        'VL_VH_opening_angle':'Heavy Opening Angle',
                        'VL_VH_distance':'Interdomain Distance',
                     }  # Nick's new coord names
        #plt.title('%s \n%s' %(coordinate, self.infiles[self.counter]), fontsize=15)
        plt.title('%s %s' %(coord_dict[coordinate], self.name[:4]), fontsize=15)
        # plt.title('%s \n/path/to/example/plot' %coordinate, fontsize=15)
        plt.xlabel(coord_dict[coordinate] + ' ' + unit, fontsize=14)
        plt.ylabel("relative frequency", fontsize = 14)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=11,numpoints=1,handlelength=1,fancybox=True)
        plt.grid(color='0.4')

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        #print self.outpath
        plt.savefig(self.outpath+'/%s%s.pdf' %(self.name, coordinate), format='pdf', bbox_inches='tight')
