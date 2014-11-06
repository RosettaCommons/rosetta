#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/calibur.py
## @brief  Functions for running and parsing calibur.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
import sys
import os
import re

class calibur:
    def __init__(self, caliburPath=set):
        """
        This is meant to be a set of functions for running and parsing calibur
        Uses system path, or specified path in calculations.
        Calibur path should be in PATH, or specified upon construction.
        """
        
        if caliburPath ==set:
            self.caliburPath = "calibur"
        else:
            self.caliburPath = caliburPath
        
    def run_calibur(self, PDBLIST_Path, chain=False, Nter=False, Cter=False, threshold=0):
        """
        Nter, Cter are not residue numbering - pretty much it is rosetta numbering as far as I can tell...
        Need to do this better with option type thing..
        """
        run = self.caliburPath
        if chain:
            run = run+" -c "+chain
        if Nter:
            run = run+" -r "+repr(Nter)+','+repr(Cter)
        run = run+' '+PDBLIST_Path
        if threshold !=0:
            run = run+ ' -t '+repr(threshold)
        self.output = os.popen(run).readlines()
        print "Calibur Completed..."
        return
    
    def ret_centers(self):
        """
        return array of top 2 clusters, and a corresponding array of sizes.
        """
        found = False
        for line in self.output:
            print line
            if re.search("Largest 2 clusters", line):
                clusline = line.strip()
                found = True
                break
        
        #First, parse by :
        if not found:
            print "Clustering failed..."
            return
            
        sp1 = clusline.split(":")
        
        sp2 = sp1[1].split(",")
        paths = []; neighbors = []
        for pdb in sp2:
            sp3 = pdb.split("(")
            paths.append(sp3[0])
            neighbors.append(int(sp3[1].split(")")[0]))
            
        return paths, neighbors
    
    def ret_threshold(self):
        """
        Returns threshold used in analysis.
        """
        
        for line in self.output:
            if re.search("Threshold = ", line):
                sp = line.split(" = ")
        
        return sp[1].strip()
        
                
    def save_neighbors(self, resultsPath, center):
        pass

    def ret_neighbors(self, resultsPath):
        pass
    
    def ret_random_num_neighbors(self):
        pass
    def save_random_num_neighbors(self):
        pass
    
    
