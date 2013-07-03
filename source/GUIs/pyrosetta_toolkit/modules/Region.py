#!/usr/bin/python
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/Region.py
## @brief  This class will eventually replace the loop_string and loops_as_strings region selection of the GUI.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.core.pack.task.operation import *
from rosetta.protocols.toolbox.task_operations import *

#Python Imports
import re

class Region:
    """
    This is a PDB numbering based region class.  Its primary purpose is to make dealing with regions/selections through the GUI and code more intuitive.
    It will eventually replace the loop_string in loops_as_strings in the main GUIinput class.
    Region is inclusive - Though, Lets think about this.
    Might eventually make it's way into C++ Rosetta in a similar form.
    """
    
    def __init__(self, chain, start=None, end=None):
        """
        What you give determines the region_type.
        Give only a chain, and the region is a chain.
        Give chain and start, and the region is a C-Terminal tail
        Give chain and end,   and the region is a N-Terminal tail
        These designations are mainly for cutpoint + Foldtree related things.
        **Please use self.check_if_region_exists(pose) before using the region.
        """
        
        assert type(chain) == str
        if type(start)== str:
            start = int(start)
            
        if type(end) ==  str:
            end = int(end)
            
        self.start = start
        self.end = end
        self.chain = chain
        
        if not self.start and not self.end:
            self.region = ":"+":"+chain.upper()
        elif not self.start:
            self.region = ":"+repr(end)+":"+chain.upper()
        elif not self.end:
            self.region = repr(start)+":"+":"+chain.upper()
        else:
            self.region = repr(start)+":"+repr(end)+":"+chain.upper()
    
    def __str__(self):
        return self.region
    
    def region_exists(self, pose):
        """
        Checks if the region exists for the given Pose.
        Returns True or False
        """
        if not rosetta.core.pose.has_chain(self.chain, pose):
            print "Chain not found in PDB"
            return False
        
        try:
            self.get_rosetta_start(pose)
            self.get_rosetta_end(pose)
        except PyRosettaException:
            return False
        
        return True
    
    def set_Loop_for_region(self, pose, cutpoint):
        #if not self.get_region_type()=='loop':
            #return
        self.loop = Loop(self.get_rosetta_start(pose), self.get_rosetta_end(pose), cutpoint)
        
    def get_loop(self):
        if not self.loop: return
        return self.loop
    
    def get_movemap(self, pose):
        movemap = MoveMap()
        for i in range(self.get_rosetta_start(pose), self.get_rosetta_end(pose)+1):
            movemap.set_bb(i, True)
            movemap.set_chi(i, True)
        return movemap
    
    def get_region(self):
        """
        Gets the region string.
        """
        return self.region
    
    def get_region_string(self):
        """
        Redundant.  Gets region string.
        """
        return self.region
    
    def get_region_string_with_all_residues(self, pose):
        """
        Returns the loops string with full residue numbers instead of blank ::
        """
        s = repr(self.get_pdb_start(pose))+":"+repr(self.get_pdb_end(pose))+":"+self.get_chain()
        return s
    
    def get_region_type(self):
        if not self.start and not self.end:
            return "chain"
        elif not self.start:
            return 'nter'
        elif not self.end:
            return 'cter'
        elif self.start == self.end:
            return 'residue'
        else:
            return 'loop'
        
    #Main Getters - These will be blank depending on region type.
    def get_start(self):
        return self.start
    def get_end(self):
        return self.end
    def get_chain(self):
        return self.chain

    #Rosetta Getters
    def get_rosetta_start(self, pose):
        
        if ((self.get_region_type()=='nter') or (self.get_region_type()=='chain')):
            #First residue of the chain - Would be useful to have pdb_info() be able to return this.
            for resnum in range(1, pose.total_residue()+1):
                chain = pose.pdb_info().chain(resnum)
                if chain == self.chain:
                    return resnum
                else:
                    continue
        else:
            return pose.pdb_info().pdb2pose(self.chain, self.start)
    

    def get_rosetta_end(self, pose):

        if (self.get_region_type()=='cter') or (self.get_region_type()=='chain'):
            #Last residue of the chain - Would be useful to have pdb_info() be able to return this.
            resnum = pose.total_residue()
            while resnum !=1:
                chain = pose.pdb_info().chain(resnum)
                #print repr(resnum); print chain; print "Should Be: "+self.chain
                if chain == self.chain:
                    return resnum
                else:
                    resnum-=1
                    continue
        else:
            return pose.pdb_info().pdb2pose(self.chain, self.end)
            
    def get_pdb_start(self, pose):
        
        return pose.pdb_info().number(self.get_rosetta_start(pose))
        
    def get_pdb_end(self, pose):
        
        return pose.pdb_info().number(self.get_rosetta_end(pose))
         
    #Useful functions
    def get_length(self, pose):
        return self.get_rosetta_end(pose)-self.get_rosetta_start(pose)+1
    
    def get_sequence(self, pose):
        if self.get_region_type()=='chain':
            c = rosetta.core.pose.get_chain_id_from_chain(self.chain, pose)
            return pose.chain_sequence(c)
        else:
            sequence = pose.sequence()
            return sequence[self.get_rosetta_start(pose)-1:self.get_rosetta_end(pose)]
            
class Regions:
    """
    This class is analogous to the Loops class in Rosetta.
    Iterable.
    Like the Region class, may eventually be ported to C++.
    """
    
    def __init__(self):
        self.regions = []
        self.current = 0
    
    def __str__(self):
        loops = ""
        for region in self.regions:
            loops = loops+" "+str(region)
        return loops
    
    def __iter__(self):
        return self
    
    def __len__(self):
        return len(self.regions)
        
    def __nonzero__(self):
        if len(self.regions)>0:
            return True
        else:
            return False
    

    
    def next(self):
        if self.current >= len(self.regions):
            self.current = 0
            raise StopIteration
        else:
            region = self.regions[self.current]
            self.current+=1
            return region
            
    def iterkeys(self):
        return self
    
    def add_region(self, region_object):
        self.regions.append(region_object)
    
    def remove_region(self, region_string):
        regions = self.regions
        for region in regions:
            if region.get_region()==region_string:
                self.regions.remove(region)
                break
            
    def get_regions(self):
        return self.regions
    
    def get_num_regions(self):
        return len(self.regions)
        
    def get_regions_of_type(self, type):
        regions = []
        for region in self.regions:
            if region.get_region_type()==type:
                regions.append(region)
        return regions
    
    def get_regions_of_types(self, array_of_types):
        regions = []
        if len(array_of_types)>1:
            for type in array_of_types:
                for region in self.regions:
                    if region.get_region_type()==type:
                        regions.append()
            return regions
        else:
            return self.get_regions_of_type(array_of_types[0])
            
    def get_regions_of_chain(self, chain):
        regions = []
        for region in self.regions:
            if region.get_chain()==chain:
                regions.append(region)
        return regions
    
    def get_cter_regions(self):
        return self.get_regions_of_type('cter')
    def get_nter_regions(self):
        return self.get_regions_of_type('nter')
    def get_residue_regions(self):
        return self.get_regions_of_type('residue')
    def get_chain_regions(self):
        return self.get_regions_of_type('chain')
    def get_loop_regions(self):
        return self.get_regions_of_type('loop')
        
    def get_movemap(self, pose, include_only_regions=False):
        """
        Create a movemap where bb and sc are on from the regions in the region object.
        """
        regions = self.regions
        movemap = MoveMap()
        
        if not self.regions:
            for i in range(1, pose.total_residue()+1):
                movemap.set_bb(i, True)
                movemap.set_chi(i, True)
            return movemap
        
        if include_only_regions:
            regions = self.get_regions_of_types(include_only_regions)
        
        for region in regions:
            #print region.get_region()
            start = region.get_rosetta_start(pose)
            end = region.get_rosetta_end(pose)
            for i in range(start, end+1):
                movemap.set_bb(i, True)
                movemap.set_chi(i, True)
        return movemap
    
    def get_packer_task(self, pose, include_only_regions = False):
        """
        Create a packer task for repacking.
        """
        regions = self.regions
        packer_task=standard_packer_task(pose)
        packer_task.restrict_to_repacking()
        
        
        if self.regions:    
            packer_task.temporarily_fix_everything()
        else:
            return packer_task
            
        if include_only_regions:
            regions = self.get_regions_of_types(include_only_regions)
        
        for region in regions:
            
            
            start = region.get_rosetta_start(pose)
            end = region.get_rosetta_end(pose)
            print "Rosetta start: "+repr(start)
            print "Rosetta end: "+repr(end)
            for i in range(start, end+1):
                packer_task.temporarily_set_pack_residue(i, True)
        print packer_task
        return packer_task
    
    def get_basic_tf(self, pose):
        """
        Returns a basic tf - no design, no neighbors.  Only regions set.
        """
        regions = self.regions
        tf = TaskFactory()
        tf.push_back(InitializeFromCommandline())
        
        tf.push_back(RestrictToRepacking())
        mmop = RestrictToMoveMapChiOperation(self.get_movemap(pose))
        tf.push_back(mmop)
        return tf
    
    #Loop/Region integration:
    def get_Loops_from_regions(self):
        """
        Must have used set_Loop_for_region for the Region object to use.
        """
        loops = Loops()
        for region in self.regions:
            if region.loop:
                loops.add_loop(region.get_loop())
        return loops
    