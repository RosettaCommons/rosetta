#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/for_design.py
## @brief  design protocols.  
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *

class design_protocols:
    def __init__(self, score_object, pose):
        self.score_object = score_object
        self.pose = pose; #If you tie the pose with an AdvancedPyMOL object, will control observing the pose.
        #Set observer
    def packDesign(self, rounds, resFile):
        rounds=int(rounds)
        
        print self.score_object.score(self.pose)
        print self.pose
        for i in range(1, rounds+1):
            task = TaskFactory.create_packer_task(self.pose)
            task.read_resfile(resFile)
            des = PackRotamersMover(self.score_object.score, task)
            des.apply(self.pose)
            print self.score_object.score(self.pose)
        print "Pack Design Complete..."
        print self.score_object.score(self.pose)
        return self.pose
    
    def pack_residue(self, res, chain):
        '''
        Packs an individual residue.  Scores before and after.
        '''
        
        res = self.pose.pdb_info().pdb2pose(chain, int(res))
        task = TaskFactory.create_packer_task(self.pose)
        task.temporarily_fix_everything()
        task.temporarily_set_pack_residue(res, True)
        pack_mover = PackRotamersMover(self.score_object.score, task)
        print self.score_object.score(self.pose)
        pack_mover.apply(self.pose)
        print self.score_object.score(self.pose)
        
    def mutateRes(self, res, chain, new_res): 
        new_res = new_res.split(":")
        new_res = new_res[2]
        res = self.pose.pdb_info().pdb2pose(chain, int(res))
        #pList.append(Pose())
        #pList[len(pList)-1].assign(self.pose)
        print self.score_object.score(self.pose)
        print "Mutating to " + new_res
        mutate_residue(self.pose, res, new_res)
        print self.score_object.score(self.pose)
        print "Mutagenesis Complete."
        return self.pose