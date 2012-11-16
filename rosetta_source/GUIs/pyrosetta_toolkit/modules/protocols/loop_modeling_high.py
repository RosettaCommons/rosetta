#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/loop_modeling_high.py
## @brief  High res loop modeling protocols. 
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *

class highRes_Loop_Modeling:
    """
    Simple class for handling HighRes loop modeling...
    """
    
    def __init__(self, score_object, pose):
        self.score_object = score_object
        self.pose = pose
    
#default_CCD        
    def default_CCD(self, rounds, loops_as_string_array):
        """
        Uses LoopMover_Refine_CCD as High Res Mover
        """
        movemap=MoveMap()
        ft = self.posefold_tree(); ft_o = FoldTree()
        ft_o.assign(ft)
        ft.clear()
        #ft.simple_tree(self.posetotal_residue())
        print len(loops_as_string_array)
        ft, movemap, loops_object=tools.loops.initLoops().InitializeLoops(p, loops_as_string_array, ft, movemap)
        print "Fold Tree Correct? " + repr(ft.check_fold_tree())
        self.posefold_tree(ft)
        score = create_score_function_ws_patch('standard', 'score12')
        ob=int(ob); rounds=int(rounds)            
        if ob == 1:
            self.obs.add_observer(self.pose)
        loop_refine = LoopMover_Refine_CCD(loops_object)
        loop_refine.set_scorefxn(score)
        print self.score_object.score(self.pose)
        for i in range(1, rounds+1):    
            loop_refine.apply(self.pose)
        ft.assign(ft_o)
        self.posefold_tree(ft)
        print self.score_object.score(self.pose)
    
#default_KIC        
    def default_KIC(self, rounds, loops):
        """
        Uses LoopMover_Refine_KIC as High Res Mover
        """
        ob=int(ob); rounds=int(rounds)
        if ob ==1:
            self.obs.add_observer(self.pose)
        movemap=MoveMap()
        ft = self.posefold_tree(); ft_o = FoldTree()
        ft_o.assign(ft)
        ft.clear()
        #ft.simple_tree(self.posetotal_residue())
        print len(loops_as_string_array)
        ft, movemap, loops_object=tools.loops.initLoops().InitializeLoops(p, loops_as_string_array, ft, movemap)

        print "Fold Tree Correct? " + repr(ft.check_fold_tree())
        self.posefold_tree(ft)
        loop_refine = LoopMover_Refine_KIC(loops_object)
        loop_refine.set_scorefxn(score)
        for i in range(1, rounds+1):    
            loop_refine.apply(self.pose)
        ft.assign(ft_o)
        self.posefold_tree(ft)
        return p
