#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/loop_modeling_low.py
## @brief  Low res loop modeling protocols
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *

class lowRes_Loop_Modeling:
     def __init__(self, score_object, pose):
        self.score_object = score_object
        if not self.score_object.centroid_score:
            self.score_object.centroid_score = create_score_function('cen_std')
        
#default_CCD        
     def default_CCD(self, rounds, loops_as_string_array, ob, fragset, fraglength):
         '''
         LoopMover_Perturb_CCD for Low Resolution Loop Modeling
         REQUIRES FRAGSET!!
         '''
         
         movemap=MoveMap()
         ft = self.pose.fold_tree(); ft_o = FoldTree()
         ft_o.assign(ft)
         ft.clear()
         #ft.simple_tree(self.pose.total_residue())
         print len(loops_as_string_array)
         ft, movemap, loopsLis=tools.loops.initLoops().InitializeLoops(p, loops_as_string_array, ft, movemap)
         loops = Loops()
         for loo in loopsLis:
             loops.add_loop(loo)
         print "Fold Tree Correct? " + repr(ft.check_fold_tree())
         self.pose.fold_tree(ft)
         #Unfortunately, still requires a damn fragset
         Frag = ConstantLengthFragSet(fraglength, fragset)
         x = LoopMover_Perturb_CCD(loops, self.score_object.centroid_score, Frag)
         for i in range(1, rounds+1):
             print self.score_object.centroid_score(self.pose); x.apply(self.pose); print self.score_object.centroid_score(self.pose)
         ft.assign(ft_o)
         self.pose.fold_tree(ft)

         
#default_KIC            
     def default_KIC(self, p, rounds, loops_as_string_array, ob, fragset, fraglength):
         '''
         LoopMover_Perturb_KIC for Low Resolution Loop Modeling
         REQUIRES FRAGSET!!
         '''
         movemap=MoveMap()
         ft = self.pose.fold_tree(); ft_o = FoldTree()
         ft_o.assign(ft)
         ft.clear()
         #ft.simple_tree(self.pose.total_residue())
         print len(loops_as_string_array)
         ft, movemap, loopsLis=tools.loops.initLoops().InitializeLoops(p, loops_as_string_array, ft, movemap)
         loops = Loops()
         for loo in loopsLis:
             loops.add_loop(loo)
         print "Fold Tree Correct? " + repr(ft.check_fold_tree())
         self.pose.fold_tree(ft)
         #Unfortunately, still requires a damn fragset
         Frag = ConstantLengthFragSet(fraglength, fragset)
         x = LoopMover_Perturb_KIC(loops, self.score_object.centroid_score)
         x.add_fragments(Frag)
         for i in range(1, rounds+1):
             print self.score_object.centroid_score(self.pose); x.apply(self.pose); print self.score_object.centroid_score(self.pose)
         
         ft.assign(ft_o)
         self.pose.fold_tree(ft)

#fragInsert_Anneal            
     def fragInsert_Anneal(self, rounds, loops_as_string_array, ob, fragset, fraglength):
         '''
         Anneals a given fragset and a given fraglength using KT from 2 to .8.
         Inner Loop is equal to 50 rounds.  kT is scaled by 1/(inner*outer) for each inner cycle down to .8
         Final kT will become an option later down the road....
         Outer Loop is chosen by user, as well as the scoring function used.
         Applies CcdLoopClosureMover after each frag insert....
         '''

         print fraglength
         print fragset
         fraglength = int(fraglength)
         movemap=MoveMap()
         ft = self.pose.fold_tree(); ft_o = FoldTree()
         ft_o.assign(ft)
         ft.clear()
         #ft.simple_tree(self.pose.total_residue())
         print len(loops_as_string_array)
         ft, movemap, loopsLis=tools.loops.initLoops().InitializeLoops(p, loops_as_string_array, ft, movemap)
         loops = Loops()

         print "Fold Tree Correct? " + repr(ft.check_fold_tree())
         self.pose.fold_tree(ft)
         self.score_object.centroid_score.set_weight(chainbreak, 100)
         frag = ConstantLengthFragSet(fraglength, fragset)
         
         #Handles Vaccinity Repacking.  Could do it each time a new loop is created and accepted.  However, this would take a lot of time.
         #Will implement it on cluster, or depending on how long it would take. (Which would still NOT be ideal....)(Should change each time a fragment is inserted, before scoring and accepting.)
         if cutoff != 0:
             #Get Vaccinity!
             vacDic = tools.input.load_vicinity(p, loops_as_string_array, cutoff)
             LisLoop = []; #Makes a new listloop
             print "Only repacking around selection...."
             for key in vacDic:
                 keySP = key.split(":")
                 newKey = keySP[0]+":"+keySP[0]+":"+keySP[1]
                 LisLoop.append(newKey)
             movemap = tools.loops.initLoops().loopChiMovemap(p, movemap, LisLoop)
            
         
         
         fragMove = ClassicFragmentMover(frag,movemap)
         
         
         #Note - Not the best implementation of ccd_closure, but its the only think I could think of...
         inner_cycles = 50
         init_temp = 2.0
         final_temp = 0.8
         gamma = math.pow((final_temp/init_temp),(1.0/(rounds*inner_cycles)))
         kT = init_temp
         mc = MonteCarlo(p, self.score_object.centroid_score, kT)
         self.score_object.centroid_score.set_weight(chainbreak, 100)
         print self.score_object.centroid_score(self.pose)
         for i in range(1,rounds+1):
             print "Low Anneal Round"
             print rounds
             mc.recover_low(self.pose)
             self.score_object.centroid_score(self.pose)
             for j in range(1, inner_cycles+1):
                 kT = kT * gamma
                 print kT
                 mc.set_temperature(kT)
                 fragMove.apply(self.pose)
                 for loo in loopsLis:
                     ccd_closure = CcdLoopClosureMover(loo, movemap)
                     ccd_closure.apply(self.pose)
                     
                 mc.boltzmann(self.pose)
                 #print self.score_object.centroid_score(self.pose)
         #mc.recover_low(self.pose)
         print self.score_object.centroid_score(self.pose)
         mc.show_counters()
                 
         ft.assign(ft_o)
         self.pose.fold_tree(ft)
             
         print "fragInsert_Anneal Complete..."    
