#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/loop_minimization.py
## @brief  main loop minimization protocols
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *
from modules.tools import loops as loop_tools
from modules.tools import general_tools as tools
from modules.tools import output as ouput_tools
from modules.tools import sequence as sequence_tools
from shutil import rmtree
from sys import platform
import os

class Loop_Min:        
    def __init__(self, score_object, pose):
        self.score_object = score_object
        self.pose = pose
        self.pwd = os.getcwd()
    
    def __exit__(self):
        self.score_object.score.set_weight(chainbreak, 0)
        
#classicMinLoop
    def classicMinLoop(self, rounds, loops_as_string_array, tolerance=0.1, movemap=0):
        """
        This is the classic MinMover, with dfpmin for a Loop.
        Not actually defining a loop, only regions in the movemap
        May have a problem if using Centroid....
        If Movemap is given, loops_as_string_array does not matter!!
        """

        if movemap ==0:
            movemap=MoveMap()
            movemap = loop_tools.loopMovemap(self.pose, movemap, loops_as_string_array)

        rounds=int(rounds)
        self.score_object.score.set_weight(chainbreak, 100); #Makes sure loop/domain does not break!
        print self.score_object.score(self.pose)
        min_type="dfpmin"
        minmover=MinMover(movemap, self.score_object.score, min_type, tolerance, True)
        print self.score_object.score(self.pose)
        for i in range(1, rounds+1):
            print "Rounds: "+repr(i)
            minmover.apply(self.pose)
            print self.score_object.score(self.pose)
            
        #OLD: ft.assign(ft_o)
        #OLD: p.fold_tree(ft)
        self.score_object.score.set_weight(chainbreak, 0)
        print "Minimization Complete"
        
#RelaxLoop    
    def RelaxLoop(self, rounds, loops_as_string_array, classic, movemap = 0):
        """
        Classic and Fast Relax for the loop.  Uses the Movemap from this class.
        ClassicRelax=0; FastRelax=1 for classic setting
        If Movemap is give, the loops_as_string_array does not matter...
        """

        rounds=int(rounds); classic=int(classic)

        #This Handles Vaccinity Relax.
        if movemap ==0:
            movemap = MoveMap()
            movemap = loop_tools.loopMovemap(self.pose, movemap, loops_as_string_array)

        self.score_object.score.set_weight(chainbreak, 100)
        
        if classic==0:
            Rel=ClassicRelax(self.score_object.score, movemap)
            for i in range(1, rounds+1):
                Rel.apply(self.pose)
                print self.score_object.score(self.pose)
        else:
            Rel=FastRelax(self.score_object.score)
            Rel.set_movemap(movemap)
            for i in range(1, rounds+1):
                print "Rounds: "+repr(i)
                Rel.apply(self.pose)
                print self.score_object.score(self.pose)
                
        self.score_object.score.set_weight(chainbreak, 0)
        print "Relax Complete"
    
    def relax_residue_and_neighbors(self, rounds, residue, chain, num_neighbors=1, bbonly=False):
        """
        Relaxes a residue and the residues on either side by Rosetta numbering!.
        Does not care if there is a jump between them.  This should change.
        """
        
        res = int(residue)
        res = self.pose.pdb_info().pdb2pose(chain, res)
        self.score_object.score.set_weight(chainbreak, 100)
        print self.score_object.score(self.pose)
        Rel = FastRelax(self.score_object.score)
        movemap = MoveMap()
        if not bbonly:
            for i in range(res-num_neighbors, res+num_neighbors+1):
                print i
                movemap.set_chi(i, True)
                movemap.set_bb(i, True)
        else:
            for i in range(res-num_neighbors, res+num_neighbors+1):
                movemap.set_bb(i, True)
        Rel.set_movemap(movemap)
        for i in range(1, rounds+1):
            print "Rounds: "+repr(i)
            Rel.apply(self.pose)
            print self.score_object.score(self.pose)
        self.score_object.score.set_weight(chainbreak, 0)
        print "Relax Complete"
        
    def relaxLoopBBonly(self, rounds, loops_as_string_array, classic, movemap = 0):
        """
        Relaxes a given loop.  Only BB.
        """
        movemap = MoveMap()
        movemap = loop_tools.loopBBMovemap(self.pose, movemap, loops_as_string_array)
        self.RelaxLoop(rounds, loops_as_string_array, classic, movemap)
        

    def optimizeRotLoop(self, rounds, loops_as_string_array):
        """
        Optomizes Loop Rotamers using the PackRotamersMover
        """
        rounds = int(rounds)
        print self.score_object.score(self.pose)
        packer_task=standard_packer_task(self.pose)
        packer_task.restrict_to_repacking()
        packer_task.temporarily_fix_everything()
        packer_task = loop_tools.loopChiPacker(self.pose, packer_task, loops_as_string_array)
        pack_mover=PackRotamersMover(self.score_object.score, packer_task)
        print packer_task
        for i in range(1, rounds+1):
            print "Rounds: "+repr(i)
            pack_mover.apply(self.pose)
            print self.score_object.score(self.pose)
        print "Optimization Complete"
    
#LoopBackrubRef    
    def LoopBackrubRef(self, rounds, loops_as_string_array):
        """
        Backrubs the loop using the LoopMover_Refine_Backrub Mover
        As far as I know, it does not work.  Though it would be nice.
        """
        movemap=MoveMap()
        ft = self.pose.fold_tree(); ft_o = FoldTree()
        ft_o.assign(ft)
        ft.clear()
        ft, movemap, loops_object=loop_tools.InitializeLoops(self.pose, loops_as_string_array, ft, movemap)
        print ft
        print "Fold Tree Correct? " + repr(ft.check_fold_tree())
        self.pose.fold_tree(ft)
        rounds=int(rounds)
        if self.score_object.score ==0:
            self.score_object.score = create_self.score_object.score_function_ws_patch('standard', 'self.score_object.score12')
        self.score_object.score.set_weight(chainbreak, 100); #Taking No Chances!
        print loops_object
        ref=LoopMover_Refine_Backrub(loops_object, self.score_object.score)
        print self.score_object.score(self.pose)
        for i in range(1, rounds+1):
            print "Rounds: "+repr(i)
            ref.apply(self.pose)
            print self.score_object.score(self.pose)
            #for loo in loopsLis:
                #ccd_closure = CcdLoopClosureMover(loo, movemap)
                #ccd_closure.apply(self.pose) 
        ft.assign(ft_o)
        self.pose.fold_tree(ft)
        self.score_object.score.set_weight(chainbreak, 0)
        print "LoopMover_Refine_Backrub Complete"
            
#SCWRL        
    def SCWRL(self, LoopsLis, rounds, seqFile=0):        
        print self.score_object.score(self.pose)
        pwd = os.path.split(os.path.abspath(__file__))[0]
        tempdir = pwd+"/temp"
        
        #This finds which platform you are using (Hopefully), and uses the info to run the correct version of scwrl.
        plat = tools.getOS()
        
        if plat=="error":
            print "Platform Unknown....Returning...."
        else:
            #name = p.pdb_info().name()
            if not os.path.exists(pwd+"/Scwrl/"+plat+"/Scwrl4"):
                print "SCWRL not compiled.  Please install scwrl to use scwrl features. (/SCWRL/platform)"
                return
            tempdir = pwd + "/temp"
            if not os.path.exists:
                os.mkdir(tempdir)

            if seqFile==0:
                for i in range(0, rounds):
                    print "Rounds: "+repr(i)
                    self.pose.dump_pdb(tempdir+"/temp.pdb")
                    filein = tempdir+"/temp.pdb"
                    #Gets info for loops, writes a sequence file, loads it into Scwrl.
                    fileout = tempdir + "/seqtemp.seq"
                    output_tools.saveSeqFile(self.pose, fileout, LoopsLis)
                    print pwd+"/Scwrl/"+plat+"/Scwrl4 -i "+filein+" -s "+fileout+" -0 "+"-o "+tempdir+"/new_temp.pdb"
                    os.system(pwd+"/Scwrl/"+plat+"/Scwrl4 -i "+filein+" -s "+fileout+" -0 "+"-o "+tempdir+"/new_temp.pdb")
                    x = Pose()
                    pose_from_pdb(x, (tempdir+"/new_temp.pdb"))
                    self.pose.assign(x)
                    #Removes temp files
                    rmtree(tempdir)
                    #os.mkdir(tempdir)
            else:
                for i in range(0, rounds):
                    print "Rounds: "+repr(i)
                    self.pose.dump_pdb(tempdir+"/temp.pdb")
                    filein = tempdir+"/temp.pdb"
                    os.system(pwd+"/Scwrl/"+plat+"/Scwrl4 -i "+filein+" -s "+seqFile+" -0 "+"-o "+tempdir+"/new_temp.pdb")
                    x = Pose()
                    pose_from_pdb(x, (tempdir+"/new_temp.pdb"))
                    self.pose.assign(x)
                    #Removes temp files
                    rmtree(tempdir)
                    #os.mkdir(pwd+"/temp")

        print self.score_object.score(self.pose)
