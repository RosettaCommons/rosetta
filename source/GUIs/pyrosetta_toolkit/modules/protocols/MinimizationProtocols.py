#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/MinimizationProtocols.py
## @brief  main minimization protocols, works with Regions class
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.basic.options import get_string_option
from rosetta.basic.options import get_real_option
from rosetta.basic.options import get_boolean_option

#Python Imports
from shutil import rmtree
from sys import platform
import os

#Toolkit Imports
from modules.tools import loops as loop_tools
from modules.tools import general_tools as tools
from modules.tools import output as ouput_tools
from modules.tools import sequence as sequence_tools
from ProtocolBaseClass import ProtocolBaseClass

class MinimizationProtocols(ProtocolBaseClass):        
    def __init__(self, pose, score_class, input_class, output_class):
        ProtocolBaseClass.__init__(self, pose, score_class, input_class, output_class)

        self.pwd = os.getcwd()
    
    def __exit__(self):
        self.score_class.score.set_weight(chainbreak, 0)
        
#minimize
    def minimize(self, tolerance=False, movemap=False, bb_only=False, sc_only=False):
        """
        This is the classic MinMover, with option relax:min_type defining the minimization type. Tolerance is set to .01
        Not actually defining a loop, only regions in the movemap
        May have a problem if using Centroid....
        """
        
        if not movemap:
            movemap = self.input_class.regions.get_movemap(self.pose)
        if not tolerance:
            tolerance = get_real_option('run:min_tolerance')
        
        if bb_only:
            movemap.set_chi(False)
        if sc_only:
            movemap.set_bb(False)
        
        #print "run:min_tolerance: "+repr(tolerance)
        print "-relax:min_type: "+get_string_option('relax:min_type')
        self.score_class.score.set_weight(chainbreak, 100); #Makes sure loop/domain does not break!
        minmover=MinMover(movemap, self.score_class.score, get_string_option('relax:min_type'), tolerance, True)
        self.run_protocol(minmover)
            
        self.score_class.score.set_weight(chainbreak, 0)

#relax    
    def relax(self, classic, movemap = False, bb_only=False, sc_only=False):
        """
        Classic and Fast Relax for the loop.  Uses the Movemap from this class.
        ClassicRelax=0; FastRelax=1 for classic setting
        If Movemap is give, the self.input_class.loops_as_strings does not matter...
        NOTE: RelaxProtocolBase DOES read from options system, so coordinate costraints will work.  Symmetry will not.
        """
        print "This application was created and documented by Mike Tyka, et al."
        print "For more options such as coordinate constraints, please use the options system.  Symmetry is not supported at this time."
        
        classic=int(classic)

        if not movemap:
            movemap = self.input_class.regions.get_movemap(self.pose)
            
        if bb_only:
            movemap.set_chi(False)
        if sc_only:
            movemap.set_bb(False)
            
        self.score_class.score.set_weight(chainbreak, 100)
        
        if classic==0:
            Rel=ClassicRelax(self.score_class.score, movemap)
            self.run_protocol(Rel)
        else:
            Rel=FastRelax(self.score_class.score)
            Rel.set_movemap(movemap)
            self.run_protocol(Rel)
                
        self.score_class.score.set_weight(chainbreak, 0)
    
    def relax_residue_and_neighbors(self, rounds, residue, chain, num_neighbors=1, bbonly=False):
        """
        Relaxes a residue and the residues on either side by Rosetta numbering.
        Does not care if there is a jump between them.
        Used in FullControlWindow
        """
        
        res = int(residue)
        res = self.pose.pdb_info().pdb2pose(chain, res)
        self.score_class.score.set_weight(chainbreak, 100)
        print self.score_class.score(self.pose)
        Rel = FastRelax(self.score_class.score)
        movemap = MoveMap()
        if not bbonly:
            for i in range(res-num_neighbors, res+num_neighbors+1):
                movemap.set_chi(i, True)
                movemap.set_bb(i, True)
        else:
            for i in range(res-num_neighbors, res+num_neighbors+1):
                movemap.set_bb(i, True)
        Rel.set_movemap(movemap)
        for i in range(1, rounds+1):
            print "Rounds: "+repr(i)
            Rel.apply(self.pose)
            print self.score_class.score(self.pose)
        self.score_class.score.set_weight(chainbreak, 0)
        print "Relax Complete"
        
    def relax_bb_only(self, classic, movemap = 0):
        """
        Relaxes a given loop.  Only BB.
        """
        movemap = MoveMap()
        movemap = loop_tools.loopBBMovemap(self.pose, movemap, self.input_class.loops_as_strings)
        self.relax(rounds, movemap)
        

    def optimize_rotamers(self):
        """
        Optimizes Loop Rotamers using the PackRotamersMover
        """
        
        tf = self.input_class.regions.get_basic_tf(self.pose)
        pack_mover=PackRotamersMover(self.score_class.score)
        pack_mover.task_factory(tf)
        pack_mover.score_function(self.score_class.score)
        
        self.run_protocol(pack_mover)
    
#LoopBackrubRef    
    def LoopBackrubRef(self):
        """
        Backrubs the loop using the LoopMover_Refine_Backrub Mover
        As far as I know, it does not work.  Though it would be nice.
        """
        movemap=MoveMap()
        ft = self.pose.fold_tree(); ft_o = FoldTree()
        ft_o.assign(ft)
        ft.clear()
        ft, movemap, loops_object=loop_tools.InitializeLoops(self.pose, self.input_class.loops_as_strings, ft, movemap)
        print ft
        print "Fold Tree Correct? " + repr(ft.check_fold_tree())
        self.pose.fold_tree(ft)
        rounds=int(rounds)
        if self.score_class.score ==0:
            self.score_class.score = create_self.score_class.score_function_ws_patch('standard', 'self.score_class.score12')
        self.score_class.score.set_weight(chainbreak, 100); #Taking No Chances!
        print loops_object
        ref=LoopMover_Refine_Backrub(loops_object, self.score_class.score)
        
        self.run_protocol(ref)
        ft.assign(ft_o)
        self.pose.fold_tree(ft)
        self.score_class.score.set_weight(chainbreak, 0)
            
#SCWRL        
    def SCWRL(self, seqFile=0):
        """
        This uses Scwrl4 to rebuild sidechains of only regions specified.
        Scwrl4 directory should be in the scwrl/[platform] directory.
        Does not work with multiprocessing
        """
        rounds = self.output_class.rounds.get()
        print self.score_class.score(self.pose)
        pwd = os.path.split(os.path.abspath(__file__))[0]
        tempdir = pwd+"/temp"
        
        #This finds which platform you are using (Hopefully), and uses the info to run the correct version of scwrl.
        plat = tools.getOS()
        
        if plat=="error":
            print "Platform Unknown....Returning...."
        else:
            #name = p.pdb_info().name()
            if not os.path.exists(pwd+"/scwrl/"+plat+"/Scwrl4"):
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
                    scwrl_cmd = ""
                    if self.input_class.regions:
                        fileout = tempdir + "/seqtemp.seq"
                        output_tools.saveSeqFile(self.pose, fileout, self.input_class.loops_as_strings)
                        scwrl_cmd = pwd+"/scwrl/"+plat+"/Scwrl4 -i "+filein+" -s "+fileout+" -0 "+"-o "+tempdir+"/new_temp.pdb"

                    else:
                        scwrl_cmd =  pwd+"/Scwrl/"+plat+"/Scwrl4 -i "+filein+" -o "+tempdir+"/new_temp.pdb"
                    
                    print scwrl_cmd
                    os.system(scwrl_cmd)
                        
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
                    os.system(pwd+"/scwrl/"+plat+"/Scwrl4 -i "+filein+" -s "+seqFile+" -0 "+"-o "+tempdir+"/new_temp.pdb")
                    x = Pose()
                    pose_from_pdb(x, (tempdir+"/new_temp.pdb"))
                    self.pose.assign(x)
                    #Removes temp files
                    rmtree(tempdir)
                    #os.mkdir(pwd+"/temp")

        print self.score_class.score(self.pose)
