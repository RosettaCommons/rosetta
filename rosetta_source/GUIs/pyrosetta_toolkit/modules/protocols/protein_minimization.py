#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/protein_minimization.py
## @brief  minimization protocols for whole pose
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *
from modules.tools import general_tools as gen_tools
import os
from shutil import rmtree
from sys import platform

class Protein_Min:
    def __init__(self, score_object, pose):
        self.score_object = score_object
        self.pose = pose
        
    def Relax(self, rounds, classic, movemap = 0):
        '''
        Relaxes the pose using either Fast Relax or Classic Relax
        If Classic is anything other then 0, then Fast Relax occurs.
        This is useful as a Tk Checkbutton. (uses self.score_object.score 12 and standard)
        '''
        rounds=int(rounds); classic=int(classic)
        self.score_object.score.set_weight(chainbreak, 100)
        print self.score_object.score(self.pose)
        if classic==0:
            Rel=ClassicRelax(self.score_object.score)
            if movemap!=0:
                ClassicRelax.set_movemap(movemap)
            for i in range(1, rounds+1):
                Rel.apply(self.pose)
                print self.self.score_object.score(self.pose)
        else:
            Rel=FastRelax(self.score_object.score)
            if movemap!=0:
                FastRelax.set_movemap(movemap)
            for i in range(1, rounds+1):
                Rel.apply(self.pose)
                print self.score_object.score(self.pose)
        print "Relax Complete"

    
    def RelaxBB(self, rounds, classic):
        movemap = MoveMap()
        movemap.set_bb(True)
        movemap.set_chi(True)
        self.Relax(rounds, classic, movemap)
        
    
    
    def classicMin(self, rounds, tolerance=0.1):
        '''
        Does a classic min using Dfpmin and the classic MinMover with a tolerance of .5
        Later will add choices if needs be
        '''
        

        print self.score_object.score(self.pose)
        min_type="dfpmin"
        movemap = MoveMap()
        movemap.set_bb(True)
        movemap.set_chi(True)
        minmover=MinMover(movemap, self.score_object.score, min_type, tolerance, True)
        
        for i in range(1, rounds+1):
            minmover.apply(self.pose)
            print self.score_object.score(self.pose)
        print "Minimization Complete"

    def Backrub(self, rounds, ob):
        print "Currently Implementing......."
        return p
    
    def optimizeRot(self, rounds):
        '''
        This optimizes the Side Chain Rotamers of your pose.  It uses the All atom self.score_object.score functions:
        :Standard and self.score_object.score12: Perhaps later, when I learn more about OOP, a self.score_object.score can be given optionally.
        '''

        print self.score_object.score(self.pose)
        packer_task=standard_packer_task(self.pose)
        packer_task.restrict_to_repacking()
        pack_mover=PackRotamersMover(self.score_object.score, packer_task)
        for i in range(1, rounds+1):
            pack_mover.apply(self.pose)
            print self.score_object.score(self.pose)
        print "Optimization Complete"
    
    def linearizePose(self):
        for i in range(1, p.total_residue()+1):
            self.pose.set_phi(i, 180)
            self.pose.set_psi(i, 180)
        for i in range(1, p.total_residue()):
            self.pose.set_omega(i, 180)
    
    def SCWRL(self, rounds):
        
        '''
        This uses Scwrl4 to rebuild sidechains of the entire protein.  Ob is not int, rounds is int.
        '''
        pwd = os.path.split(os.path.abspath(__file__))[0]
        print self.score_object.score(self.pose)
        tempdir = pwd + "/temp/"
        if not os.path.exists:
            os.mkdir(tempdir)

        #This finds which platform you are using (Hopefully), and uses the info to run the correct version of scwrl.
        plat = gen_tools.getOS()
        if plat=="error":
            print "Platform Unknown....Returning...."
        else:
            if not os.path.exists(pwd+"/Scwrl/"+plat+"/Scwrl4"):
                print "SCWRL not compiled.  Please install scwrl to use scwrl features. (/SCWRL/platform)"
                return
            for i in range(0, rounds):
                #This may need to be changed to a call function in the future.
                self.pose.dump_pdb(tempdir+"/temp.pdb")
                filein = tempdir +"/temp.pdb"
                ScwrlInput = pwd+"/Scwrl/"+plat+"/Scwrl4 -i "+filein+" -o "+tempdir+"/new_temp.pdb"
                print ScwrlInput
                os.system(ScwrlInput)
                x = Pose()
                pose_from_pdb(x, (tempdir+"/new_temp.pdb"))
                self.pose.assign(x)
                #Removes temp files
                rmtree(tempdir)
                #os.mkdir(tempdir)
        print self.score_object.score(self.pose)

