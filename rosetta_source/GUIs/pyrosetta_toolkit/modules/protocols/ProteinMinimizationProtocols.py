#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/protein_minimization.py
## @brief  minimization protocols for whole pose
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *

#Python Imports
import os
from shutil import rmtree
from sys import platform

#Toolkit Imports
from modules.tools import general_tools as gen_tools
from ProtocolBaseClass import ProtocolBaseClass


class ProteinMinimizationProtocols(ProtocolBaseClass):
    def __init__(self, pose, score_class, input_class, output_class):
        ProtocolBaseClass.__init__(self, pose, score_class, input_class, output_class)
    
    def __exit__(self):
        self.score_class.score.set_weight(chainbreak, 0)
        
    def Relax(self, classic, movemap = 0):
        """
        Relaxes the pose using either Fast Relax or Classic Relax
        If Classic is anything other then 0, then Fast Relax occurs.
        """
        classic=int(classic)
        self.score_class.score.set_weight(chainbreak, 100)
        print self.score_class.score(self.pose)
        if classic==0:
            Rel=ClassicRelax(self.score_class.score)
            if movemap!=0:
                ClassicRelax.set_movemap(movemap)
            self.run_protocol(Rel)
        else:
            Rel=FastRelax(self.score_class.score)
            if movemap!=0:
                FastRelax.set_movemap(movemap)
            self.run_protocol(Rel)
                
        self.score_class.score.set_weight(chainbreak, 0)

    def RelaxBB(self, classic):
        movemap = MoveMap()
        movemap.set_bb(True)
        movemap.set_chi(False)
        self.Relax(classic, movemap)
        
    
    
    def classicMin(self, tolerance=0.1):
        """
        Does a classic min using Dfpmin and the classic MinMover with a tolerance of .5
        """
        
        min_type="dfpmin"
        movemap = MoveMap()
        movemap.set_bb(True)
        movemap.set_chi(True)
        minmover=MinMover(movemap, self.score_class.score, min_type, tolerance, True)
        self.run_protocol(minmover)

    def Backrub(self):
        print "Currently Implementing......."
        return p
    
    def optimizeRot(self):
        """
        This optimizes the Side Chain Rotamers of your pose.  It uses the All atom self.score_class.score functions:
        :Standard and self.score_class.score12: Perhaps later, when I learn more about OOP, a self.score_class.score can be given optionally.
        """

        packer_task=standard_packer_task(self.pose)
        packer_task.restrict_to_repacking()
        pack_mover=PackRotamersMover(self.score_class.score, packer_task)
        self.run_protocol(pack_mover)
    
    def linearizePose(self):
        for i in range(1, p.total_residue()+1):
            self.pose.set_phi(i, 180)
            self.pose.set_psi(i, 180)
        for i in range(1, p.total_residue()):
            self.pose.set_omega(i, 180)
    
    def SCWRL(self, rounds):
        """
        This uses Scwrl4 to rebuild sidechains of the entire protein.
        Scwrl4 should be in the scwrl/[platform] directory.
        """
        pwd = os.path.split(os.path.abspath(__file__))[0]
        print self.score_class.score(self.pose)
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
            for i in range(0, self.output_class.rounds):
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
        print self.score_class.score(self.pose)

