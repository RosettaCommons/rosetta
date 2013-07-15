#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/RegionalScoring.py
## @brief  Class/Functions for Regional Scoring, term switching, and weight/scoretype functions
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.core.scoring import ScoreType

class RegionalScoring():
    def __init__(self, pose, scorefxn):
        """
        Class for quickly scoring more complex things - Regions (loops, domains, cores?), switching energy types and (possibly) rotamer libraries on-the-fly
        Not used as a container class.
        Also for switching scoreterms (rama->rama2b, pair->fa_elec, etc.)
        Some functions may eventually be ported to C++ Rosetta if they are deemed useful enough.
        Not for General Use.  Soon to be deprecated.
        """
        
        self.pose = pose
        self.score = scorefxn
        #self.score.show(pose)
        self.scoreEnumeration = dict()
        self.scoretypes = []
        self.updateweights()
        
        self.cd_hbond_terms = ["hbond_sr_bb", "hbond_lr_bb", "hbond_sc", "hbond_bb_sc"]
        self.disulfide_terms= ["dslf_ss_dst", "dslf_cs_ang", "dslf_ss_dih", "dslf_ca_dih"]
        self.rama_terms = ["rama", "p_aa_pp"] #Could put omega on here too!
        
    def ret_pose_number(self, chain, resNum):
        """
        Returns Pose Numbering.  Takes Chain and Resnum.
        """
        return self.pose.pdb_info().pdb2pose(chain, resNum)
        
    def ret_loop_energy(self, chain, Nter, Cter):
        """
        Returns the total energy for each loop by each residue resNum.
        Weighted! ;)
        """
        
        e = 0
        poseNter = self.ret_pose_number(chain, Nter)
        poseCter = self.ret_pose_number(chain, Cter)
        for i in range(poseNter, poseCter+1):
            #print i
            for type in self.scoretypes:
                e = e+self.weights[type]*self.ret_residue_energy(i)[type]
        #print "Loop energy: "+repr(e)
        return e
    
    def ret_loop_energy_by_type(self, chain, Nter, Cter, type):
        """
        Return weighted energy by type
        Weighted!
        """
        poseNter = self.ret_pose_number(chain, Nter)
        poseCter = self.ret_pose_number(chain, Cter)  
        e=0
        for i in range(poseNter, poseCter+1):
            #print "Unweighted: " +repr(self.ret_residue_energy(i)[type])
            #print "Weighted  : "+repr(self.weights [type]*self.ret_residue_energy(i)[type])
            e = e+self.weights [type]*self.ret_residue_energy(i)[type]
        #print "Loop Energy:"+repr(e)
        return e
    
    def ret_loop_neighbor_energy(self, chain, Nter, Cter):
        """
        Returns the neighbor energy of the loop in question.
        neighbor is defined as any residue pair having a non-zero ci_2b energy.
        Divides the sum by two to account for double counting.  This has been checked in the Energies object and is true.
        Weighted!
        """
        
        emap = core.scoring.EMapVector()
        e = 0
        for i in range(Nter, Cter+1):
            poseNum = self.ret_pose_number(chain, i)
            r1 = self.pose.residue(i)
            for x in range(1, self.pose.total_residue()+1):
                for loop_residue in range(Nter, Cter+1):
                    if x==loop_residue:
                        continue
                emap.zero()
                r2 = self.pose.residue(x)
                self.score.eval_ci_2b(r1, r2, self.pose, emap)
                
                for type in self.scoretypes:
                    e = e+self.weights[type]*emap[type]/2 #Is this unweighted??? From the tutorials/manual No.
                emap.zero()
                self.score.eval_cd_2b(r1, r2, self.pose, emap)
                for type in self.scoretypes:    
                    e = e+self.weights[type]*emap[type]/2
        #print "Neighbor Energy: "+repr(e)
        return e
    
    def ret_loop_neighbor_energy_by_type(self, chain, Nter, Cter, type):
        """
        Returns the neighbor energy of the loop in question. Type is not a string!
        Note that hack_elec has hack_elec_bb_bb, hack_elec_sc_sc, and hack_elec_bb_sc.  These add up to hack_elec.  May want to see these...
        Weighted!
        """
        #print "Neighbor: "
        emap = core.scoring.EMapVector()
        e = 0
        for i in range(Nter, Cter+1):
            poseNum = self.ret_pose_number(chain, i)
            r1 = self.pose.residue(i)
            for x in range(1, self.pose.total_residue()+1):
                for loop_residue in range(Nter, Cter+1):
                    if x==loop_residue:
                        continue
                emap.zero()
                r2 = self.pose.residue(x)
                self.score.eval_ci_2b(r1, r2, self.pose, emap)
                """
                if emap[type]!=0:
                    print "ci_2b "+repr(i)+':'+repr(x)+' '+ repr(self.weights[type]*(emap[type]/2))
                    print "  unweighted "+repr(emap[type]/2)
                """
                
                e = e+self.weights[type]*(emap[type]/2)
                emap.zero()
                self.score.eval_cd_2b(r1, r2, self.pose, emap)
                """
                if emap[type]!=0:
                    print "cd_2b "+repr(i)+':'+repr(x)+' '+ repr(self.weights[type]*(emap[type]/2))
                    print "  unweighted "+repr(emap[type]/2)
                """
                
                e = e+self.weights[type]*(emap[type]/2)
        #print "Total neighbor energy: "+repr(e)
        return e
    
    def ret_total_weighted_residue_energy(self, poseNum):
        """
        Returns the total energy of the residue, weighted according to scorefunction.
        """
        self.updateweights()
        emap = self.pose.energies().residue_total_energies(poseNum)
        e = 0.0
        for type in self.scoretypes:
            e = e+self.weights[type]*emap[type]
        return e

    def ret_total_unweighted_residue_energy(self, poseNum):
        emap = self.ret_residue_energy(poseNum)
        e = 0.0
        for type in self.scoretypes:
            e = e+emap[type]
        return e
    
    def ret_residue_energy(self, poseNum):
        """
        Returns individual residue energies. By EMAP.
        UNWEIGHTED On Purpose
        """
        
        #poseNum = self.ret_pose_number(chain, resNum)
        emap = self.pose.energies().residue_total_energies(poseNum) #UNWEIGHTED!!!
        return emap
        
    def ret_n_hbonds(self):
        """
        Returns Hbonds FOR WHOLE PROTEIN
        """
        self.score(self.pose)
        set = core.scoring.hbonds.HBondSet()
        self.pose.update_residue_neighbors()
        rosetta.core.scoring.hbonds.fill_hbond_set(self.pose, False, set)
        return set.nhbonds()
        
    def ret_energy_string(self):
        """
        Recreates ScoreFunctions Show method to use return a string with output.  Used for insertion into textboxes.
        However, unsolved problem is that the spaces with textbox font size makes each line uneven.
        """
        self.score(self.pose)
        out =    "------------------------------------------------------------\n"
        out=out+ " Scores                       Weight   Raw Score Wghtd.Score\n"
        out=out+ "------------------------------------------------------------\n"
        sum_weighted=0.0;
        self.updateweights()
        for type in self.scoretypes:
            if ( self.weights[ type ] != 0.0 ):
                out=out+ ' ' +  rosetta.core.scoring.name_from_score_type(type).ljust(24,) +(" %.3f"%self.weights[ type ]).ljust(9)+"   "
                out=out+ ("%.3f"%self.pose.energies().total_energies()[ type ]).ljust(9)+ "   "
                out=out+ ("%.3f"%(self.weights[ type ] * self.pose.energies().total_energies()[ type ] )).ljust(9)+"\n"
                sum_weighted += self.weights[ type ] * self.pose.energies().total_energies()[ type ]
        
        out = out+ "---------------------------------------------------\n"
        out = out+ " Total weighted score:                    %.3f\n"%sum_weighted
        return out
    
#### Switches Don't use these - deprecated####
    
    def switch_to_fa_pair(self):
        self.score.set_weight(fa_elec, 0)
        self.score.set_weight(fa_pair, .56)
        self.updateweights()
        self.elec_id = "fa_pair"
        
    def switch_to_fa_elec(self):
        self.score.set_weight(fa_pair, 0)
        self.score.set_weight(fa_elec, .56)
        self.updateweights()
        self.elec_id = "fa_elec"
        
    def switch_to_rama(self):
        self.score.set_weight(p_aa_pp, 0)
        self.score.set_weight(rama, .54)
        self.updateweights()
        self.rama_id = "rama"
        
    def switch_to_paapp(self):
        self.score.set_weight(rama, 0)
        self.score.set_weight(p_aa_pp, .54)
        self.updateweights()
        self.rama_id = "p_aa_pp"
    
    def enable_rama2b(self):
        self.score.set_weight(rama2b, .2)
        self.updateweights()
        
    def switch_on_orbitals(self):
        #Don't know the correct weights, or anything!
        self.hbond_id = "orbitals"
        pass
    
    def switch_off_orbitals(self):
        pass
    
    def updateweights(self):
        self.score(self.pose)
        self.weights = self.pose.energies().weights()
        self.scoretypes = []
        for i in range(1, rosetta.core.scoring.end_of_score_type_enumeration+1):
            ii = rosetta.core.scoring.ScoreType(i)
            if self.weights[ii] != 0: self.scoretypes.append(ii)
