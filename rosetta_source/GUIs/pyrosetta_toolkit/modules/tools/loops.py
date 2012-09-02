#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/tools/loops.py
## @brief  general loop functions for the toolkit
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


from rosetta import *
import math
import random
import os
import general_tools
import re

#### LOOP CLOSURE ####
def closeLoop_ccd( rounds, loo, movemap):
    """
    Closes Loop using CCD. Loo is a Rosetta Loops object.
    """
    
    ccd_closure = CcdLoopClosureMover(loo, movemap)
    for i in range(1, rounds+1):
        ccd_closure.apply(self.pose)
        self.score_object.score(self.pose)
    print "CCD complete"
    
#### LOOP INITIALIZATION ####

def return_residue_array( p, loop_string):
    """
    Returns an array in Rosetta numbering from the loop_string (start:end:chain)
    Works with individual residues, loops, termini, and chains.
    Use this to interact with movemap, packer, etc and loops_as_strings from GUI.
    """
    
    residue_array = []
    loop_stringSP = loop_string.split(":")
    if (loop_stringSP[0] == "" and loop_stringSP[1]==""):
        #print "Getting Sequence for Chain:"
        for i in range(1, (p.total_residue())):
            x = (p.pdb_info().pose2pdb(i)).split()
            if x[1]==loop_stringSP[2]:
                residue_array.append(i)
                #print "added residue " +x[0]+ " to the array"
    elif loop_stringSP[0]=="":
        #print "Getting N-terminus of chain "+loop_stringSP[2] + " From" + loop_stringSP[1] + " to end of chain..."
        for i in range(1, (p.total_residue())):
            x = (p.pdb_info().pose2pdb(i)).split()
            if (x[1] == loop_stringSP[2] and i<= p.pdb_info().pdb2pose(loop_stringSP[2], int(loop_stringSP[1]))):
                residue_array.append(i)
                #print "added residue " +x[0] +" To the array"
    elif loop_stringSP[1]=="":
        #print "Getting C-terminus of chain "+loop_stringSP[2]+ " From start of chain to" + loop_stringSP[0]
        for i in range(1, (p.total_residue())):
            x = (p.pdb_info().pose2pdb(i)).split()
            if (x[1] == loop_stringSP[2] and i>= p.pdb_info().pdb2pose(loop_stringSP[2], int(loop_stringSP[0]))):
                residue_array.append(i)
                #print "added residue " + x[0] + " To the array"
    else:
        #print "Getting Loop Region"
        start=p.pdb_info().pdb2pose(loop_stringSP[2], int(loop_stringSP[0]));
        end=p.pdb_info().pdb2pose(loop_stringSP[2], int(loop_stringSP[1]));
        for z in range(start, end+1):
            residue_array.append(z)
    return residue_array

def InitializeLoop( start, end, p, chain):
    """
    This Initializes the loop, uses the Pose numbering, returns a loop with a cutpoint in the middle
    Sort of pointless, I know.  Really don't know why it's still here.
    """
    
    #Initializes loop, and makes a cut in middle#
    #start = p.pdb_info().pdb2pose(chain, start); end = p.pdb_info().pdb2pose(chain,end)
    
    loo=Loop(start,end)
    loo.auto_choose_cutpoint(p)
    return loo


def loopMovemap( p, movemap, loops_as_strings):
    """
    Opens Both Chi and BB for LisLoop.  Does so for Chains, Termini, Loops, and Individual Residues
    """
    for loop_string in loops_as_strings:
        residue_array = return_residue_array(p, loop_string)
        for residue in residue_array:
            movemap.set_chi(residue, True)
            movemap.set_bb(residue, True)
    print repr(movemap)
    return movemap


def loopChiMovemap( p, movemap, loops_as_strings):
    """
    Opens Chi for LisLoop.  Works with Chains, Loops, Termini, Individuals.
    """
    
    for loop_string in loops_as_strings:
        residue_array = return_residue_array(p, loop_string)
        for residue in residue_array:
            movemap.set_chi(i, True)
    print repr(movemap)
    return movemap
    
def loopBBMovemap( p, movemap, loops_as_strings):
    """
    Opens BB up for LisLoops.  Works with Chains, Termini, Loops, and Individual Residues
    """
    for loop_string in loops_as_strings:
        residue_array = return_residue_array(p, loop_string)
        for residue in residue_array:
            movemap.set_bb(i, True)
    print repr(movemap)
    return movemap
    
def loopChiPacker( p, pack, loops_as_strings):
    """
    Returns a Packer_task open with the chi specified.  Takes in a Packer Task.  Works with chains, loops, termini, individuals.
    """
    for loop_string in loops_as_strings:
        residue_array = return_residue_array(p, loop_string)
        for i in residue_array:
            if re.search("ProteinFull", p.residue(i).name())==None:
                pack.temporarily_set_pack_residue(i, True)
    return pack


        
def setLoopBreak( p, start, end, chain, cut):
    start = p.pdb_info().pdb2pose(chain, int(start)); end = p.pdb_info().pdb2pose(chain,int(end))
    loo=Loop(start,end, int(cut))
    self.InitializeLoopMovemap(start, end)
    set_single_loop_fold_tree(p, loo)
    print "Break Setup"
    print p
    return p

def loopArea( p, loops_as_strings, rosetta=1):
    """
    Takes in the raw Loops List, Returns a new loop list of list.  rosetta=1; list of resdidue_arrays for each loop., or regular of res:chain (rosetta=0)
    Used for creating a sequence file for scwrl.
    """
    
    newList = []
    for loop_string in loops_as_strings:
        residue_array = return_residue_array(p, loop_string)
        newList.append(residue_array)
  
    if rosetta == 1:
        return newList
    if rosetta == 0:
        fulllist = [] #FullList - List of residues res:chain
        for lists in newList:
            for res in lists:
                x = p.pdb_info().pose2pdb(res)
                x = x.split()
                new = x[0]+":"+x[1]
                fulllist.append(new)
        return fulllist

def return_rosetta_Loop(p, loop_string):
    """
    Determines start and end of piece, returns a Rosetta Loop object.
    """
    residue_array = return_residue_array(p, loop_string)
    start = residue_array[0]; end = residue_array[-1]
    cut = start+((end-start)/2)
    rosetta_loop = Loop(start, end, cut)
    return rosetta_loop

def InitializeLoops( p, loops_as_strings, ft=0, movemap=0):
    """
    DEPRECATED? Returns loop information giving a list filled with loops, a movemap, and a new foldtree to apply
    """
    #Movemap for each type works...
    begin = 1
    loopsLis = []
    i =1
    for x in loops_as_strings:
        LoopFull = x.split(":")
        if LoopFull[0]=="" and LoopFull[1] =="":
            print "Modeling chain "+LoopFull[2]
            
            chain = LoopFull[2]
            start = p.pdb_info().pdb2pose(chain, 1)
        elif LoopFull[0]=="":
            print "Modeling N-terminus of chain "+LoopFull[2]
            print "Not Yet Implemented"
            return
        elif LoopFull[1]=="":
            print "Modeling C-terminus of chain "+LoopFull[2]
            print "Not Yet Implemented..."
            return
        else:
            start=p.pdb_info().pdb2pose(LoopFull[2], int(LoopFull[0]));
            end=p.pdb_info().pdb2pose(LoopFull[2], int(LoopFull[1]));
            cut = start+((end-start)/2)
            loo = Loop(start, end, cut)
            loopsLis.append(loo)
            #ft.new_jump(start-2, end+2, cut)
            if ft!=0:
                if i ==1:
                    ft.add_edge(begin, start-2, -1)
                else:
                    ft.add_edge(begin+2, start-2, -1)
                ft.add_edge(start-2, cut, -1)
                ft.add_edge(start-2, end+2, i)
                ft.add_edge(end+2, cut+1, -1)
                begin = end
                i+=1
                ft.add_edge(begin+2, int(p.total_residue()), -1)
        if movemap !=0:
            movemap = self.loopMovemap(p, movemap, loops_as_strings)
    if ft ==0 and movemap==0:
        return loopsLis
    elif ft ==0 and movemap!=0:
        return(movemap, loopsLis)
    elif ft!=0 and movemap ==0:
        return(ft, loopsLis)
    else:
        return(ft, movemap, loopsLis)
#### Random Loop Tools ####


def RandLoop( p, LisLoops, ob):
    """
    Randomizes the loop.  Completely.  It can get pretty ugly....
    
    """
    for loop in LisLoops:
        loopSP = loop.split(":")
        start=int(loopSP[0]); end= int(loopSP[1]); chain = loopSP[2]; ob = int(ob)
        start = p.pdb_info().pdb2pose(loopSP[2], start); end = p.pdb_info().pdb2pose(loopSP[2],end)
        if ob ==1:
            self.obs.add_observer(p)
        start=int(start); end=int(end)
        loo=self.InitializeLoop(start, end, p, chain)
        ft = p.fold_tree(); ft_o = FoldTree()
        ft_o.assign(ft)
        set_single_loop_fold_tree(p, loo)
        
        for i in range(start, end+1):
            print "Completely Randomizing Loop"
            print i
            phi_angle=p.phi(i); phi_perturb=(phi_angle-90+random.random()*180)
            psi_angle=p.psi(i); psi_perturb=(psi_angle-90+random.random()*180)
            p.set_phi(i, phi_perturb); p.set_psi(i, psi_perturb)
        movemap= self.InitializeLoopMovemap(start, end)
        closer = Close_Loop_Algorithms()
        p.assign(closeLoop_ccd(p, 10, loo, movemap))
        ft.assign(ft_o)
        p.fold_tree(ft)
    return p

def linearizeLoops( p, LisLoop, ob):
    """
    Linearizes the LisLoop; Only works for real loops!
    """
    for loop in LisLoop:
        loopSp = loop.split(":")
        start = int(loopSp[0]); end = int(loopSp[1]); chain = loopSp[2]
        start = p.pdb_info().pdb2pose(chain, start); end = p.pdb_info().pdb2pose(chain,end)
        loo=self.InitializeLoop(start, end, p, chain)
        set_single_loop_fold_tree(p, loo)
        for i in range(start, end+1):
            print "Completely Randomizing Loop"
            print i
            p.set_phi(i, 180); p.set_psi(i, 180)
        movemap= self.InitializeLoopMovemap(start, end)
        p.assign(closeLoop_ccd(p, 10, loo, movemap))
    return p

#### Loop Editing ####

def delLoop( p, start, end, chain):
    """
    Deletes the Loop from the pdb.
    Soon will output the original residue numbering.
    """
    
    self.obs.add_observer(p)
    start = int(start); end = int(end); print start; print end;
    start = p.pdb_info().pdb2pose(chain, start); end = p.pdb_info().pdb2pose(chain,end)
    length=end-start
    print length
    dump_pdb(p, "temp_del.pdb")
    pstart = Pose(); pstart.assign(p)
    
    
    for i in range(1, length+2):
        p.delete_polymer_residue(start)
    return p

def delCDRLoops( p):
    #L1
    start = p.pdb_info().pdb2pose(L, 24); end = p.pdb_info().pdb2pose(L,42)
    length=end-start
    print length
    for i in range(1, length+1):
        p.delete_plolymer_residue(start)


def insertNewLoop( p, seq, start, end, chain):
    """
    Warning:  Does not work well.
    Start-end signifies where between the loop you want to put something.
    """
    
    p2=Pose(); #end=int(end)
    seq = "A"+seq+"A"
    make_pose_from_sequence(p2, seq, "fa_standard")
    print p2
    length=p2.total_residue()
    TrueLength = length-2
    norm=180
    res_num=1
    while res_num<(p2.total_residue()+1):
        p2.set_phi(res_num, norm)
        p2.set_psi(res_num, norm)
        p2.set_omega(res_num, norm)
        res_num +=1
        Ores=1
    
    end = p.pdb_info().pdb2pose(chain, (int(end)))
    cut=p.pdb_info().pdb2pose(chain, int(start))
    for i in range(1, length-1):
        print i
        x = length-i; print x
        res=p2.residue(x); print p2.residue(x)
        p.prepend_polymer_residue_before_seqpos(res, end, True)
        loo = Loop(cut-1, cut+i, cut)
        movemap = initLoops().InitializeLoopMovemap(cut-1, cut+i)
        set_single_loop_fold_tree(p, loo)
        ccd_closure = CcdLoopClosureMover(loo, movemap)
        ccd_closure.apply(p)
        #p.set_phi(end, norm); p.set_psi(end, norm); p.set_omega(end, norm)
    
    
    
    #This tries to set the residues inserted to a good starting position before CCD...
    
    dL = tools.general_tools.getDist(p, int(start), int(end), "CA", "CA")
    print "Distance between cuts:  " +repr(dL)
    VarSres=round(dL/3.75)
    VarStartRes = (TrueLength-VarSres)/2
    print "Residues needed to span length:  " + repr(VarSres)
    pwd = os.getcwd(); name = pwd + "/"+"temp.pdb"
    dump_pdb(p, name)
    pose_from_pdb(p, name)
    os.remove(name)
    print "Sequence Inserted....."
    print p
    print "New Sequence Numbering: Residues " +repr(cut+1) +" To " + repr(cut+TrueLength) + " Chain A"
    return p