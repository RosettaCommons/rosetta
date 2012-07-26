#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/interface_and_vicinity.py
## @brief  interface and vicinity protocols.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


from rosetta import *
from modules.tools import interfaces as interface_tools
from modules.protocols.loop_minimization import *
from shutil import rmtree
from sys import platform
import os





#### INTERFACE PROTOCOLS ####

#packVicinity
class Interface_and_Vicinity(Loop_Min):
    '''
    Needs a lot of work.  Almost abandoned.
    '''
    
    def __init__(self, score_object, pose):
        Loop_Min.__init__(score_object, pose)
        
    def packVicinity(self, rounds, LisLoop, ob, type = 0, vaccinity = 5.0, onlyaround=0):
        '''
        Packs The Vacinity using either SCWRL, ROSETTA, or by adding to a design file.
        '''
        
    #Types: SCWRL, ROSETTA, DESIGN
    
        #Sets Defaults:
        if vaccinity ==0:
            vaccinity = 5.0
        if self.score_object.score ==0:
            self.score_object.score = create_score_function_ws_patch('standard', 'score12')
        if type == 0:
            type = "SCWRL"
        loops = Loops()    
        #Dumps PDB so we can read it and find out what is in the vaccinity through good old python.
        tempdir = pwd + "/temp/"
        if not os.path.exists(tempdir):
            os.mkdir(tempdir)
        self.pose.dump_pdb(tempdir+"/temp.pdb"); temp = tempdir + "/temp.pdb"
        
        #Gets Residues that are in Vicinity using python. ([res:chain]=atomic contact #)
        vacDic = interface_tools.around().getVicinity(self.pose, LisLoop, temp, vaccinity)
        
        #Sets up a new LisLoop (Takeing Res:chain, and making Res:Res:Chain so that it is understandable to loop_tools...)
        if onlyaround == 1:
            LisLoop = []; #This Clears the Listloop if onlyaround = 1
            print "Only repacking around selection...."
        for key in vacDic:
            keySP = key.split(":")
            newKey = keySP[0]+":"+keySP[0]+":"+keySP[1]
            LisLoop.append(newKey)
        
        #Kicks either SCWRL or ROSETTA to pack residues
        
        if type == "Rosetta":
            self.pose.assign(self.optimizeRotLoop(rounds, loops, LisLoop))
        if type == "SCWRL":
            #Take the VacDic, and make list of list for loops.
            self.pose.assign(self.SCWRL(LisLoop, rounds, 0))
        if type == "Design":
            #Take the VacDic, and edit the design file.
            pass
        
        #Removes temp files
        rmtree(tempdir)
        #os.mkdir(pwd+"/temp")
    
#relaxVicinity        
    def relaxVicinity(self, rounds, LisLoop, vaccinity, targetSet, vaccinitySet):
        '''
        Relaxes The Vacinity Takes options (onlyaround, fixloopbb, fixvacbb).
        '''
        
    #Types: SCWRL, ROSETTA, DESIGN
    
        #Sets Defaults:

        loops = Loops()    
        #Dumps PDB so we can read it and find out what is in the vaccinity through good old python.
        t = time.clock()
        tempdir = pwd + "/temp2/" + "_"+repr(t)
        os.mkdir(tempdir)
        self.pose.dump_pdb(tempdir+"/temp.pdb"); temp = tempdir + "/temp.pdb"
        
        #Gets Residues that are in Vicinity using python. ([res:chain]=atomic contact #)
        vacDic = interface_tools.around().getVicinity(self.pose, LisLoop, temp, vaccinity)
        
        
        #Settings = (targetSet[0], targetSet[1], vaccinitySet[0], vaccinitySet[1])
        movemap = MoveMap()
        vacLoop = []
        for key in vacDic:
            keySP = key.split(":")
            newKey = keySP[0]+":"+keySP[0]+":"+keySP[1]
            vacLoop.append(newKey)
                
                
        #First, we open up BB:
        if targetSet[0]==0:#Open Loop BB
            movemap = loop_tools.loopBBMovemap(self.pose, movemap, LisLoop)
        if vaccinitySet[0]==0:#Open Vac BB
            movemap = loop_tools.loopBBMovemap(self.pose, movemap, vacLoop)
        
        #Next, we open up Chi:
        if targetSet[1]==0:#Open Loop Chi
            movemap = loop_tools.loopChiMovemap(self.pose, movemap, LisLoop)
        if vaccinitySet[1]==0:#Open Vac Chi
            movemap = loop_tools.loopChiMovemap(self.pose, movemap, vacLoop)
            
            
        #Run Classic Min Protocol:
        self.pose.assign(self.RelaxLoop(rounds, LisLoop, 1, movemap))
        #Removes temp files
        rmtree(tempdir)
        #os.mkdir(pwd+"/temp")
        
        return p
    
#minVicinity        
    def minVicinity(self, rounds, LisLoop, vaccinity, tolerance, targetSet, vaccinitySet):
        '''
        Minimizes The Vacinity
        '''
    
        #Sets Defaults:

        loops = Loops()    
        #Dumps PDB so we can read it and find out what is in the vaccinity through good old python.
        t = time.clock()
        tempdir = pwd + "/temp2/" + "_"+repr(t)
        os.mkdir(tempdir)
        self.pose.dump_pdb(tempdir+"/temp.pdb"); temp = tempdir + "/temp.pdb"
        
        #Gets Residues that are in Vicinity using python. ([res:chain]=atomic contact #)
        vacDic = interface_tools.around().getVicinity(self.pose, LisLoop, temp, vaccinity)
        
        
        
        #Settings = (targetSet[0], targetSet[1], vaccinitySet[0], vaccinitySet[1])
        movemap = MoveMap()
        vacLoop = []
        for key in vacDic:
            keySP = key.split(":")
            newKey = keySP[0]+":"+keySP[0]+":"+keySP[1]
            vacLoop.append(newKey)
                
                
        #First, we open up BB:
        if targetSet[0]==0:#Open Loop BB
            movemap = loop_tools.loopBBMovemap(self.pose, movemap, LisLoop)
        if vaccinitySet[0]==0:#Open Vac BB
            movemap = loop_tools.loopBBMovemap(self.pose, movemap, vacLoop)
        
        #Next, we open up Chi:
        if targetSet[1]==0:#Open Loop Chi
            movemap = loop_tools.loopChiMovemap(self.pose, movemap, LisLoop)
        if vaccinitySet[1]==0:#Open Vac Chi
            movemap = loop_tools.loopChiMovemap(self.pose, movemap, vacLoop)
            
            
        #Run Classic Min Protocol:
        self.pose.assign(self.classicMinLoop(self.pose, rounds, loops, LisLoop, ob, tolerance, score, movemap))
                      

        #Removes temp files
        rmtree(tempdir)
        #os.mkdir(pwd+"/temp")
    
#packInterface        
    def packInterface(self, rounds, interface, type = 0, distance = 0, fixedRes = 0):
    #Types: SCWRL, ROSETTA, DESIGN
        
        #Grabs the interface dictionary: [::chain] = (res:chain, res:chain, etc.)
        #Represents the interfaces for each individual chain given in LisLoop.  Does not include given chain.
        #If you want that info, specify BOTH chains or ALL chains.
        #If only ONE chain is given, will pack all other chains.
        #Options are given to only pack specified interfaces (fixed - (A-B-C)).
        
        
        
        #Sets Defaults:
        if distance ==0:
            distance = 5.0
        if type == 0:
            type = "SCWRL"
        
        loops = Loops()    
        print "Distance that we should be using to calculate= "+repr(distance)
        interfaceDic = interface_tools.around().getInterface(self.pose, interface, distance)
        
        #This Creates the looplist passed onto Scwrl and Rosetta side chain packers.
        LisLoop = []
        if fixedRes!=0:
            fixedSP = fixedRes.split("-")
        
        for key in interfaceDic:
            for residues in interfaceDic[key]:
                fixed = False
                residuesSP = residues.split(":")
                if fixedRes!=0:
                    for chain in fixedSP:
                        if chain == residuesSP[1]:
                            fixed = True
                if fixed == False:
                    newRes = residuesSP[0]+":"+residuesSP[0]+":"+residuesSP[1]
                    LisLoop.append(newRes)    
        
        
        #Kicks either SCWRL or ROSETTA to pack residues
        
        if type == "Rosetta":
            self.pose.assign(self.optimizeRotLoop(self.pose, rounds, loops, LisLoop, ob, self.score_object.score))
        if type == "SCWRL":
            #Take the VacDic, and make list of list for loops.
            self.pose.assign(self.SCWRL(self.pose, LisLoop, rounds, ob, 0, self.score_object.score))
        if type == "Design":
            #Take the VacDic, and edit the design file.
            pass    
    
#relaxInterface        
    def relaxInterface(self, rounds, interface, distance = 0, fixedRes = 0, interfacepack=0):
        
        #Grabs the interface dictionary: [::chain] = (res:chain, res:chain, etc.)
        #Represents the interfaces for each individual chain given in LisLoop.  Does not include given chain.
        #If you want that info, specify BOTH chains or ALL chains.
        #If only ONE chain is given, will pack all other chains.
        #Options are given to only Relax specified interface chains (fixed - (A-B-C)).  This should be smarter.  Should be literally interface to interface.
        #However, this is for future work.
        #Currently - You can't Only Pack certain interfaces, and then relax the other.  This will be changed.
        
        loops = Loops()    
        #Dumps PDB so we can read it and find out what is in the vaccinity through good old python.
        
        interfaceDic = interface_tools.around().getInterface(self.pose, interface, distance)
        
        #This Creates the looplist passed onto Scwrl and Rosetta side chain packers.
        LisLoop = []
        if fixedRes!=0:
            fixedSP = fixedRes.split("-")
        
        for key in interfaceDic:
            for residues in interfaceDic[key]:
                fixed = False
                residuesSP = residues.split(":")
                if fixedRes!=0:
                    for chain in fixedSP:
                        if chain == residuesSP[1]:
                            fixed = True
                if fixed == False:
                    newRes = residuesSP[0]+":"+residuesSP[0]+":"+residuesSP[1]
                    LisLoop.append(newRes)    
        
        
        #Kicks either SCWRL or ROSETTA to pack residues
        
        movemap = MoveMap()
        if interfacepack ==0:
            movemap = loop_tools.loopMovemap(self.pose, movemap, LisLoop)
        else:
            packLisLoop = []
            fullLisLoop = []
            packSP = interfacepack.split(":")
            for loop in LisLoop:
                loopSP = loop.split(":")
                for chain in packSP:
                    if chain == loopSP[2]:
                        packLisLoop.append(loop)
                    else:
                        fullLisLoop.append(loop)
            movemap = loop_tools.loopMovemap(self.pose, movemap, fullLisLoop); #Creates movemap that relaxes both BB and CHI
            movemap = loop_tools.loopChiMovemap(self.pose, movemap, packLisLoop); #Creates movemap for BB only.
            
        self.pose.assign(self.RelaxLoop(rounds, LisLoop, 1, movemap))

    
    
#minInterface        
    def minInterface(self, rounds, interface, distance, fixedRes = 0, tolerance=0, interfacepack = 0):
        
        #Grabs the interface dictionary: [::chain] = (res:chain, res:chain, etc.)
        #Represents the interfaces for each individual chain given in LisLoop.  Does not include given chain.
        #If you want that info, specify BOTH chains or ALL chains.
        #If only ONE chain is given, will pack all other chains.
        #Options are given to only Relax specified interface chains (fixed - (A-B-C)).  This should be smarter.  Should be literally interface to interface.
        #However, this is for future work.
        #Currently - You can't Only Pack certain interfaces, and then relax the other.  This will be changed.
        
        loops = Loops()    
        interfaceDic = interface_tools.around().getInterface(self.pose, interface, distance)
        
        #This Creates the looplist passed onto Scwrl and Rosetta side chain packers.
        LisLoop = []
        if fixedRes!=0:
            fixedSP = fixedRes.split("-")
        
        for key in interfaceDic:
            for residues in interfaceDic[key]:
                fixed = False
                residuesSP = residues.split(":")
                if fixedRes!=0:
                    for chain in fixedSP:
                        if chain == residuesSP[1]:
                            fixed = True
                if fixed == False:
                    newRes = residuesSP[0]+":"+residuesSP[0]+":"+residuesSP[1]
                    LisLoop.append(newRes)    
        
        
        #Kicks either SCWRL or ROSETTA to pack residues
        
        movemap = MoveMap()
        if interfacepack ==0:
            movemap = loop_tools.loopMovemap(self.pose, movemap, LisLoop)
        else:
            #Keep some BB open, and some Chi open
            pass
        self.pose.assign(self.classicMinLoop(rounds, loops, LisLoop,tolerance, movemap))
